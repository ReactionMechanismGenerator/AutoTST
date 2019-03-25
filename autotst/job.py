import rdkit.DistanceGeometry
import rdkit.Chem.rdDistGeom
import rdkit
from autotst.calculators.statmech import StatMech
from autotst.calculators.vibrational_analysis import VibrationalAnalysis, percent_change
from autotst.calculators.gaussian import Gaussian
from ase.calculators.gaussian import Gaussian as ASEGaussian
from ase.atoms import Atom, Atoms
from autotst.calculators.calculator import Calculator
from ase.io.gaussian import read_gaussian_out
from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
from autotst.molecule import Molecule
from autotst.reaction import Reaction, TS
import os, time
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.reaction import Reaction as RMGReaction, ReactionError
from rmgpy.species import Species as RMGSpecies
from rmgpy.molecule import Molecule as RMGMolecule
import rmgpy
import ase
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Chem import AllChem
from rdkit import Chem
import numpy as np
import logging
import subprocess
from shutil import move
import cclib
import time
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


# To perform TS search


class Job():
    """
    A class to deal with the input and output of calculations
    """

    def __init__(
            self,
            reaction=None,
            dft_calculator=None,
            conformer_calculator=None):

        self.reaction = reaction
        if isinstance(reaction, Reaction):
            self.label = reaction.label

        self.dft_calculator = dft_calculator
        self.conformer_calculator = conformer_calculator

    def __repr__(self):
        return "< Job '{}'>".format(self.calculator, self.label)

    def read_log(self, file_path=None):
        """
        A helper method that allows one to easily parse log files
        """

        symbol_dict = {
            8: "O",
            6: "C",
            7: "N",
            1: "H",
        }

        atoms = []
        from cclib.io import ccread
        parser = ccread(file_path)

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):

            atoms.append(Atom(symbol=symbol_dict[atom_num], position=coords))

        return Atoms(atoms)

    def write_input(self, conformer, ase_calculator):
        """
        A helper method that will write an input file and move it to the correct scratch directory
        """
        ase_calculator.write_input(conformer.ase_molecule)

        move(
            ase_calculator.label + ".com", 
            os.path.join(
                ase_calculator.scratch,
                ase_calculator.label + ".com"
            ))

        move(
            ase_calculator.label + ".ase", 
            os.path.join(
                ase_calculator.scratch,
                ase_calculator.label + ".ase"
            ))

    def submit(self, calculator, partition):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        scratch = calculator.scratch
        file_path = os.path.join(scratch, calculator.label)
        complete_file_path = file_path.replace(
            "left", "(").replace("right", ")")
        command = calculator.command.split()[0]
        
        if isinstance(calculator, ase.calculators.gaussian.Gaussian):
            os.environ["COMMAND"] = "g16"
        else:
            logging.info("Assuming this is a Gaussian Calculator...")
            os.environ["COMMAND"] = "g16"

        os.environ["FILE_PATH"] = os.path.join(scratch, file_path)

        if os.path.exists(complete_file_path + ".log"):
            logging.info("It looks like this has already been attempted")
            try:
                self.read_log(complete_file_path + ".log")
            except BaseException:
                logging.info("This file failed... running it again")
                subprocess.call(
                    """sbatch --exclude=c5003 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 submit.sh""".format(
                        calculator.label, partition), shell=True)
        elif os.path.exists(file_path + ".log"):
            # a workaround for if a file already exists
            logging.info("It looks like this has already been attempted")
            try:
                self.read_log(file_path + ".log")

            except BaseException:
                logging.info("This file failed... running it again")
                subprocess.call(
                    """sbatch --exclude=c5003 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 submit.sh""".format(
                        calculator.label, partition), shell=True)
        else:
            subprocess.call(
                """sbatch --exclude=c5003 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 submit.sh""".format(
                    calculator.label, partition), shell=True)

    def check_complete(self, ase_calculator=None):
        """
        A method to determine if a job is still running
        """
        command = """squeue -n "{}" """.format(ase_calculator.label)
        output = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE).communicate()[0]
        if len(output.split("\n")) <= 2:
            return True
        else:
            return False

    # This portion is specific for Species / Conformers
    def submit_conformer(self, conformer=None, calculator=None):
        """
        A method that uses submit to run the geometry optimization of a single conformer
        """
        calc = calculator.get_species_calc(conformer=conformer,
                                        mem=calculator.mem,
                                        nprocshared=calculator.nprocshared,
                                        scratch=calculator.scratch,
                                        method=calculator.method,
                                        basis=calculator.basis
                                        )
        logging.info("Submitting {} conformer".format(calc.label))
        logging.info("{} has the following torsions".format(calc.label))
        for torsion in conformer.torsions:
            logging.info(torsion)
            logging.info("\t-{}".format(torsion.mask))

            
        self.write_input(conformer, calc)
        self.submit(calc, "general")
        return calc

    def submit_species(self, species=None, calculator=None):
        """
        A method that uses submit to run the geometry optimization of an entire species
        """
        calcs = []
        for smiles, confs in species.conformers.items():
            for conf in confs:
                calc = self.submit_conformer(
                    conformer=conf, calculator=calculator)
                calcs.append(calc)
        return calcs

    def calculate_species(self, species=None, conformer_calculator=None, calculator=None):
        """
        A method that submits all conformers for a species and checks if the jobs are complete
        """

        logging.info("Submitting {}".format(species))
        # This doesn't work rn because RDKit is acting up on discovery
        if conformer_calculator:
            species.generate_conformers(calculator=conformer_calculator)
        
        calcs = self.submit_species(species=species, calculator=calculator)
        complete = {}
        first_try = {}
        for calc in calcs:
            complete[calc.label] = False
            first_try[calc.label] = False

        while not all(complete.values()):
            for smiles, confs in species.conformers.items():
                for conf in confs:
                    calc = calculator.get_species_calc(conf)
                    if complete[calc.label]:
                        continue
                    if self.check_complete(ase_calculator=calc):
                        if not first_try[calc.label]:
                            logging.info("{} is complete!".format(calc.label))
                            try:  # reading in results after the first attempt
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                conf.update_coords() ### Broken cuz of RDKit
                                conf.energy = conf.ase_molecule.get_potential_energy()
                                complete[calc.label] = True
                                first_try[calc.label] = True
                            except BaseException:
                                logging.info(
                                    "{} failed, trying it again...".format(
                                        calc.label))
                                self.submit_conformer(
                                    conformer=conf, calculator=calculator)
                                first_try[calc.label] = True
                        else:  # look into trying something
                            logging.info("{} second attempt is complete")
                            try:
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                conf.update_coords() ### Broken cuz of RDKit
                                conf.energy = conf.ase_molecule.get_potential_energy()
                                complete[calc.label] = True
                            except BaseException:
                                logging.info(
                                    "{} failed a second time...".format(
                                        calc.label))
                                complete[calc.label] = True
    
    def run_rotor(self, conformer, torsion, steps, step_size):
        """
        a method to run a hindered rotor calculation for a single rotor
        """
        calc = calculator.get_rotor_calc(conformer=conformer, torsion=torsion, steps=steps, step_size=step_size)
        logging.info("Running hindered rotor calculation for {}".format(torsion))
        self.write_input(conformer, calc)
        self.submit(calc, "general")

    def run_rotors(self, conformer, steps, step_size):
        """
        A method to run hindered rotor scans for all torsions in a conformer
        """

        for torsion in conformer.torsions:
            self.run_rotor(conformer=conformer, torsion=torsion, steps=steps, step_size=step_size)

    def calculate_rotors(self, conformer, calculator, steps, step_size):
        """
        A method to submit and verify all hindered rotor scans for a conformer
        """
        self.run_rotors(conformer=conformer, steps=steps, step_size=step_size)
        complete = {}
        verified = {}
        for torsion in conformer.torsions:
            complete[torsion.index] = False
            verified[torsion.index] = False

        while not all(complete.values()):
            for torsion in conformer.torsions:
                calc = calculator.get_rotor_calc(conformer=conformer, torsion=torsion)
                if self.check_complete(calc):
                    complete[torsion.index] = True
        
        for torsion in conformer.torsions:
            calc = calculator.get_rotor_calc(conformer=conformer, torsion=torsion)
            if self.verify_rotor(conformer=conformer, ase_calculator=calc, steps=steps, step_size=step_size):
                verified[torsion.index] = True
            else:
                logging.info("Something went wrong with {}".format(torsion))


    def verify_rotor(self, steps, step_size, ase_calculator=None, file_name=None, parser=None):
        
        assert (ase_calculator is not None) or (file_name is not None) or (parser is not None)
        
        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

            parser = cclib.io.ccread(file_name)
        
        continuous = self.check_rotor_continuous(steps, step_size, ase_calculator=ase_calculator, file_name=file_name, parser=parser)
        
        good_slope = self.check_rotor_slope(steps, step_size, ase_calculator=ase_calculator, file_name=file_name, parser=parser)
        
        [lowest_conf, energy, atomnos, atomcoords] = self.check_rotor_lowest_conf(ase_calculator=ase_calculator, file_name=file_name, parser=parser)
        
        opt_count_check = self.check_rotor_opts(steps, ase_calculator=ase_calculator, file_name=file_name, parser=parser)
        
        return [lowest_conf, continuous, good_slope, opt_count_check] 

    def check_rotor_opts(self, steps, ase_calculator=None, file_name=None, parser=None):
        
        assert (ase_calculator is not None) or (file_name is not None) or (parser is not None)
        
        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

            parser = cclib.io.ccread(file_name)

        opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]
        
        n_opts_check = (steps + 1)==len(opt_SCFEnergies)
        
        return n_opts_check
    
    
    
    def check_rotor_slope(self, steps, step_size, ase_calculator=None, file_name=None, parser=None, tol=0.35):
        
        assert (ase_calculator is not None) or (file_name is not None) or (parser is not None)
        
        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

            parser = cclib.io.ccread(file_name)

        opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]
        
        max_energy = max(opt_SCFEnergies)
        min_energy = min(opt_SCFEnergies)
        
        max_slope = (max_energy - min_energy) / step_size
        slope_tol = tol*max_slope
        
        for i, energy in enumerate(opt_SCFEnergies):
            prev_energy = opt_SCFEnergies[i-1]
            slope = np.absolute((energy-prev_energy)/float(step_size))
            if slope > slope_tol:
                return False
        
        return True
        
        
    def check_rotor_continuous(self, steps, step_size, ase_calculator=None, file_name=None, parser=None, tol=0.05):
        
        assert isinstance(step_size, float)
        
        assert (ase_calculator is not None) or (file_name is not None) or (parser is not None)
        
        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

            assert os.path.isfile(file_name)
            parser = cclib.io.ccread(file_name)

            
        opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]
        
        max_energy = max(opt_SCFEnergies)
        min_energy = min(opt_SCFEnergies)
        energy_tol = np.absolute(tol*(max_energy - min_energy))
        
        checked = [None for angle in range(0,360)]

        continuous = True

        for step, energy in enumerate(opt_SCFEnergies):
            abs_theta = int(step*step_size)
            theta = abs_theta%360

            mismatch = False

            if checked[theta] is None:
                checked[theta] = energy

            else:
                checked_energy = checked[theta]

                abs_diff = np.absolute(energy - checked_energy)

                if abs_diff > energy_tol:
                    mismatch = True
                    continuous = False
                    return False

                #print abs_theta, theta, checked_energy, energy, mismatch, abs_diff, energy_tol

        return continuous
    
    
    
    def check_rotor_lowest_conf(self, ase_calculator=None, file_name=None, parser=None, tol=0.03):
        
        assert (ase_calculator is not None) or (file_name is not None) or (parser is not None)
        
        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

            assert os.path.isfile(file_name)
            parser = cclib.io.ccread(file_name)

        opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]
        
        max_energy = max(opt_SCFEnergies)
        min_energy = min(opt_SCFEnergies)
        energy_tol = tol*(max_energy - min_energy)
        
        
        first_is_lowest = True #Therefore...
        min_idx = 0
        min_energy = opt_SCFEnergies[min_idx]
        
        for i, energy in enumerate(opt_SCFEnergies):
            if min_energy - energy > energy_tol:
                min_energy = energy
                min_idx = i
        
        if min_idx != 0:
            first_is_lowest = False
        
        min_opt_idx = opt_indices[min_idx]
        
        atomnos = parser.atomnos
        atomcoords = parser.atomcoords[min_opt_idx]
        
        return [first_is_lowest, min_energy, atomnos, atomcoords]


    def submit_shell(self, ts=None, calculator=None):
        """
        A method to submit a calculation for the reaction shell
        """
        scratch = calculator.scratch
        calc = calculator.get_shell_calc(ts,
                                         mem=calculator.mem,
                                         nprocshared=calculator.nprocshared,
                                         scratch=calculator.scratch,
                                         method=calculator.method,
                                         basis=calculator.basis
                                         )
        logging.info("Submitting {} TS".format(calc.label))
        self.write_input(ts, calc)
        self.submit(calc, "general")

    def submit_center(self, ts=None, calculator=None):
        """
        A method to submit a calculation for the reaction center
        """
        scratch = calculator.scratch
        calc = calculator.get_center_calc(ts,
                                          mem=calculator.mem,
                                          nprocshared=calculator.nprocshared,
                                          scratch=calculator.scratch,
                                          method=calculator.method,
                                          basis=calculator.basis
                                          )
        logging.info("Submitting {} TS".format(calc.label))
        self.write_input(ts, calc)
        self.submit(calc, "general")

    def submit_overall(self, ts=None, calculator=None):
        """
        A method to submit a calculation for the entire ts
        """
        scratch = calculator.scratch
        calc = calculator.get_overall_calc(ts,
                                           mem=calculator.mem,
                                           nprocshared=calculator.nprocshared,
                                           scratch=calculator.scratch,
                                           method=calculator.method,
                                           basis=calculator.basis
                                           )
        logging.info("Submitting {} TS".format(calc.label))
        self.write_input(ts, calc)
        self.submit(calc, "general")

    def submit_irc(self, ts=None, calculator=None):
        scratch = calculator.scratch
        calc = calculator.get_irc_calc(ts,
                                       mem=calculator.mem,
                                       nprocshared=calculator.nprocshared,
                                       scratch=calculator.scratch,
                                       method=calculator.method,
                                       basis=calculator.basis
                                       )
        logging.info("Submitting {} TS".format(calc.label))
        self.write_input(ts, calc)
        self.submit(calc, "west")

    def submit_ts(self, ts=None, calculator=None, vibrational_analysis=True):
        """
        A method to perform the entire ts optimization for a transition state
        """

        assert isinstance(ts, TS), "The object provided is not a TS object"

        # For shell
        self.submit_shell(ts=ts, calculator=calculator)
        calc = calculator.get_shell_calc(ts, scratch=calculator.scratch)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))
        try:
            file_name = os.path.join(
                calc.scratch,
                calc.label + ".log"
            )
            if all(calculator.verify_output_file(file_name)):
                ts.ase_molecule = self.read_log(file_name)
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
            else:
                logging.info("FAILED: {} calculation".format(calc.label))
                return False
        except BaseException:
            logging.info("FAILED: {} calculation".format(calc.label))
            return False

        # For center
        self.submit_center(ts=ts, calculator=calculator)
        calc = calculator.get_center_calc(ts, scratch=calculator.scratch)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))
        try:
            file_name = os.path.join(
                calc.scratch,
                calc.label + ".log"
            )
            if all(calculator.verify_output_file(file_name)):
                ts.ase_molecule = self.read_log(file_name)
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
            else:
                logging.info("FAILED: {} calculation".format(calc.label))
                return False
        except BaseException:
            logging.info("FAILED: {} calculation".format(calc.label))
            return False

        # For overall
        self.submit_overall(ts=ts, calculator=calculator)
        calc = calculator.get_overall_calc(ts, scratch=calculator.scratch)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))
        try:
            file_name = os.path.join(
                calc.scratch,
                calc.label + ".log"
            )
            if all(calculator.verify_output_file(file_name)):
                ts.ase_molecule = self.read_log(file_name)
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
            else:
                logging.info("FAILED: {} calculation".format(calc.label))
                return False
        except BaseException:
            logging.info("FAILED: {} calculation".format(calc.label))
            return False

        # For validation
        if vibrational_analysis:  # If we're running vibrational analysis
            from autotst.calculators.vibrational_analysis import VibrationalAnalysis, percent_change
            vib = VibrationalAnalysis(ts=ts, scratch=calculator.scratch)
            result = vib.validate_ts()

            if not result:
                logging.info(
                    "Vibrational analysis not definitive. Running IRC instead")
                self.submit_irc(ts=ts, calculator=calculator)
                calc = calculator.get_irc_calc(ts=ts, scratch=calculator.scratch)
                while not self.check_complete(ase_calculator=calc):
                    time.sleep(30)
                logging.info("{} is complete!".format(calc.label))
                logging.info("The scratch directory is {}".format(calc.scratch))
                result = calculator.validate_irc(calc=calc)
                if not result:
                    logging.info("Cannot validate via IRC... Job failed")
                    return False
                return True

        else:
            self.submit_irc(ts=ts, calculator=calculator)
            calc = calculator.get_irc_calc(ts, scratch=calculator.scratch)
            while not self.check_complete(ase_calculator=calc):
                time.sleep(30)
            logging.info("{} is complete!".format(calc.label))

            result = calculator.validate_irc(calc=calc)

            if not result:
                logging.info("Cannot validate via IRC... Job failed")
                return False
            return True

    def submit_reaction(self, reaction=None, calculator=None, vibrational_analysis=True):
        """
        A method that uses submit to run the geometry optimization of an entire species
        """
        processes = []
        for smiles, transitionstates in reaction.ts.items():
            for ts in transitionstates:
                p = Process(target=self.submit_ts, args=(ts, calculator, vibrational_analysis))
                p.start()
                processes.append(p)
        return processes

    def calculate_reaction(self, reaction=None, calculator=None):

        logging.info("Submitting {}".format(reaction))
        # This doesn't work rn because RDKit is acting up on discovery
        if conformer_calculator:
            reaction.generate_conformers(calculator=conformer_calculator)
        
        calcs = self.submit_species(species=species, calculator=calculator)
        complete = {}
        first_try = {}
        for calc in calcs:
            complete[calc.label] = False
            first_try[calc.label] = False

        while not all(complete.values()):
            for smiles, confs in species.conformers.items():
                for conf in confs:
                    calc = calculator.get_species_calc(conf)
                    if complete[calc.label]:
                        continue
                    if self.check_complete(ase_calculator=calc):
                        if not first_try[calc.label]:
                            logging.info("{} is complete!".format(calc.label))
                            try:  # reading in results after the first attempt
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                conf.update_coords() ### Broken cuz of RDKit
                                conf.energy = conf.ase_molecule.get_potential_energy()
                                complete[calc.label] = True
                                first_try[calc.label] = True
                            except BaseException:
                                logging.info(
                                    "{} failed, trying it again...".format(
                                        calc.label))
                                self.submit_conformer(
                                    conformer=conf, calculator=calculator)
                                first_try[calc.label] = True
                        else:  # look into trying something
                            logging.info("{} second attempt is complete")
                            try:
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                conf.update_coords() ### Broken cuz of RDKit
                                conf.energy = conf.ase_molecule.get_potential_energy()
                                complete[calc.label] = True
                            except BaseException:
                                logging.info(
                                    "{} failed a second time...".format(
                                        calc.label))
                                complete[calc.label] = True


