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
import os
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

        os.rename(
            ase_calculator.label + ".com", 
            os.path.join(
                ase_calculator.scratch,
                ase_calculator.label + ".com"
            ))

        os.rename(
            ase_calculator.label + ".ase", 
            os.path.join(
                ase_calculator.scratch,
                ase_calculator.label + ".ase"
            ))

    def fix_io_file(self, calc=None):
        """
        A method that removes the `left` and `right` text from a log, ase, and
        com files and turns it back into a smiles structure
        """
        if calc:
            old_log_file = calc.label + ".log"
            old_log_path = os.path.join(calc.scratch, old_log_file)
            if os.path.exists(old_log_path):
                new_log_path = old_log_path.replace(
                    "left", "(").replace("right", ")")
                os.rename(old_log_path, new_log_path)

            old_ase_file = calc.label + ".ase"
            old_ase_path = os.path.join(calc.scratch, old_ase_file)
            if os.path.exists(old_ase_path):
                new_ase_path = old_ase_path.replace(
                    "left", "(").replace("right", ")")
                os.rename(old_ase_path, new_ase_path)

            old_com_file = calc.label + ".com"
            old_com_path = os.path.join(calc.scratch, old_com_file)
            if os.path.exists(old_com_path):
                new_com_path = old_com_path.replace(
                    "left", "(").replace("right", ")")
                os.rename(old_com_path, new_com_path)

        else:
            logging.info("No calculator object provided... not doing anything")

    def submit(self, calculator, partition):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        scratch = calculator.scratch
        file_path = os.path.join(scratch, calculator.label)
        complete_file_path = file_path.replace(
            "left", "(").replace("right", ")")
        command = calculator.command.split()[0]

        if os.path.exists(complete_file_path + ".log"):
            logging.info("It looks like this has already been attempted")
            try:
                self.read_log(complete_file_path + ".log")
            except BaseException:
                logging.info("This file failed... running it again")
                os.environ["COMMAND"] = "g16"
                os.environ["FILE_PATH"] = file_path
                subprocess.call(
                    "sbatch --exclude=c5003 --job-name={0} --output={0}.slurm.log --error={0}.slurm.log -p {1} -N 1 -n 20 --mem=100000 submit.sh".format(
                        calculator.label, partition), shell=True)
        elif os.path.exists(file_path + ".log"):
            # a workaround for if a file already exists
            logging.info("It looks like this has already been attempted")
            try:
                self.read_log(file_path + ".log")

            except BaseException:
                logging.info("This file failed... running it again")
                os.environ["COMMAND"] = "g16"
                os.environ["FILE_PATH"] = file_path
                subprocess.call(
                    "sbatch --exclude=c5003 --job-name={0} --output={0}.slurm.log --error={0}.slurm.log -p {1} -N 1 -n 20 --mem=100000 submit.sh".format(
                        calculator.label, partition), shell=True)
        else:
            os.environ["COMMAND"] = "g16"
            os.environ["FILE_PATH"] = file_path
            subprocess.call(
                "sbatch --exclude=c5003 --job-name={0} --output={0}.slurm.log --error={0}.slurm.log -p {1} -N 1 -n 20 --mem=100000 submit.sh".format(
                    calculator.label,
                    partition),
                shell=True)

    def check_complete(self, ase_calculator=None):
        """
        A method to determine if a job is still running
        """
        command = "squeue -n {}".format(ase_calculator.label)
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
        calc = calculator.get_species_calc(conformer=conformer)
        logging.info("Submitting {} conformer".format(calc.label))
        self.write_input(calc, conformer)
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

    def calculate_species(self, species=None, calculator=None):
        """
        A method that submits all conformers for a species and checks if the jobs are complete
        """

        logging.info("Submitting {}".format(species))
        # This doesn't work rn because RDKit is acting up on discovery
        species.generate_conformers(calculator=Hotbit())
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
                        # Skipping in the case of complete
                        continue
                    if check_complete(ase_calculator=calc):
                        if not first_try[calc.label]:
                            logging.info("{} is complete!".format(calc.label))
                            try:  # reading in results after the first attempt
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                # conf.update_coords() ### Broken cuz of RDKit
                                conf.energy = conf.ase_molecule.get_potential_energy()
                                complete[calc.label] = True
                                first_try[calc.label] = True
                            except BaseException:
                                logging.info(
                                    "{} failed, trying it again...".format(
                                        calc.label))
                                submit_conformer(
                                    conformer=conf, calculator=calculator)
                                first_try[calc.label] = True
                        else:  # look into trying something
                            logging.info("{} second attempt is complete")
                            try:
                                conf.ase_molecule = self.read_log(
                                    os.path.join(calc.scratch, calc.label + ".log"))
                                # conf.update_coords() ### Broken cuz of RDKit
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
        self.write_input(calc, conformer)
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
            if self.verify_rotor(conformer=conformer, ase_calculator=calc):
                verified[torsion.index] = True
            else:
                logging.info("Something went wrong with {}".format(torsion))


    def verify_rotor(self, conformer, ase_calculator):

        if (self.check_rotor_slope(ase_calculator) and self.check_rotor_continuous(ase_calculator)):
            return True
        else:
            return False

    def check_rotor_slope(self, ase_calculator):

        file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

        parser = ccread(file_name)


    def check_rotor_continuous(self, ase_calculator):
        file_name = os.path.join(ase_calculator.scratch, ase_calculator.label + ".log")

        parser = ccread(file_name)



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
        self.write_input(calc, ts)
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
        self.write_input(calc, ts)
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
        self.write_input(calc, ts)
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
        self.write_input(calc, ts)
        self.submit(calc, "west")

    def submit_ts(self, ts=None, calculator=None, vibrational_analysis=True):
        """
        A method to perform the entire ts optimization for a transition state
        """

        assert isinstance(ts, TS), "The object provided is not a TS object"

        # For shell
        self.submit_shell(ts=ts, calculator=calculator)
        calc = calculator.get_shell_calc(ts)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))
        try:
            ts.ase_molecule = self.read_log(os.path.join(
                calc.scratch,
                calc.label + ".log"
            ))
            ts.update_coords()
            logging.info("{} was successful!".format(calc.label))
            self.fix_io_file(calc)
        except BaseException:
            try:
                ts.ase_molecule = self.read_log(
                    os.path.join(
                        calc.scratch,
                        calc.label.replace(
                            "left",
                            "(").replace(
                            "right",
                            ")") +
                        ".log"))
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
                self.fix_io_file(calc)
            except BaseException:
                logging.info("FAILED: {} calculation".format(calc.label))
                self.fix_io_file(calc)
                return False

        # For center
        self.submit_center(ts=ts, calculator=calculator)
        calc = calculator.get_center_calc(ts)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))

        try:
            ts.ase_molecule = self.read_log(os.path.join(
                calc.scratch,
                calc.label + ".log"
            ))
            ts.update_coords()
            logging.info("{} was successful!".format(calc.label))
            self.fix_io_file(calc)
        except BaseException:
            try:
                ts.ase_molecule = self.read_log(
                    os.path.join(
                        calc.scratch,
                        calc.label.replace(
                            "left",
                            "(").replace(
                            "right",
                            ")") +
                        ".log"))
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
                self.fix_io_file(calc)
            except BaseException:
                logging.info("FAILED: {} calculation".format(calc.label))
                self.fix_io_file(calc)
                return False

        # For overall
        self.submit_overall(ts=ts, calculator=calculator)
        calc = calculator.get_overall_calc(ts)
        while not self.check_complete(ase_calculator=calc):
            time.sleep(30)
        logging.info("{} is complete!".format(calc.label))

        try:
            ts.ase_molecule = self.read_log(os.path.join(
                calc.scratch,
                calc.label + ".log"
            ))
            ts.update_coords()
            logging.info("{} was successful!".format(calc.label))
            self.fix_io_file(calc)
        except BaseException:
            try:
                ts.ase_molecule = self.read_log(
                    os.path.join(
                        calc.scratch,
                        calc.label.replace(
                            "left",
                            "(").replace(
                            "right",
                            ")") +
                        ".log"))
                ts.update_coords()
                logging.info("{} was successful!".format(calc.label))
                self.fix_io_file(calc)
            except BaseException:
                logging.info("FAILED: {} calculation".format(calc.label))
                self.fix_io_file(calc)
                return False

        # For validation
        if vibrational_analysis:  # If we're running vibrational analysis
            from autotst.calculators.vibrational_analysis import VibrationalAnalysis
            vib = VibrationalAnalysis(ts=ts, scratch=calculator.scratch)
            result = vib.validate_ts()

            if not result:
                logging.info(
                    "Vibrational analysis not definitive. Running IRC instead")
                self.submit_irc(ts=ts, calculator=calculator)
                calc = calculator.get_irc_calc(ts)
                while not self.check_complete(ase_calculator=calc):
                    time.sleep(30)
                logging.info("{} is complete!".format(calc.label))

                result = calculator.validate_irc(calc=calc)
                if not result:
                    logging.info("Cannot validate via IRC... Job failed")
                    return False
                return True

        else:
            self.submit_irc(ts=ts, calculator=calculator)
            calc = calculator.get_irc_calc(ts)
            while not self.check_complete(ase_calculator=calc):
                time.sleep(30)
            logging.info("{} is complete!".format(calc.label))

            result = calculator.validate_irc(calc=calc)

            if not result:
                logging.info("Cannot validate via IRC... Job failed")
                return False
            return True

#    def calculate_ts(self, ts=None, calculator=None):
