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
from autotst.species import Species, Conformer
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
from cclib.io import ccread
import time
import yaml
import multiprocessing
from multiprocessing import Process, Manager
from shutil import copyfile
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
        elif reaction is None:
            self.label = None
        else:
            assert False, "Reaction provided was not a reaction object"

        self.dft_calculator = dft_calculator
        self.conformer_calculator = conformer_calculator

    def __repr__(self):
        return "< Job '{}'>".format(self.dft_calculator, self.label)

    def read_log(self, file_path=None):
        """
        A helper method that allows one to easily parse log files
        """
        symbol_dict = {
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }
        atoms = []

        parser = ccread(file_path)

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
            atoms.append(Atom(symbol=symbol_dict[atom_num], position=coords))

        return Atoms(atoms)

    def write_input(self, conformer, ase_calculator):
        """
        A helper method that will write an input file and move it to the correct scratch directory
        """

        ase_calculator.write_input(conformer.ase_molecule)
        if not os.path.exists(ase_calculator.scratch):
            os.makedirs(ase_calculator.scratch)
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


    def submit_conformer(self, conformer, ase_calculator, partition="general"):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        assert conformer, "Please provide a conformer to submit a job"

        self.write_input(conformer, ase_calculator)

        label = conformer.smiles + "_{}".format(conformer.index) 
        scratch = ase_calculator.scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16" #only using gaussian for now 
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info("It appears that this job has already been run, not running it a second time.")
                
        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 $AUTOTST/autotst/submit.sh""".format(
                    label, partition), shell=True)

        return label
    
    def submit_transitionstate(self, transitionstate, ase_calculator, partition="general"):
        """
        A methods to submit a job for a TS object based on a single calculator
        """
        assert transitionstate, "Please provide a transitionstate to submit a job"

        self.write_input(transitionstate, ase_calculator)
        label = ase_calculator.label
        scratch = ase_calculator.scratch
        print scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16" #only using gaussian for now 
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info("It appears that this job has already been run")
                
        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 $AUTOTST/autotst/submit.sh""".format(
                    label, partition), shell=True)

        return label

    def calculate_transitionstate(self, transitionstate, calculator):
        """
        A method to perform the partial optimizations for a transitionstate and arrive
        at a final geometry. Returns True if we arrived at a final geometry, returns false
        if there is an error along the way.
        """
        #######################
        ##### For shells  #####
        #######################
        direction = transitionstate.direction

        calc = calculator.get_shell_calc(transitionstate,
                                direction=direction,
                                mem=calculator.mem,
                                nprocshared=calculator.nprocshared,
                                scratch=calculator.scratch,
                                method=calculator.method,
                                basis=calculator.basis
                                )

        self.write_input(transitionstate, calc)
        label = self.submit_transitionstate(transitionstate, calc, "west")

        while not self.check_complete(label):
            time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info("{} failed the shell optimization".format(calc.label))
            return False

        transitionstate.ase_molecule = self.read_log(os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        #######################
        ##### For centers #####
        #######################
        calc = calculator.get_center_calc(transitionstate,
                                direction=direction,
                                mem=calculator.mem,
                                nprocshared=calculator.nprocshared,
                                scratch=calculator.scratch,
                                method=calculator.method,
                                basis=calculator.basis
                                )

        self.write_input(transitionstate, calc)
        label = self.submit_transitionstate(transitionstate, calc, "general")

        while not self.check_complete(label):
            time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info("{} failed the center optimization".format(calc.label))
            return False 

        transitionstate.ase_molecule = self.read_log(os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        #######################
        ##### For overall #####
        #######################
        calc = calculator.get_overall_calc(transitionstate,
                                direction=direction,
                                mem=calculator.mem,
                                nprocshared=calculator.nprocshared,
                                scratch=calculator.scratch,
                                method=calculator.method,
                                basis=calculator.basis
                                )

        self.write_input(transitionstate, calc)
        label = self.submit_transitionstate(transitionstate, calc, "general")

        while not self.check_complete(label):
            time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info("{} failed the overall optimization".format(calc.label))
            return False 

        transitionstate.ase_molecule = self.read_log(os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        logging.info("Calculations for {}_{}_{} were successful!".format(transitionstate.reaction_label, direction, transitionstate.index ))
        return True

    def calculate_reaction(self, reaction=None, conformer_calculator=None, calculator=None):
        """
        A method to run calculations for all tranitionstates for a reaction
        """

        logging.info("Calculating geometries for {}".format(reaction))

        if conformer_calculator:
            reaction.generate_conformers(calculator=conformer_calculator)


        currently_running = [] 
        processes = {}
        for direction, transitionstates in reaction.ts.items():

            for transitionstate in transitionstates:

                process = Process(target=self.calculate_transitionstate, args=(transitionstate, calculator))
                processes[process.name] = process
                

        for name, process in processes.items():
            print process
            while len(currently_running) >= 50:
                for running in currently_running:
                    if not running.is_alive():
                        currently_running.remove(name)
            process.start()
            #process.join()
            currently_running.append(name)

        complete = False
        print len(currently_running), currently_running
        while len(currently_running) > 0:
            for name, process in processes.items():
                if not (name in currently_running):
                    continue
                if not process.is_alive():
                    currently_running.remove(name)

        print "Here"
        lowest_energy_f = None
        lowest_energy = 1e5
        for direction, transitionstates in reaction.ts.items():
            for transitionstate in transitionstates:
                f = "{}_{}_{}.log".format(reaction.label, direction, transitionstate.index)
                path = os.path.join(calculator.scratch, "ts", reaction.label, "conformers", f)
                if not os.path.exists(path):
                    logging.info("It appears that {} failed...".format(f))
                    continue
                parser = ccread(f)
                energy = parser.scfenergies[-1]
                if energy < lowest_energy:
                    lowest_energy = energy
                    lowest_energy_f = f
        
        if lowest_energy_f is None:
            logging.info("No transition state for {} was successfully calculated :(".format(reaction))
            return False

        copyfile(
            os.path.join(calculator.scratch, "ts", reaction.label, "conformers", lowest_energy_f),
            os.path.join(calculator.scratch, "ts", reaction.label, lowest_energy_f[:-6] + ".log")
        )
        logging.info("The lowest energy file is {}".format(lowest_energy_f))


    def check_complete(self, label):
        """
        A method to determine if a job is still running
        """
        command = """squeue -n "{}" """.format(label)
        output = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE).communicate()[0]
        if len(output.split("\n")) <= 2:
            return True
        else:
            return False

    def calculate_species(self, species=None, conformer_calculator=None, calculator=None):

        logging.info("Calculating geometries for {}".format(species))

        if conformer_calculator:
            species.generate_conformers(calculator=conformer_calculator)

        to_calculate = []
        results = {}
        for smiles, conformers in species.conformers.items():
            results[smiles] = {}
            for conformer in conformers:
                calc = calculator.get_conformer_calc(
                    conformer=conformer,
                    mem=calculator.mem,
                    nprocshared=calculator.nprocshared,
                    scratch=calculator.scratch,
                    method=calculator.method,
                    basis=calculator.basis
                )
                to_calculate.append([conformer, calc])

        currently_running = []
        for conformer, calc in to_calculate:

            while len(currently_running) >= 50:
                for running_label in currently_running:
                    if self.check_complete(running_label):
                        currently_running.remove(running_label)

            label = self.submit_conformer(conformer, calc, "general")
            currently_running.append(label)
        scratch = calculator.scratch
        for conformer, calc in to_calculate:
            
            
            starting_molecule = RMGMolecule(SMILES=conformer.smiles)
            starting_molecule = starting_molecule.toSingleBonds()

            scratch_dir = os.path.join(
                calculator.scratch,
                "species",
                conformer.smiles,
                "conformers"
            )
            f = calc.label + ".log"
            if not os.path.exists(os.path.join(scratch_dir, f)):
                logging.info("It seems that {} was never run...".format(calc.label))
                result = False

            complete, converged = calculator.verify_output_file(
                os.path.join(scratch_dir, f)
            )

            if not complete:
                logging.info("It appears that {} was killed prematurely".format(calc.label))
                result = False
                

            elif not converged:
                logging.info("{} failed QM optimization".format(calc.label))
                result = False
                
            else:
                atoms = self.read_log(
                    os.path.join(scratch_dir,f)
                )

                starting_molecule = RMGMolecule(SMILES=conformer.smiles)
                starting_molecule = starting_molecule.toSingleBonds()

                test_molecule = RMGMolecule()
                test_molecule.fromXYZ(
                    atoms.arrays["numbers"],
                    atoms.arrays["positions"]
                )
                
                if not starting_molecule.isIsomorphic(test_molecule):
                    logging.info("Output geometry of {} is not isomorphic with input geometry".format(calc.label))
                    result = False
                else:
                    logging.info("{} was successful and was validated!".format(calc.label))
                    result = True

            if not result:
                fail_dir = os.path.join(scratch_dir, "failures")
                if not os.path.exists(fail_dir):
                    os.makedirs(os.path.join(scratch_dir, "failures"))
                move(
                    os.path.join(scratch_dir, f),
                    os.path.join(scratch_dir, "failures", f)
                )
            results[conformer.smiles][f] = result



        for smiles, smiles_results in results.items():

            scratch_dir = os.path.join(
                calculator.scratch,
                "species",
                conformer.smiles,
                "conformers"
            )

            log_file = os.path.join(scratch_dir, "validations.yaml")

            with open(log_file, "w") as to_write:
                yaml.dump(smiles_results, to_write, default_flow_style=False)

            lowest_energy_f = None
            lowest_energy = 1e5

            for f, result in smiles_results.items():
                if result:
                    parser = ccread(os.path.join(scratch_dir, f))
                    if parser.scfenergies[-1] < lowest_energy:
                        lowest_energy = parser.scfenergies[-1]
                        lowest_energy_f = f

            logging.info("The lowest energy conformer is {}".format(lowest_file_name))

            copyfile(
                os.path.join(scratch_dir, lowest_energy_f),
                os.path.join(
                    calculator.scratch,
                    "species",
                    label,
                    lowest_energy_f[:-6] + ".log")
            )

        

    """ Skipping over this for now
    def run_rotor(self, conformer, torsion, steps, step_size):
        
        #a method to run a hindered rotor calculation for a single rotor
        
        calc = calculator.get_rotor_calc(conformer=conformer, torsion=torsion, steps=steps, step_size=step_size)
        logging.info("Running hindered rotor calculation for {}".format(torsion))
        self.write_input(conformer, calc)
        self.submit(calc, "general")

    def run_rotors(self, conformer, steps, step_size):
        ""
        #A method to run hindered rotor scans for all torsions in a conformer
        

        for torsion in conformer.torsions:
            self.run_rotor(conformer=conformer, torsion=torsion, steps=steps, step_size=step_size)

    def calculate_rotors(self, conformer, calculator, steps, step_size):
        
        #A method to submit and verify all hindered rotor scans for a conformer
        
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
"""


