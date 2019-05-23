from autotst.calculator.gaussian import Gaussian
from autotst.calculator.vibrational_analysis import VibrationalAnalysis, percent_change
from autotst.calculator.statmech import StatMech
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
from cclib.io import ccread
import cclib
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.species import Species as RMGSpecies
from rmgpy.reaction import Reaction as RMGReaction, ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase
import rmgpy
from ase.calculators.gaussian import Gaussian as ASEGaussian
from ase.atoms import Atom, Atoms
import ase
import rdkit.Chem.rdDistGeom
import rdkit.DistanceGeometry
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Chem import AllChem
from rdkit import Chem
import rdkit
import os
import time
import yaml
from shutil import move, copyfile
import numpy as np
import pandas as pd
import subprocess
import multiprocessing
from multiprocessing import Process, Manager
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


class Job():
    """
    A class to deal with the input and output of calculations
    """

    def __init__(
            self,
            reaction=None,
            calculator=None, # An AutoTST Gaussian calculator with proper directory settings
            conformer_calculator=None,
            partition="general" # The partition to run calculations on
            ):

        self.reaction = reaction
        if isinstance(reaction, Reaction):
            self.label = self.reaction.label
        elif reaction is None:
            self.label = None
        else:
            assert False, "Reaction provided was not a reaction object"

        self.calculator = calculator
        if self.calculator:
            self.directory = self.calculator.directory
        self.conformer_calculator = conformer_calculator
        self.partition = partition

    def __repr__(self):
        return "< Job '{}'>".format(self.label)

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
        try:
            os.makedirs(ase_calculator.scratch)
        except OSError:
            pass

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

#################################################################################

    def submit_conformer(self, conformer):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        assert conformer, "Please provide a conformer to submit a job"

        self.calculator.conformer = conformer
        ase_calculator = self.calculator.get_conformer_calc()
        self.write_input(conformer, ase_calculator)

        label = conformer.smiles + "_{}".format(conformer.index)
        file_path = os.path.join(ase_calculator.scratch, label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path
        
        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info(
                "It appears that this job has already been run, not running it a second time.")

        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 -t 12:00:00 --mem=60GB $AUTOTST/autotst/job/submit.sh""".format(
                    label, self.partition), shell=True)

        return label

    def calculate_species(self, species):


        assert isinstance(species, Species), "Please provide a species object for this type of calculation"

        logging.info("Calculating geometries for {}".format(species))

        if self.conformer_calculator:
            species.generate_conformers(calculator=self.conformer_calculator)

        to_calculate = []
        results = {}
        for smiles, conformers in list(species.conformers.items()):
            results[smiles] = {}
            for conformer in conformers:
                to_calculate.append([conformer])

        currently_running = []
        for conformer in to_calculate:
            while len(currently_running) >= 50:
                for running_label in currently_running:
                    if self.check_complete(running_label):
                        currently_running.remove(running_label)

            label = self.submit_conformer(conformer)
            currently_running.append(label)

        while len(currently_running) != 0:
            for running_label in currently_running:
                if self.check_complete(running_label):
                    currently_running.remove(running_label)

        for conformer in to_calculate:

            starting_molecule = RMGMolecule(SMILES=conformer.smiles)
            starting_molecule = starting_molecule.toSingleBonds()

            scratch_dir = os.path.join(
                self.directory,
                "species",
                conformer.smiles,
                "conformers"
            )
            f = "{}_{}.log".format(conformer.smiles, conformer.index)
            if not os.path.exists(os.path.join(scratch_dir, f)):
                logging.info(
                    "It seems that {} was never run...".format(calc.label))
                result = False

            complete, converged = calculator.verify_output_file(
                os.path.join(scratch_dir, f)
            )

            if not complete:
                logging.info(
                    "It appears that {} was killed prematurely".format(f))
                result = False

            elif not converged:
                logging.info("{} failed QM optimization".format(f))
                result = False

            else:
                atoms = self.read_log(
                    os.path.join(scratch_dir, f)
                )

                starting_molecule = RMGMolecule(SMILES=conformer.smiles)
                starting_molecule = starting_molecule.toSingleBonds()

                test_molecule = RMGMolecule()
                test_molecule.fromXYZ(
                    atoms.arrays["numbers"],
                    atoms.arrays["positions"]
                )

                if not starting_molecule.isIsomorphic(test_molecule):
                    logging.info(
                        "Output geometry of {} is not isomorphic with input geometry".format(calc.label))
                    result = False
                else:
                    logging.info(
                        "{} was successful and was validated!".format(calc.label))
                    result = True

            if not result:
                fail_dir = os.path.join(scratch_dir, "failures")
                try:
                    os.makedirs(os.path.join(scratch_dir, "failures"))
                except OSError:
                    logging.info("{} already exists...".format(fail_dir))
                move(
                    os.path.join(scratch_dir, f),
                    os.path.join(scratch_dir, "failures", f)
                )
            results[conformer.smiles][f] = result

        for smiles, smiles_results in list(results.items()):
            scratch_dir = os.path.join(
                self.directory,
                "species",
                conformer.smiles,
                "conformers"
            )

            log_file = os.path.join(scratch_dir, "validations.yaml")

            with open(log_file, "w") as to_write:
                yaml.dump(smiles_results, to_write, default_flow_style=False)

            lowest_energy_f = None
            lowest_energy = 1e5

            for f, result in list(smiles_results.items()):
                if result:
                    parser = ccread(os.path.join(scratch_dir, f))
                    if parser.scfenergies[-1] < lowest_energy:
                        lowest_energy = parser.scfenergies[-1]
                        lowest_energy_f = f

            logging.info(
                "The lowest energy conformer is {}".format(lowest_energy_f))

            copyfile(
                os.path.join(scratch_dir, lowest_energy_f),
                os.path.join(
                    self.directory,
                    "species",
                    conformer.smiles,
                    conformer.smiles + ".log")
            )

#################################################################################

    def submit_transitionstate(self, transitionstate, opt_type, restart=False):
        """
        A methods to submit a job for a TS object based on a single calculator
        """
        assert transitionstate, "Please provide a transitionstate to submit a job"
        
        if opt_type.lower() == "shell":
            ase_calculator = self.calculator.get_shell_calc()
            time = "12:00:00"
        elif opt_type.lower() == "center":
            ase_calculator = self.calculator.get_center_calc()
            time = "12:00:00"
        elif opt_type.lower() == "overall":
            ase_calculator = self.calculator.get_overall_calc()
            time = "12:00:00"
        elif opt_type.lower() == "irc":
            ase_calculator = self.calculator.get_irc_calc()
            time = "24:00:00"

        self.write_input(transitionstate, ase_calculator)

        label = ase_calculator.label
        scratch = ase_calculator.scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info("It appears that {} has already been attempted...".format(label))

        if (not attempted) or restart:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=60GB -t {2} $AUTOTST/autotst/job/submit.sh""".format(
                    label, self.partition, time), shell=True)

        return label

    def calculate_transitionstate(self, transitionstate, vibrational_analysis=True):
        """
        A method to perform the partial optimizations for a transitionstate and arrive
        at a final geometry. Returns True if we arrived at a final geometry, returns false
        if there is an error along the way.
        """

        ts_identifier = "{}_{}_{}".format(
            transitionstate.reaction_label, transitionstate.direction, transitionstate.index)

        for opt_type in ["shell", "center", "overall"]:
            self.calculator.conformer = transitionstate

            if opt_type == "overall":
                 file_path = "{}_{}_{}.log".format(transitionstate.reaction_label, transitionstate.direction, transitionstate.index)
            else:
                 file_path = "{}_{}_{}_{}.log".format(transitionstate.reaction_label, transitionstate.direction, opt_type, transitionstate.index)

            file_path = os.path.join(
                self.directory, 
                "ts", 
                transitionstate.reaction_label, 
                "conformers", 
                file_path
            )


            if not os.path.exists(file_path):
                logging.info(
                    "Submitting {} calculations for {}".format(opt_type.upper(),ts_identifier))
                label = self.submit_transitionstate(
                    transitionstate, opt_type=opt_type.lower())
                print label
                while not self.check_complete(label):
                    time.sleep(15)

            else:
                logging.info(
                    "It appears that we already have a complete {} log file for {}".format(opt_type.upper(), ts_identifier))

                complete, converged = self.calculator.verify_output_file(file_path)
                
                if not complete:
                    logging.info(
                        "It seems that the {} file never completed for {} never completed, running it again".format(opt_type.upper(), ts_identifier))
                    label = self.submit_transitionstate(
                        transitionstate, opt_type=opt_type.lower(), restart=True)
                    while not self.check_complete(label):
                        time.sleep(15)

            complete, converged = self.calculator.verify_output_file(file_path)

            if not (complete and converged):
                logging.info(
                    "{} failed the {} optimization".format(ts_identifier, opt_type.upper()))
                results[ts_identifier] = False
                return False
            logging.info(
                "{} successfully completed the {} optimization!".format(ts_identifier, opt_type.upper()))
            transitionstate.ase_molecule = self.read_log(file_path)
            transitionstate.update_coords_from("ase")

        logging.info(
            "Calculations for {} are complete and resulted in a normal termination!".format(ts_identifier))

        got_one = self.validate_transitionstate(
                transitionstate=transitionstate, vibrational_analysis=vibrational_analysis)
        if got_one:
            results[ts_identifier] = True
            return True
        else:
            results[ts_identifier] = False
            return False

    def calculate_reaction(self, vibrational_analysis=True):
        """
        A method to run calculations for all tranitionstates for a reaction
        """

        logging.info("Calculating geometries for {}".format(self.reaction))

        if self.conformer_calculator:
            self.reaction.generate_conformers(ase_calculator=self.conformer_calculator)

        currently_running = []
        processes = {}
        manager = multiprocessing.Manager()
        results = manager.dict()
        global results

        for direction, transitionstates in list(self.reaction.ts.items()):

            for transitionstate in transitionstates:

                process = Process(target=self.calculate_transitionstate, args=(
                    transitionstate,))
                processes[process.name] = process

        for name, process in list(processes.items()):
            while len(currently_running) >= 50:
                for running in currently_running:
                    if not running.is_alive():
                        currently_running.remove(name)
            process.start()
            currently_running.append(name)

        while len(currently_running) > 0:
            for name, process in list(processes.items()):
                if not (name in currently_running):
                    continue
                if not process.is_alive():
                    currently_running.remove(name)

        energies = []
        for label, result in results.items():
            if not result:
                logging.info("Calculations for {} FAILED".format(label))
            f = "{}.log".format(label)
            path = os.path.join(self.calculator.directory, "ts",
                    self.reaction.label, "conformers", f)
            if not os.path.exists(path):
                logging.info("It appears that {} failed...".format(f))
                continue
            try:
                parser = ccread(path)
                if parser is None:
                    logging.info(
                        "Something went wrong when reading in results for {}...".format(f))
                    continue
                energy = parser.scfenergies[-1]
            except:
                logging.info(
                    "The parser does not have an scf energies attribute, we are not considering {}".format(f))
                energy = 1e5

            energies.append([energy, transitionstate, f])

        energies = pd.DataFrame(
            energies, columns=["energy", "transitionstate", "file"]).sort_values("energy")

        if energies.shape[0] == 0:
            logging.info(
                "No transition state for {} was successfully calculated... :(".format(self.reaction))
            return False

        energies.reset_index(inplace=True)
        lowest_energy_label = energies.iloc[0].file
        logging.info("The lowest energy transition state is {}".format(lowest_energy_label))

        copyfile(
            os.path.join(self.calculator.directory, "ts", self.reaction.label,
                         "conformers", lowest_energy_label + ".log"),
            os.path.join(self.calculator.directory, "ts",
                         self.reaction.label, self.reaction.label + ".log")
        )
        logging.info("The lowest energy file is {}! :)".format(
            lowest_energy_label + ".log"))
        return True

    def validate_transitionstate(self, transitionstate, vibrational_analysis=True):

        validated = False
        if vibrational_analysis:
            vib = VibrationalAnalysis(
                transitionstate=transitionstate, directory=self.directory)
            validated = vib.validate_ts()
        if not validated:
            logging.info("Could not validate with Vibrational Analysis... Running an IRC to validate instead...")
            label = self.submit_transitionstate(
                transitionstate, opt_type="irc")
            while not self.check_complete(label):
                time.sleep(15)
            result = calculator.validate_irc(calc=irc_calc)
            if result:
                logging.info("Validated via IRC")
                return True
            else:
                logging.info(
                    "Could not validate this conformer... trying the next lowest energy conformer")
                return False
        else:
            logging.info("Validated via Vibrational Analysis")
            return True

#################################################################################

    def submit_rotor(self, conformer, torsion_index):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        assert conformer, "Please provide a conformer to submit a job"

        ase_calculator = self.calculator.get_rotor_calc(conformer, torsion_index)

        self.write_input(conformer, ase_calculator)

        file_path = os.path.join(ase_calculator.scratch, ase_calculator.label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info(
                "It appears that this job has already been run, not running it a second time.")

        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 -t 8:00:00 --mem=100000 $AUTOTST/autotst/job/submit.sh""".format(
                    label, self.partition), shell=True)

        return label

    def calculate_rotors(self, conformer, steps=36, step_size=10.0):

        complete = {}
        calculators = {}
        verified = {}
        if len(conformer.torsions) == 0:
            logging.info("No torsions to run scans on.")
            return {}

        for torsion in conformer.torsions:
            label = self.submit_rotor(
                conformer, torsion.index)
            logging.info(label)
            complete[label] = False
            verified[label] = False

        done = False
        lowest_energy_label = None
        conformer_error = False

        while not done:
            for label in list(complete.keys()):
                if not self.check_complete(label):
                    continue
                if done:
                    continue
                complete[label] = True
                lowest_conf, continuous, good_slope, opt_count_check = self.verify_rotor( ##################################
                    conformer, label)
                if all([lowest_conf, continuous]):
                    verified[label] = True
                else:
                    verified[label] = False

                if not lowest_conf:
                    done = True
                    lowest_energy_label = label
                    conformer_error = True
                    continue
                elif all(complete.values()):
                    done = True

        if conformer_error:
            logging.info(
                "A lower energy conformer was found... Going to optimize this insted")
            for label in list(complete.keys()):
                subprocess.call("""scancel -n '{}'""".format(label), shell=True)
            if isinstance(conformer, TS):
                file_name = os.path.join(
                    self.directory, "ts", conformer.reaction_label, "rotors", lowest_energy_label + ".log")
            else:
                file_name = os.path.join(
                    self.directory, "species",conformer.smiles , "rotors", lowest_energy_label + ".log")
            parser = ccread(file_name)
            first_is_lowest, min_energy, atomnos, atomcoords = self.check_rotor_lowest_conf(
                parser=parser)
            symbol_dict = {
                17: "Cl",
                9:  "F",
                8:  "O",
                7:  "N",
                6:  "C",
                1:  "H",
            }
            atoms = []
            for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
                atoms.append(
                    Atom(symbol=symbol_dict[atom_num], position=coords))
            conformer.ase_molecule = Atoms(atoms)
            conformer.update_coords_from("ase")
            for index in ["X", "Y", "Z"]:
                if index != conformer.index:
                    logging.info("Setting index of {} to {}...".format(conformer, index))
                    conformer.index = index
                    break

            label = self.submit_conformer(conformer)

            while not self.check_complete(label):
                time.sleep(15)

            logging.info(
                "Reoptimization complete... performing hindered rotors scans again")
            return self.calculate_rotors(conformer, steps, step_size)

        else:
            for label, boolean in list(verified.items()):
                if not boolean:
                    try:
                        if isinstance(conformer, TS):
                            file_path = os.path.join(
                                self.directory, "ts", conformer.reaction_label, "rotors")
                        else:
                            file_path = os.path.join(
                                self.directory, "species",conformer.smiles , "rotors")

                        os.mkdirs(os.path.join(file_path, failures))
                    except:
                        pass
                    move(
                        os.path.join(file_path, label + ".log"),
                        os.path.join(file_path, "failures",
                                     label + ".log")
                    )
            return verified

    def verify_rotor(self, conformer, label, steps=36, step_size=10.0):

        if isinstance(conformer, TS):
            file_name = os.path.join(
                self.directory, "ts", conformer.reaction_label, "rotors", label  + ".log")
        elif isinstance(conformer, Conformer):
             file_name = os.path.join(
                self.directory, "species", conformer.smiles, "rotors", label  + ".log")           
        parser = cclib.io.ccread(file_name)

        continuous = self.check_rotor_continuous(
            steps, step_size, parser=parser)
        good_slope = self.check_rotor_slope(steps, step_size, parser=parser)
        [lowest_conf, energy, atomnos,
            atomcoords] = self.check_rotor_lowest_conf(parser=parser)
        opt_count_check = self.check_rotor_opts(steps, parser=parser)

        return [lowest_conf, continuous, good_slope, opt_count_check]

    def check_rotor_opts(self, steps, parser):


        #opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_indices = [i for i, status in enumerate(
            parser.optstatus) if status > 1]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]

        n_opts_check = (steps + 1) == len(opt_SCFEnergies)

        return n_opts_check

    def check_rotor_slope(self, steps, step_size, parser, tol=0.1):


        opt_indices = [i for i, status in enumerate(
            parser.optstatus) if status in [2, 4]]
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

    def check_rotor_continuous(self, steps, step_size, parser, tol=0.1):

        assert isinstance(step_size, float)

        opt_indices = [i for i, status in enumerate(
            parser.optstatus) if status in [2, 4]]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]

        max_energy = max(opt_SCFEnergies)
        min_energy = min(opt_SCFEnergies)
        energy_tol = np.absolute(tol*(max_energy - min_energy))

        checked = [None for angle in range(0, 360)]

        continuous = True

        for step, energy in enumerate(opt_SCFEnergies):
            abs_theta = int(step*step_size)
            theta = abs_theta % 360

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

        return continuous

    def check_rotor_lowest_conf(self, parser, tol=0.1):

        opt_indices = [i for i, status in enumerate(
            parser.optstatus) if status in [2, 4]]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]

        max_energy = max(opt_SCFEnergies)
        min_energy = min(opt_SCFEnergies)
        energy_tol = tol*(max_energy - min_energy)

        first_is_lowest = True  # Therefore...
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
