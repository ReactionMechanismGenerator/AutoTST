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
            calculator=None,
            conformer_calculator=None):

        self.reaction = reaction
        if isinstance(reaction, Reaction):
            self.label = reaction.label
        elif reaction is None:
            self.label = None
        else:
            assert False, "Reaction provided was not a reaction object"

        self.calculator = calculator
        self.conformer_calculator = conformer_calculator

    def __repr__(self):
        return "< Job '{}'>".format(self.label)

    def read_log(self, file_path=None):
        """
        A helper method that allows one to easily parse log files
        """
        symbol_dict = {
            35: "Br",
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

    def submit_conformer(self, conformer, ase_calculator, partition="general"):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        assert conformer, "Please provide a conformer to submit a job"

        self.write_input(conformer, ase_calculator)

        label = conformer.smiles + "_{}".format(conformer.index)
        scratch = ase_calculator.scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info(
                "It appears that this job has already been run, not running it a second time.")

        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=60GB $AUTOTST/autotst/job/submit.sh""".format(
                    label, partition), shell=True)

        return label

    def calculate_conformer(self, conformer, calculator):

        if isinstance(calculator,Gaussian):
            calc = calculator.get_conformer_calc(conformer=conformer,convergence='Tight')

        else:
            calc = calculator.get_conformer_calc(conformer=conformer)

        label = self.submit_conformer(conformer, calc, "general")

        scratch = calculator.scratch

        scratch_dir = os.path.join(
            calculator.scratch,
            "species",
            conformer.smiles,
            "conformers")

        f = calc.label + ".log"
        if not os.path.exists(os.path.join(scratch_dir, f)):
            logging.info(
                "It seems that {} was never run...".format(calc.label))
            result = False

        complete, converged = calculator.verify_output_file(os.path.join(scratch_dir, f))

        if (complete and converged):
            logging.info(
            "{} was successful and was validated!".format(calc.label))
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

        if not complete:
            logging.info(
                "It appears that {} was killed prematurely".format(calc.label))
            result = False

        elif not converged:
            logging.info("{} did not converge".format(calc.label))
            result = False
            if isinstance(calculator,Gaussian):
                if not calculator.convergence.lower() in ["tight", "verytight", "loose"]:
                    logging.info("{} failed QM optimization".format(calc.label))
                else:
                    logging.info("Resubmitting {} with default convergence criteria".format(calc.label))
                    atoms = self.read_log(os.path.join(scratch_dir, f))
                    conformer.ase_molecule = atoms
                    conformer.update_coords_from("ase")
                    calc = calculator.get_conformer_calc(
                        conformer,
                        convergence=""
                    )
                    logging.info("Removing the old log file that didn't converge, restarting from last geometry")
                    os.remove(os.path.join(scratch_dir, f))

                    label = self.submit_conformer(conformer, calc, "general")

                    scratch_dir = os.path.join(
                        calculator.scratch,
                        "species",
                        conformer.smiles,
                        "conformers"
                    )
                    f = calc.label + ".log"
                    if not os.path.exists(os.path.join(scratch_dir, f)):
                        logging.info(
                        "It seems that {} was never run...".format(calc.label))
                        result = False

                    complete, converged = calculator.verify_output_file(
                        os.path.join(scratch_dir, f)
                    )

                    if not complete:
                        logging.info(
                        "It appears that {} was killed prematurely".format(calc.label))
                        result = False

                    elif not converged:
                        logging.info("{} failed second QM optimization".format(calc.label))
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
            return False

        return True


    def calculate_species(self, species=None, conformer_calculator=None, calculator=None):

        logging.info("Calculating geometries for {}".format(species))

        if conformer_calculator:
            species.generate_conformers(calculator=conformer_calculator)

        currently_running = []
        processes = {}
        for smiles, conformers in list(species.conformers.items()):

            for conformer in conformers:

                process = Process(target=self.calculate_conformer, args=(
                conformer, calculator))
                processes[process.name] = process

        for name, process in list(processes.items()):
            while len(currently_running) >= 50:
                for running in currently_running:
                    if not running.is_alive():
                        currently_running.remove(name)
            process.start()
            # process.join()
            currently_running.append(name)

        complete = False
        while len(currently_running) > 0:
            for name, process in list(processes.items()):
                if not (name in currently_running):
                    continue
                if not process.is_alive():
                    currently_running.remove(name)

        results = []
        for smiles, conformers in list(species.conformers.items()):
            for conformer in conformers:
                scratch_dir = os.path.join(
                    calculator.scratch,
                    "species",
                    conformer.smiles,
                    "conformers"
                )
                f = "{}_{}.log".format(conformer.smiles, conformer.index)
                path = os.path.join(scratch_dir, f)
                if not os.path.exists(os.path.join(scratch_dir, f)):
                    logging.info(
                        "It seems that {} was never run...".format(calc.label))
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

                results.append([energy, conformer, f])

        results = pd.DataFrame(
            results, columns=["energy", "conformer", "file"]).sort_values("energy").reset_index()

        if results.shape[0] == 0:
            logging.info(
                "No conformer for {} was successfully calculated... :(".format(reaction))
            return False

        for index in range(results.shape[0]):
            conformer = results.conformer[index]
            lowest_energy_file = results.file[index]
            break

        logging.info(
            "The lowest energy conformer is {}".format(lowest_energy_file))

        dest = os.path.join(calculator.scratch,"species",conformer,lowest_energy_file + ".log")

        try:
            copyfile(os.path.join(scratch_dir, lowest_energy_file),dest)
        except IOError:
            os.makedirs(os.path.dirname(dest))
            copyfile(os.path.join(scratch_dir, lowest_energy_file),dest)

        logging.info("The lowest energy file is {}! :)".format(
            lowest_energy_file))

        return True
#################################################################################

    def submit_transitionstate(self, transitionstate, ase_calculator, partition="general"):
        """
        A methods to submit a job for a TS object based on a single calculator
        """
        assert transitionstate, "Please provide a transitionstate to submit a job"

        self.write_input(transitionstate, ase_calculator)
        label = ase_calculator.label
        scratch = ase_calculator.scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info("It appears that this job has already been run")

        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=60GB $AUTOTST/autotst/job/submit.sh""".format(
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

        ts_identifier = "{}_{}_{}".format(
            transitionstate.reaction_label, direction, transitionstate.index)

        calc = calculator.get_shell_calc(transitionstate,
                                         direction=direction
                                         )

        if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
            logging.info(
                "Submitting SHELL calculations for {}".format(ts_identifier))
            self.write_input(transitionstate, calc)
            label = self.submit_transitionstate(
                transitionstate, calc, "general")
            while not self.check_complete(label):
                time.sleep(15)

        else:
            logging.info(
                "It appears that we already have a complete SHELL log file for {}".format(ts_identifier))

            complete, converged = calculator.verify_output_file(
                os.path.join(calc.scratch, calc.label + ".log")
            )
            if not complete:
                logging.info(
                    "It seems that the SHELL file never completed for {} never completed, running it again".format(ts_identifier))
                self.write_input(transitionstate, calc)
                label = self.submit_transitionstate(
                    transitionstate, calc, "general")
                while not self.check_complete(label):
                    time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info(
                "{} failed the SHELL optimization".format(ts_identifier))
            return False
        logging.info(
            "{} successfully completed the SHELL optimization!".format(ts_identifier))
        transitionstate.ase_molecule = self.read_log(
            os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        #######################
        ##### For centers #####
        #######################
        calc = calculator.get_center_calc(transitionstate,
                                          direction=direction
                                          )

        if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
            logging.info(
                "Submitting CENTER calculations for {}".format(ts_identifier))
            self.write_input(transitionstate, calc)
            label = self.submit_transitionstate(
                transitionstate, calc, "general")
            while not self.check_complete(label):
                time.sleep(15)

        else:
            logging.info(
                "It appears that we already have a complete CENTER log file for {}".format(ts_identifier))

            complete, converged = calculator.verify_output_file(
                os.path.join(calc.scratch, calc.label + ".log")
            )
            if not complete:
                logging.info(
                    "It seems that the CENTER file never completed for {} never completed, running it again".format(ts_identifier))
                self.write_input(transitionstate, calc)
                label = self.submit_transitionstate(
                    transitionstate, calc, "general")
                while not self.check_complete(label):
                    time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info(
                "{} failed the CENTER optimization".format(ts_identifier))
            return False
        logging.info(
            "{} successfully completed the CENTER optimization!".format(ts_identifier))
        transitionstate.ase_molecule = self.read_log(
            os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        #######################
        ##### For overall #####
        #######################
        calc = calculator.get_overall_calc(transitionstate,
                                           direction=direction
                                           )

        if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
            logging.info(
                "Submitting OVERALL calculations for {}".format(ts_identifier))
            self.write_input(transitionstate, calc)
            label = self.submit_transitionstate(
                transitionstate, calc, "general")
            while not self.check_complete(label):
                time.sleep(15)

        else:
            logging.info(
                "It appears that we already have a complete OVERALL log file for {}".format(ts_identifier))

            complete, converged = calculator.verify_output_file(
                os.path.join(calc.scratch, calc.label + ".log")
            )
            if not complete:
                logging.info(
                    "It seems that the OVERALL file never completed for {} never completed, running it again".format(ts_identifier))
                self.write_input(transitionstate, calc)
                label = self.submit_transitionstate(
                    transitionstate, calc, "general")
                while not self.check_complete(label):
                    time.sleep(15)

        complete, converged = calculator.verify_output_file(
            os.path.join(calc.scratch, calc.label + ".log")
        )

        if not (complete and converged):
            logging.info(
                "{} failed the OVERALL optimization".format(ts_identifier))
            return False
        logging.info(
            "{} successfully completed the OVERALL optimization!".format(ts_identifier))
        transitionstate.ase_molecule = self.read_log(
            os.path.join(calc.scratch, calc.label + ".log"))
        transitionstate.update_coords_from("ase")

        logging.info(
            "Calculations for {} are complete and resulted in a normal termination!".format(ts_identifier))
        return True

    def calculate_reaction(self, reaction=None, conformer_calculator=None, calculator=None, vibrational_analysis=True):
        """
        A method to run calculations for all tranitionstates for a reaction
        """

        logging.info("Calculating geometries for {}".format(reaction))

        if conformer_calculator:
            reaction.generate_conformers(calculator=conformer_calculator)

        currently_running = []
        processes = {}
        for direction, transitionstates in list(reaction.ts.items()):

            for transitionstate in transitionstates:

                process = Process(target=self.calculate_transitionstate, args=(
                    transitionstate, calculator))
                processes[process.name] = process

        for name, process in list(processes.items()):
            while len(currently_running) >= 50:
                for running in currently_running:
                    if not running.is_alive():
                        currently_running.remove(name)
            process.start()
            # process.join()
            currently_running.append(name)

        complete = False
        while len(currently_running) > 0:
            for name, process in list(processes.items()):
                if not (name in currently_running):
                    continue
                if not process.is_alive():
                    currently_running.remove(name)

        results = []
        for direction, transitionstates in list(reaction.ts.items()):
            for transitionstate in transitionstates:
                f = "{}_{}_{}.log".format(
                    reaction.label, direction, transitionstate.index)
                path = os.path.join(calculator.scratch, "ts",
                                    reaction.label, "conformers", f)
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

                results.append([energy, transitionstate, f])

        results = pd.DataFrame(
            results, columns=["energy", "transitionstate", "file"]).sort_values("energy").reset_index()

        if results.shape[0] == 0:
            logging.info(
                "No transition state for {} was successfully calculated... :(".format(reaction))
            return False
        got_one = False
        for index in range(results.shape[0]):
            ts = results.transitionstate[index]
            got_one = self.validate_transitionstate(
                transitionstate=ts, calculator=calculator, vibrational_analysis=vibrational_analysis)
            if got_one:
                lowest_energy_file = results.file[index]
                break

        if not got_one:
            logging.info(
                "None of the transition states could be validated... :(")
            return False

        copyfile(
            os.path.join(calculator.scratch, "ts", reaction.label,
                         "conformers", lowest_energy_file),
            os.path.join(calculator.scratch, "ts",
                         reaction.label, reaction.label + ".log")
        )
        logging.info("The lowest energy file is {}! :)".format(
            lowest_energy_file))
        return True

    def validate_transitionstate(self, transitionstate, calculator, vibrational_analysis=True):

        validated = False
        if vibrational_analysis:
            vib = VibrationalAnalysis(
                ts=transitionstate, scratch=calculator.scratch)
            validated = vib.validate_ts()
        if not validated:
            logging.info("Running an IRC to validate")
            irc_calc = calculator.get_irc_calc(
                ts=transitionstate
            )
            label = self.submit_transitionstate(
                transitionstate, irc_calc, "general")
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

    def submit_rotor(self, conformer, ase_calculator, partition="general"):
        """
        A methods to submit a job based on the calculator and partition provided
        """
        assert conformer, "Please provide a conformer to submit a job"

        self.write_input(conformer, ase_calculator)

        label = ase_calculator.label
        scratch = ase_calculator.scratch
        file_path = os.path.join(scratch, label)

        os.environ["COMMAND"] = "g16"  # only using gaussian for now
        os.environ["FILE_PATH"] = file_path

        attempted = False
        if os.path.exists(file_path + ".log"):
            attempted = True
            logging.info(
                "It appears that this job has already been run, not running it a second time.")

        if not attempted:
            subprocess.call(
                """sbatch --exclude=c5003,c3040 --job-name="{0}" --output="{0}.slurm.log" --error="{0}.slurm.log" -p {1} -N 1 -n 20 --mem=100000 $AUTOTST/autotst/job/submit.sh""".format(
                    label, partition), shell=True)

        return label

    def calculate_rotors(self, conformer, calculator, steps=36, step_size=10.0):

        complete = {}
        calculators = {}
        verified = {}
        for torsion in conformer.torsions:
            calc = calculator.get_rotor_calc(
                conformer=conformer,
                torsion=torsion,
                steps=steps,
                step_size=step_size,
            )
            label = self.submit_rotor(
                conformer=conformer, ase_calculator=calc, partition="general")
            logging.info(label)
            complete[label] = False
            calculators[label] = calc
            verified[label] = False

        done = False
        lowest_energy_label = None
        conformer_error = False
        if len(conformer.torsions) == 0:
            logging.info("No torsions to run scans on.")
            return {}
        while not done:
            for label in list(complete.keys()):
                if not self.check_complete(label):
                    continue
                if done:
                    continue
                complete[label] = True
                ase_calc = calculators[label]
                lowest_conf, continuous, good_slope, opt_count_check = self.verify_rotor(
                    steps=steps, step_size=step_size, ase_calculator=ase_calc)
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
                command = """scancel -n '{}'""".format(label)
            ase_calculator = calculators[label]
            file_name = os.path.join(
                ase_calculator.scratch, lowest_energy_label + ".log")
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

            if isinstance(conformer, TS):
                calc = calculator.get_overall_calc(conformer,
                                                   direction=conformer.direction
                                                   )

                calc.scratch = calc.scratch.strip("/conformers")
                conformer.direction = "forward"
                conformer.index = "X"
                label = self.submit_transitionstate(
                    transitionstate=conformer, ase_calculator=calc)
            else:
                calc = calculator.get_conformer_calc(conformer
                                                     )
                calc.scratch = calc.scratch.strip("/conformers")
                conformer.index = "X"
                label = self.submit_conformer(conformer, calc, "general")

            while not self.check_complete(label):
                os.sleep(15)

            logging.info(
                "Reoptimization complete... performing hindered rotors scans again")
            return self.calculate_rotors(conformer, calculator, steps, step_size)

        else:
            for label, boolean in list(verified.items()):
                if not boolean:
                    calc = calculators[label]
                    try:
                        os.mkdir(os.path.join(calc.scratch, "failures"))
                    except:
                        pass
                    move(
                        os.path.join(calc.scratch, calc.label + ".log"),
                        os.path.join(calc.scratch, "failures",
                                     calc.label + ".log")
                    )
            return verified

    def verify_rotor(self, steps=36, step_size=10.0, ase_calculator=None):

        assert (ase_calculator is not None)

        file_name = os.path.join(
            ase_calculator.scratch, ase_calculator.label + ".log")
        parser = cclib.io.ccread(file_name)

        continuous = self.check_rotor_continuous(
            steps, step_size, parser=parser)
        good_slope = self.check_rotor_slope(steps, step_size, parser=parser)
        [lowest_conf, energy, atomnos,
            atomcoords] = self.check_rotor_lowest_conf(parser=parser)
        opt_count_check = self.check_rotor_opts(steps, parser=parser)

        return [lowest_conf, continuous, good_slope, opt_count_check]

    def check_rotor_opts(self, steps, parser=None):

        assert (parser is not None)

        if parser is None:
            if file_name is None:
                assert ase_calculator is not None
                file_name = os.path.join(
                    ase_calculator.scratch, ase_calculator.label + ".log")

            parser = cclib.io.ccread(file_name)

        #opt_indices = [i for i, status in enumerate(parser.optstatus) if status==2]
        opt_indices = [i for i, status in enumerate(
            parser.optstatus) if status > 1]
        opt_SCFEnergies = [parser.scfenergies[index] for index in opt_indices]

        n_opts_check = (steps + 1) == len(opt_SCFEnergies)

        return n_opts_check

    def check_rotor_slope(self, steps, step_size, parser=None, tol=0.1):

        assert (parser is not None)

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

    def check_rotor_continuous(self, steps, step_size, parser=None, tol=0.1):

        assert isinstance(step_size, float)
        assert (parser is not None)

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

    def check_rotor_lowest_conf(self, parser=None, tol=0.1):

        assert (parser is not None)
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
