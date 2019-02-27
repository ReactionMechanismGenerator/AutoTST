#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2018 Prof. Richard H. West (r.west@northeastern.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
##########################################################################

import os
import itertools
import logging
import numpy as np
from cclib.io import ccread
from rdkit import Chem

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.reaction import Reaction as RMGReaction

import autotst
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.calculators.calculator import Calculator


from cclib.io import ccread

from ase.io.gaussian import read_gaussian, read_gaussian_out
from ase.calculators.gaussian import Gaussian as ASEGaussian


class Gaussian(Calculator):

    def __init__(self,
                 conformer=None,
                 mem="5GB",
                 nprocshared=20,
                 scratch=".",
                 method="m062x",
                 basis="6-311+g(2df,2p)",
                 save_directory="."):

        self.command = "g16"

        assert isinstance(conformer, (type(None), Conformer)
                          ), "Please provide a Conformer object"
        self.conformer = conformer
        self.mem = mem
        self.nprocshared = nprocshared
        self.scratch = scratch
        self.method = method
        self.basis = basis
        self.save_directory = save_directory


    def __repr__(self):
        if not self.conformer:
            return '<Gaussian Calculator "">'.format(None)
        return '<Gaussian Calculator "{0}">'.format(self.conformer.smiles)

    @property
    def label(self):
        if isinstance(self.conformer, TS):
            return self.conformer.reaction_label + \
                "_{}".format(self.conformer.index)
        elif isinstance(self.conformer, Conformer):
            label = Chem.rdinchi.InchiToInchiKey(Chem.MolToInchi(
                Chem.MolFromSmiles(self.conformer.smiles))).strip("-N")
            label += "_{}".format(self.conformer.index)
            return label
        else:
            return None

    def get_rotor_calc(self,
                       conformer=None,
                       torsion=None,
                       mem="5GB",
                       nprocshared=20,
                       scratch=".",
                       method="m062x",
                       basis="6-311+g(2df,2p)",
                       steps=36,
                       step_size=10.0):
        """
        A method to create all of the calculators needed to perform hindered rotor calculations given a conformer and a torsion
        """

        assert (torsion and (isinstance(torsion, Torsion))
                ), "To create a rotor calculator, you must provide a Torsion object."

        if not conformer:
            if not self.conformer:
                return None
            conformer = self.conformer

        assert isinstance(
            conformer, Conformer), "A Conformer object was not provided..."

        string = ""
        for bond in conformer.bonds:
            i, j = bond.atom_indices
            string += "B {} {}\n".format(i + 1, j + 1)

        i, j, k, l = torsion.atom_indices
        string += "D {} {} {} {} S {} {}".format(i + 1, j + 1, k + 1, l + 1, steps, float(step_size))

        if isinstance(conformer, TS):
            label = self.label + "_tor_{}_{}".format(j, k)
            t = "ts"
            d = self.label
        elif isinstance(conformer, Conformer):
            label = d = Chem.rdinchi.InchiToInchiKey(Chem.MolToInchi(
                Chem.MolFromSmiles(conformer.smiles))).strip("-N")
            label += "_tor{}{}".format(j, k)
            t = "species"
            
        conformer.rmg_molecule.updateMultiplicity()
        mult = conformer.rmg_molecule.multiplicity

        new_scratch = os.path.join(
                scratch,
                t,
                d,
                "rotors"
            )

        calc = ASEGaussian(mem=mem,
                           nprocshared=nprocshared,
                           label=label,
                           scratch=new_scratch,
                           method=method,
                           basis=basis,
                           extra="Opt=(CalcFC,ModRedun)",
                           multiplicity=mult,
                           addsec=[string])

        calc.atoms = conformer.ase_molecule
        del calc.parameters['force']

        return calc

    def get_species_calc(self,
                         conformer=None,
                         mem="5GB",
                         nprocshared=20,
                         scratch=".",
                         method="m062x",
                         basis="6-311+g(2df,2p)"):
        "A method that creates a calculator for a reactant or product that will perform a geometry optimization"

        if not conformer:
            if not self.conformer:
                return None
            conformer = self.conformer

        if isinstance(conformer, TS):
            logging.info(
                "TS object provided, cannot obtain a species calculator for a TS")
            return None

        assert isinstance(
            conformer, Conformer), "A Conformer object was not provided..."

        conformer.rmg_molecule.updateMultiplicity()

        # using this round about way of doing stuff because rmg's
        # `toAugumentedInChIKey` method doesn't work on our cluster

        smiles = conformer.rmg_molecule.toSMILES()
        short_label = Chem.rdinchi.InchiToInchiKey(Chem.MolToInchi(
            Chem.MolFromSmiles(smiles))).strip("-N") 
        label = short_label + "_{}".format(conformer.index)

        new_scratch = os.path.join(
                scratch,
                "species",
                short_label,
                "conformers"
            )

        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)

        calc = ASEGaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=method,
            basis=basis,
            extra="opt=(calcfc,verytight,gdiis,maxcycles=900) freq IOP(2/16=3)",
            multiplicity=conformer.rmg_molecule.multiplicity)
        calc.atoms = conformer.ase_molecule
        del calc.parameters['force']

        return calc

    def get_shell_calc(self,
                       ts=None,
                       mem="5GB",
                       nprocshared=20,
                       scratch=".",
                       method="m062x",
                       basis="6-311+g(2df,2p)"):
        "A method to create a calculator that optimizes the reaction shell"

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        assert isinstance(ts, TS), "A TS object was not provided..."

        indicies = []
        for i, atom in enumerate(ts.rmg_molecule.atoms):
            if not (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a + 1, b + 1)

        ts.rmg_molecule.updateMultiplicity()

        label = ts.reaction_label.replace(
            "(", "left").replace(")", "right") + "_shell_" + str(ts.index)

        new_scratch = os.path.join(
                scratch,
                "ts",
                ts.reaction_label,
                "conformers"
            )

        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)

        calc = ASEGaussian(mem=mem,
                           nprocshared=nprocshared,
                           label=label,
                           scratch=new_scratch,
                           method=method,
                           basis=basis,
                           extra="Opt=(ModRedun,Loose,maxcycles=900) Int(Grid=SG1)",
                           multiplicity=ts.rmg_molecule.multiplicity,
                           addsec=[combos[:-1]])
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def get_center_calc(self,
                        ts=None,
                        mem="5GB",
                        nprocshared=20,
                        scratch=".",
                        method="m062x",
                        basis="6-311+g(2df,2p)"):
        "A method to create a calculator that optimizes the reaction shell"

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        assert isinstance(ts, TS), "A TS object was not provided..."

        indicies = []
        for i, atom in enumerate(ts.rmg_molecule.atoms):
            if not (atom.label != ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a + 1, b + 1)

        ts.rmg_molecule.updateMultiplicity()

        label = ts.reaction_label.replace(
            "(", "left").replace(")", "right") + "_center_" + str(ts.index)

        new_scratch = os.path.join(
                scratch,
                "ts",
                ts.reaction_label,
                "conformers"
            )

        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)

        calc = ASEGaussian(mem=mem,
                           nprocshared=nprocshared,
                           label=label,
                           scratch=new_scratch,
                           method=method,
                           basis=basis,
                           extra="Opt=(ModRedun,Loose,maxcycles=900) Int(Grid=SG1)",
                           multiplicity=ts.rmg_molecule.multiplicity,
                           addsec=[combos[:-1]])
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def get_overall_calc(self,
                         ts=None,
                         mem="5GB",
                         nprocshared=20,
                         scratch=".",
                         method="m062x",
                         basis="6-311+g(2df,2p)"):
        "A method to create a calculator that optimizes the reaction shell"

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        assert isinstance(ts, TS), "A TS object was not provided..."

        ts.rmg_molecule.updateMultiplicity()

        label = ts.reaction_label.replace(
            "(", "left").replace(")", "right") + "_" + str(ts.index)

        new_scratch = os.path.join(
                scratch,
                "ts",
                ts.reaction_label,
                "conformers"
            )

        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)

        calc = ASEGaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=method,
            basis=basis,
            extra="opt=(ts,calcfc,noeigentest,maxcycles=900) freq",
            multiplicity=ts.rmg_molecule.multiplicity)
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def get_irc_calc(self,
                     ts=None,
                     mem="5GB",
                     nprocshared=20,
                     scratch=".",
                     method="m062x",
                     basis="6-311+g(2df,2p)"):
        "A method to create the IRC calculator object"

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        ts.rmg_molecule.updateMultiplicity()
        label = ts.reaction_label.replace(
            "(", "left").replace(")", "right") + "_irc_" + str(ts.index)

        new_scratch = os.path.join(
                scratch,
                "ts",
                ts.reaction_label,
                "irc"
            )
        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)

        calc = ASEGaussian(mem=mem,
                           nprocshared=nprocshared,
                           label=label,
                           scratch=new_scratch,
                           method=method,
                           basis=basis,
                           extra="irc=(calcall)",
                           multiplicity=ts.rmg_molecule.multiplicity)
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def calculate(self, conformer=None, calc=None):
        """
        A method to perform a calculation given a calculator and an AutoTST
        object. If the corresponding log file already exists, we will skip it

        :params:
        autotst_object: (Molecule, TS, Reaction) an
        AutoTST object that you want to run calculations on
        calc: (ase.calculators.calculator) the calculator that you want to run
        """

        assert conformer, "A Conformer or TS object needs to be provided to run calculate..."
        assert calc, "An ASECalculator object must be provided to run calculate..."

        current_path = os.getcwd()
        scratch_path = os.path.expanduser(
            calc.scratch)

        new_file_name = calc.label.replace(
            "left", "(").replace("right", ")") + ".log"
        old_file_name = calc.label + ".log"

        os.chdir(scratch_path)

        if os.path.exists(new_file_name):
            # We found a finished file file... it should be fixed
            logging.info(
                "Found previous file for {}, verifying it...".format(new_file_name))
            complete, success = self.verify_output_file(new_file_name)
            if success:
                logging.info("Old output file verified, reading it in...")
                conformer.ase_molecule = read_gaussian_out(new_file_name)
                conformer.energy = conformer.ase_molecule.get_potential_energy()
                conformer.update_coords()
                os.chdir(current_path)
                return conformer, True

            elif complete:
                logging.info(
                    "Output file did not converge, attempting to run one last time...")
                try:
                    calc.calculate(conformer.ase_molecule)
                    conformer.ase_molecule = read_gaussian_out(
                        old_file_name)
                    conformer.energy = conformer.ase_molecule.get_potential_energy()
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True

                except BaseException:  # TODO: add error for seg fault
                    logging.info("{} failed... again...".format(new_file_name))
                    os.chdir(current_path)
                    return conformer, False

            elif (new_file_name == old_file_name) and (not complete):
                # The file names are identical and the job isn't complete yet

                logging.info(
                    "Job appears to be running for this calculation, waiting for it to complete...")

                from time import sleep

                f = open(old_file_name)
                lines = f.readlines()[:5]
                num = ""
                for line in lines:
                    if "Entering Link" in line:
                        num = line.split()[-1][:-1]
                scratch_file = "Gau-" + num + ".int"
                while os.path.exists(scratch_file):
                    sleep(60)
                logging.info(
                    "Job complete, reading in results now by running calculate again...")

                # waiting a lil while to make sure that the file is fixed...
                # just in case...
                sleep(15)
                try:
                    conformer.ase_molecule = read_gaussian_out(
                        old_file_name)
                    conformer.energy = conformer.ase_molecule.get_potential_energy()
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True
                except IndexError:
                    logging.info(
                        "It appears that the previous log file wasn't finished... removing the files and rerunning")
                    os.remove(old_file_name)
                    os.remove(old_file_name.replace(".log", ".ase"))
                    os.remove(old_file_name.replace(".log", ".com"))
                    return self.calculate(conformer, calc)

            else:
                logging.info(
                    "Something went wrong... File is neither complete nor successful...")

                return conformer, False

        elif os.path.exists(old_file_name):
            complete, success = self.verify_output_file(old_file_name)

            if not complete:
                logging.info(
                    "Job appears to be running already, waiting for it to complete...")

                from time import sleep

                f = open(old_file_name)
                lines = f.readlines()[:5]
                num = None
                for line in lines:
                    if "Entering Link" in line:
                        num = line.split()[-1][:-1]

                if not num:
                    logging.info(
                        "Something is wrong... it seems this run was interupted...")
                    logging.info(
                        "deleting {} and recalculating...".format(old_file_name))
                    os.remove(old_file_name)
                    return self.calculate(conformer, calc)

                scratch_file = "Gau-" + num + ".int"
                while os.path.exists(scratch_file):
                    sleep(60)
                logging.info(
                    "Job complete, reading in results now by running calculate again...")

                # waiting a lil while to make sure that the file is fixed...
                # just in case...
                sleep(30)

                return self.calculate(conformer, calc)

            else:
                logging.info(
                    "Found previous file for {}, verifying it...".format(old_file_name))
                if success:
                    logging.info("Old output file verified, reading it in...")
                    conformer.ase_molecule = read_gaussian_out(old_file_name)
                    conformer.energy = conformer.ase_molecule.get_potential_energy()
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True
                else:
                    logging.info(
                        "Something went wrong... File is neither complete nor successful...")

                    return conformer, False
        # Seeing if the file exists
        else:
            # File doesn't exist, running calculations
            logging.info(
                "Starting calculation for {}...".format(new_file_name))
            try:
                calc.calculate(conformer.ase_molecule)
                conformer.ase_molecule = read_gaussian_out(old_file_name)
                conformer.energy = conformer.ase_molecule.get_potential_energy()
                conformer.update_coords()
                os.chdir(current_path)
                return conformer, True
            except BaseException:  # TODO: add error for seg fault
                # first calc failed, trying it once more
                logging.info(
                    "Failed first attempt for {}. Trying it once more...".format(new_file_name))
                try:
                    calc.calculate(conformer.ase_molecule)
                    conformer.ase_molecule = read_gaussian_out(old_file_name)
                    conformer.energy = conformer.ase_molecule.get_potential_energy()
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True
                except BaseException:  # TODO: add error for seg fault
                    logging.info(
                        "{} failed first and second attempt...".format(new_file_name))
                    os.chdir(current_path)
                    return conformer, False

    def verify_output_file(self, path):
        """
        A method to verify output files and make sure that they successfully converged, if not, re-running them

        Returns a tuple where the first entry indicates if the file is complete, the second indicates if it was successful
        """

        if not os.path.exists(path):
            logging.info("Not a valid path, cannot be verified...")
            return (False, False)

        f = open(path, "r")
        file_lines = f.readlines()[-5:]
        verified = (False, False)
        for file_line in file_lines:
            if " Normal termination" in file_line:
                verified = (True, True)
            if " Error termination" in file_line:
                verified = (True, False)

        return verified

    def run_irc(self, conformer=None, calc=None):
        "A method to run the IRC calculation"

        assert "irc" in calc.label, "The calculator provided is not an IRC calculator"
        logging.info("Running IRC calculation")

        current_path = os.getcwd()
        scratch_path = os.path.expanduser(
            calc.scratch)

        new_file_name = calc.label.replace(
            "left", "(").replace("right", ")") + ".log"
        old_file_name = self.irc_calc.label + ".log"

        os.chdir(scratch_path)
        if os.path.exists(new_file_name):
            logging.info(
                "It seems that an old IRC has been run, seeing if it's complete...")
            complete, converged = self.verify_output_file(new_file_name)
            if complete and converged:
                logging.info(
                    "Previous IRC complete and resulted in Normal Termination, verifying it...")
                os.chdir(current_path)

            else:
                logging.info(
                    "Previous IRC was not successful or incomplete... Rerunning it...")
                try:
                    calc.calculate(conformer.ase_molecule)
                except BaseException:
                    # This normally fails because of an issue with ase's
                    # `read_results` method.
                    os.chdir(current_path)
                    pass
                logging.info("IRC calc complete!")
        else:
            logging.info(
                "No previous IRC clac has been run, starting a new one...")
            try:
                calc.calculate(conformer.ase_molecule)
            except BaseException:
                # This normally fails because of an issue with ase's
                # `read_results` method.
                os.chdir(current_path)
                pass
            logging.info("IRC calc complete!")

    def validate_irc(self, calc=None):
        """
        A method to verify an IRC calc
        """
        assert "irc" in calc.label, "The calculator provided is not an IRC calculator"

        reaction_label = calc.label.split("_irc")[0]

        reaction_label = reaction_label.replace("left", "(").replace("right", ")")

        logging.info("Validating IRC file...")
        irc_path = os.path.join(calc.scratch,
                                calc.label + ".log")
        if not os.path.exists(irc_path):
            logging.info(
                "It seems that the file was `fixed`, reading in the `fixed` version.")
            irc_path = irc_path.replace("left", "(").replace("right", ")")

            if not os.path.exists(irc_path):
                logging.info(
                    "It seems that the IRC claculation has not been run.")
                return False

        f = open(irc_path, "r")
        file_lines = f.readlines()[-5:]

        completed = False
        for file_line in file_lines:
            if "Normal termination" in file_line:
                logging.info("IRC successfully ran")
                completed = True
        if not completed:
            logging.info("IRC failed... could not be validated...")
            return False

        pth1 = list()
        steps = list()
        with open(irc_path) as outputFile:
            for line in outputFile:
                line = line.strip()

                if line.startswith('Point Number:'):
                    if int(line.split()[2]) > 0:
                        if int(line.split()[-1]) == 1:
                            ptNum = int(line.split()[2])
                            pth1.append(ptNum)
                        else:
                            pass
                elif line.startswith('# OF STEPS ='):
                    numStp = int(line.split()[-1])
                    steps.append(numStp)
        # This indexes the coordinate to be used from the parsing
        if steps == []:
            logging.error('No steps taken in the IRC calculation!')
            return False
        else:
            pth1End = sum(steps[:pth1[-1]])
            # Compare the reactants and products
            ircParse = ccread(irc_path)
            # cf.
            # http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information

            atomcoords = ircParse.atomcoords
            atomnos = ircParse.atomnos
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = RMGMolecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = RMGMolecule()
            mol2.fromXYZ(atomnos, atomcoords[-1])

            testReaction = RMGReaction(
                reactants=mol1.split(),
                products=mol2.split(),
            )

            r, p = reaction_label.split("_")

            reactants = []
            products = []

            for react in r.split("+"):
                react = RMGMolecule(SMILES=react)
                react.toSingleBonds()
                reactants.append(react)

            for prod in p.split("+"):
                prod = RMGMolecule(SMILES=prod)
                prod.toSingleBonds()
                products.append(prod)

            targetReaction = RMGReaction(
                reactants=reactants,
                products=products,
            )

            if targetReaction.isIsomorphic(testReaction):
                return True
            else:
                return False

    def run(self,
            conformer=None,
            vibrational_analysis=True,
            hindered_rotors=True,):
        """
        A method to perform all the necessary calculations required for a particular conformer
        """

        if not conformer:
            conformer = self.conformer

        assert isinstance(conformer, (Conformer, TS)
                          ), "`conformer` provided not a Conformer type..."

        if isinstance(conformer, TS):
            # Performing the TS optimizations
            logging.info("Conformer provided is a TS object")

            shell = self.get_shell_calc(conformer)
            logging.info("Running optimization of reaction shell")
            conformer, result = self.calculate(conformer, shell)
            self.fix_io_file(shell)
            if not result:
                logging.info("FAILED SHELL CALCULATION")
                return result

            center = self.get_center_calc(conformer)
            logging.info("Running optization of reaction center")
            conformer, result = self.calculate(conformer, center)
            self.fix_io_file(center)
            if not result:
                logging.info("FAILED CENTER CALCULATION")
                return result

            overall = self.get_overall_calc(conformer)
            logging.info("Running overall optimization of TS")
            conformer, result = self.calculate(conformer, overall)
            self.fix_io_file(overall)
            if not result:
                logging.info("FAILED OVERALL CALCULATION")

            if not vibrational_analysis:
                logging.info(
                    "Running without vibrational analysis. \nRunning IRC instead")
                irc = self.get_irc_calc(conformer)
                self.run_irc(conformer, irc)
                result = self.validate_irc(irc)
                self.fix_io_file(irc)

            else:
                from autotst.calculators.vibrational_analysis import VibrationalAnalysis
                vib = VibrationalAnalysis(ts=conformer, scratch=self.scratch)
                result = vib.validate_ts()

                if not result:
                    logging.info(
                        "Vibrational Analysis not conclusive...\n Running IRC instead")
                    irc = self.get_irc_calc(conformer)
                    self.run_irc(conformer, irc)
                    result = self.validate_irc(irc)
                    self.fix_io_file(irc)

            if result:
                logging.info(
                    "TS validated, now running hindered rotor calculations")
                # Add hindered rotor work here
                logging.info("jk, this feature hasn't been added just yet")

            if result:
                logging.info("Arrived at a TS!")
                return result

            else:
                logging.info("Could not arrive at a TS!")
                return result

        elif isinstance(conformer, Conformer):
            logging.info("Conformer provided is NOT a TS object")

            calc = self.get_species_calc(conformer)
            conformer, result = self.calculate(conformer, calc)
            self.fix_io_file(calc)

            if result:
                logging.info(
                    "TS validated, now running hindered rotor calculations")
                # Add hindered rotor work here
                logging.info("jk, this feature hasn't been added just yet")

            if result:
                logging.info("Conformer species successfully optimized")
                return result

            else:
                logging.info("Could not optimize species geometry")
                return result


"""
def calculate_rotor(self, conformer, calculator):

    A method to run hindered rotor calculations


    try:
        calc.calculate(conformer.ase_molecule)
    except:
        pass

    path = os.path.join(calc.scratch, calc.label + ".log")

    if not (self.verify_rotor(path) and self.verify_output_file(path)):
        logging.info(
            "Could not verify the rotor, this file will not be included in calculations.")
        logging.info("File {} renamed as {}...".format(
            path, path.replace(".log", "-failed.log")))
        os.rename(path, path.replace(".log", "-failed.log"))

def verify_rotor(self, path):
    "This could be extrapolated to the general calculators class...?"

    parser = ccread(path)

    smallest = max(parser.scfenergies) + 1
    results = []
    for i in parser.scfenergies:
        if i < smallest:
            smallest = i
        else:
            results.append(smallest)
            smallest = max(parser.scfenergies) + 1
    # adding the last one which should be a converged geometry
    results.append(smallest)

    if ((results[0] - results[-1] < 1e-5) and # The energy difference is less than 1e-5 eV
        (((parser.converged_geometries[0] - parser.converged_geometries[i]) ** 2).mean() < 0.01)): # the RMSE between initial and final geometries is less than 1%
        return True

    else:
        return False"""
