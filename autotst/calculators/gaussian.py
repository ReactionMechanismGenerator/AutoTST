#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
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
################################################################################

import os
import itertools
import logging
import numpy as np

import rmgpy

import autotst
from autotst.reaction import AutoTST_Reaction, AutoTST_TS
from autotst.molecule import AutoTST_Molecule
from autotst.calculators.vibrational_analysis import Vibrational_Analysis
from autotst.calculators.calculator import AutoTST_Calculator

from rdkit import Chem
from cclib.io import ccread

from ase.io.gaussian import read_gaussian, read_gaussian_out
from ase.calculators.gaussian import Gaussian


def update_from_ase(autotst_obj, ase_object):
    """
    A function designed to update all objects based off of their ase objects
    """
    if isinstance(autotst_obj, autotst.molecule.AutoTST_Molecule):
        autotst_obj.ase_molecule = ase_object
        autotst_obj.update_from_ase_mol()

    if isinstance(autotst_obj, autotst.reaction.AutoTST_Reaction):
        autotst_obj.ts.ase_ts = ase_object
        autotst_obj.ts.update_from_ase_ts()

    if isinstance(autotst_obj, autotst.reaction.AutoTST_TS):
        autotst_obj.ase_ts = ase_object
        autotst_obj.update_from_ase_ts()

    return autotst_obj


class AutoTST_Gaussian(AutoTST_Calculator):

    def __init__(self,
                 autotst_reaction,
                 mem="5GB",
                 nprocshared=20,
                 scratch=".",
                 method="m062x",
                 basis="6-311+g(2df,2p)",
                 save_directory="."):
        """
        A method to create all of the calculators needed for AutoTST

        :params:
        autotst_reaction: (AutoTST_Reaction) The reaction of interest
        scratch: (str) The directory that you would like to use for calculations
        """

        self.reaction = autotst_reaction
        self.mem = mem
        self.nprocshared = nprocshared
        self.scratch = scratch
        self.method = method
        self.basis = basis
        self.save_directory = save_directory

        self.get_reactant_and_product_calcs(
            self.mem, self.nprocshared, self.scratch, self.method, self.basis)

        self.shell_calc = self.get_shell_calc(
            self.mem, self.nprocshared, self.scratch, self.method, self.basis)
        self.center_calc = self.get_center_calc(
            self.mem, self.nprocshared, self.scratch, self.method, self.basis)
        self.overall_calc = self.get_overall_calc(
            self.mem, self.nprocshared, self.scratch, self.method, self.basis)
        self.irc_calc = self.get_irc_calc(
            self.mem, self.nprocshared, self.scratch, self.method, self.basis)

        self.completed_irc = False

    def __repr__(self):
        return '<AutoTST Gaussian Calculators "{0}">'.format(self.reaction.label)

    def reactants_or_products_calc(self, autotst_mol, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method that creates a calculator for a reactant or product"

        autotst_mol.rmg_molecule.updateMultiplicity()

        # using this round about way of doing stuff because rmg's `toAugumentedInChIKey` method doesn't work on our cluster

        smiles = autotst_mol.rmg_molecule.toSMILES()
        label = Chem.rdinchi.InchiToInchiKey(
            Chem.MolToInchi(Chem.MolFromSmiles(smiles))).strip("-N")

        calc = Gaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="opt=(verytight,gdiis) freq IOP(2/16=3)",
                        multiplicity=autotst_mol.rmg_molecule.multiplicity
                        )
        del calc.parameters['force']

        return calc

    def get_rotor_calcs(self, autotst_object, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        """ A method to create all of the calculators needed to perform hindered rotor calculations"""

        calculators = {}
        for torsion in autotst_object.torsions:
            string = ""
            for bond in mol.bonds:
                i, j = bond.indices
                string += "B {} {} F\n".format(i+1, j+1)

            for angle in mol.angles:
                i, j, k = angle.indices
                string += "A {} {} {} F\n".format(i+1, j+1, k+1)

            i, j, k, l = torsion.indices
            string += "D {} {} {} {} S 36 10.0".format(i+1, j+1, k+1, l+1)

            smiles = autotst_object.rmg_molecule.toSMILES()

            label = Chem.rdinchi.InchiToInchiKey(
                Chem.MolToInchi(Chem.MolFromSmiles(smiles))).strip("-N")
            label += "_tor{}{}".format(j, k)
            calc = Gaussian(mem=mem,
                            nprocshared=nprocshared,
                            label=label,
                            scratch=scratch,
                            method=method,
                            basis=basis,
                            extra="Opt=(CalcFC,ModRedun)",
                            multiplicity=autotst_object.rmg_molecule.multiplicity,
                            addsec=[string])

            del calc.parameters['force']
            calculators[(i, j)] = calc

        return calculators

    def get_reactant_and_product_calcs(self, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method that collects all of the calculators for reactants and prods"

        self.reactant_calcs = {}
        self.product_calcs = {}

        for reactant in self.reaction.reactant_mols:
            calc = self.reactants_or_products_calc(
                reactant, mem, nprocshared, scratch, method, basis)
            self.reactant_calcs[reactant] = calc

        for product in self.reaction.product_mols:
            calc = self.reactants_or_products_calc(
                product, mem, nprocshared, scratch, method, basis)
            self.product_calcs[product] = calc

    def get_shell_calc(self, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create a calculator that optimizes the reaction shell"

        indicies = []
        for i, atom in enumerate(self.reaction.ts.rmg_ts.atoms):
            if not (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a+1, b+1)

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace(
            "(", "left").replace(")", "right") + "_shell"

        calc = Gaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="Opt=(ModRedun,Loose) Int(Grid=SG1)",
                        multiplicity=self.reaction.ts.rmg_ts.multiplicity,
                        addsec=[combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_center_calc(self, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the calculator to perform the reaction center opt"

        indicies = []
        for i, atom in enumerate(self.reaction.ts.rmg_ts.atoms):
            if (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a+1, b+1)

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace(
            "(", "left").replace(")", "right") + "_center"

        calc = Gaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="opt=(ts,calcfc,noeigentest,ModRedun)",
                        multiplicity=self.reaction.ts.rmg_ts.multiplicity,
                        addsec=[combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_overall_calc(self, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the calculator to perform the full TS optimization"

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace("(", "left").replace(")", "right")

        calc = Gaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="opt=(ts,calcfc,noeigentest) freq",
                        multiplicity=self.reaction.ts.rmg_ts.multiplicity)

        del calc.parameters['force']
        return calc

    def get_irc_calc(self, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the IRC calculator object"

        self.reaction.ts.rmg_ts.updateMultiplicity()
        label = self.reaction.label.replace(
            "(", "left").replace(")", "right") + "_irc"

        calc = Gaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="irc=(calcall)",
                        multiplicity=self.reaction.ts.rmg_ts.multiplicity)

        del calc.parameters['force']
        return calc

    def calculate(self, autotst_object, calc):
        """
        A method to perform a calculation given a calculator and an AutoTST
        object. If the corresponding log file already exists, we will skip it

        :params:
        autotst_object: (AutoTST_Molecule, AutoTST_TS, AutoTST_Reaction) an
        AutoTST object that you want to run calculations on
        calc: (ase.calculators.calculator) the calculator that you want to run
        """

        def update_from_ase(autotst_obj, ase_object):
            """
            A function designed to update all objects based off of their ase objects
            """
            if isinstance(autotst_obj, autotst.molecule.AutoTST_Molecule):
                autotst_obj.ase_molecule = ase_object
                autotst_obj.update_from_ase_mol()

            if isinstance(autotst_obj, autotst.reaction.AutoTST_Reaction):
                autotst_obj.ts.ase_ts = ase_object
                autotst_obj.ts.update_from_ase_ts()

            if isinstance(autotst_obj, autotst.reaction.AutoTST_TS):
                autotst_obj.ase_ts = ase_object
                autotst_obj.update_from_ase_ts()

            return autotst_obj

        current_path = os.getcwd()
        scratch_path = os.path.expanduser(
            calc.scratch).replace(".", os.getcwd())

        new_file_name = calc.label.replace(
            "left", "(").replace("right", ")") + ".log"
        old_file_name = calc.label + ".log"

        if isinstance(autotst_object, AutoTST_Molecule):
            ase_object = autotst_object.ase_molecule
        elif isinstance(autotst_object, AutoTST_Reaction):
            ase_object = autotst_object.ts.ase_ts
        elif isinstance(autotst_object, AutoTST_TS):
            ase_object = autotst_object.ase_ts

        os.chdir(scratch_path)

        # Seeing if the file exists
        if not os.path.exists(new_file_name):
            # File doesn't exist, running calculations
            logging.info(
                "Starting calculation for {}...".format(new_file_name))
            try:
                calc.calculate(ase_object)
                ase_object = read_gaussian_out(old_file_name)
                autotst_object = update_from_ase(autotst_object, ase_object)
                os.chdir(current_path)
                return autotst_object, True
            except:  # TODO: add error for seg fault
                # first calc failed, trying it once more
                logging.info(
                    "Failed first attempt for {}. Trying it once more...".format(new_file_name))
                try:
                    calc.calculate(ase_object)
                    ase_object = read_gaussian_out(old_file_name)
                    autotst_object = update_from_ase(
                        autotst_object, ase_object)
                    os.chdir(current_path)
                    return autotst_object, True
                except:  # TODO: add error for seg fault
                    logging.info(
                        "{} failed first and second attempt...".format(new_file_name))
                    os.chdir(current_path)
                    return autotst_object, False

        else:
            # We found an old file... it should be fixed
            logging.info(
                "Found previous file for {}, verifying it...".format(new_file_name))
            if self.verify_output_file(new_file_name):
                logging.info("Old output file verified, reading it in...")
                ase_object = read_gaussian_out(new_file_name)
                autotst_object = update_from_ase(autotst_object, ase_object)
                os.chdir(current_path)
                return autotst_object, True

            else:
                logging.info(
                    "Could not verify output file, attempting to run one last time...")
                try:
                    calc.calculate(ase_object)
                    autotst_object.ase_molecule = read_gaussian_out(
                        old_file_name)
                    autotst_object = update_from_ase(
                        autotst_object, ase_object)
                    os.chdir(current_path)
                    return autotst_object, True

                except:  # TODO: add error for seg fault
                    logging.info("{} failed... again...".format(new_file_name))
                    os.chdir(current_path)
                    return autotst_object, False

    def verify_output_file(self, path):
        "A method to verify output files and make sure that they successfully converged, if not, re-running them"

        if not os.path.exists(path):
            print "Not a validat path, cannot be verified..."
            return False

        f = open(path, "r")
        file_lines = f.readlines()[-5:]
        verified = False
        for file_line in file_lines:
            if " Normal termination" in file_line:
                verified = True

        return verified

    def run_rotors(self, calculators, autotst_object):

        assert len(calculators) == len(
            autotst_object.torsions), "Incorrectly matched calculators to molecule..."

        if isinstance(autotst_object, AutoTST_Molecule):
            ase_object = autotst_object.ase_molecule
        elif isinstance(autotst_object, AutoTST_Reaction):
            ase_object = autotst_object.ts.ase_ts
        elif isinstance(autotst_object, AutoTST_TS):
            ase_object = autotst_object.ase_ts

        for torsion in autotst_object.torsions:
            i, j, k, l = torsion.indices
            calc = calculators[(j, k)]
            try:
                calc.calculate(ase_object)
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
            return False

    def run_reactants_and_products(self):
        "A method to run the calculations for all reactants and products"

        bools = []
        for mol, calc in self.reactant_calcs.iteritems():
            mol, b = self.calculate(mol, calc)
            self.fix_io_file(calc)
            bools.append(b)

        for mol, calc in self.product_calcs.iteritems():
            mol, b = self.calculate(mol, calc)
            self.fix_io_file(calc)
            bools.append(b)

        return np.array(bools).all()

    def run_shell(self):
        "A method to run the shell optimization with the reaction center frozen"
        logging.info("Running shell optimization with center frozen...")
        self.reaction, bool = self.calculate(self.reaction, self.shell_calc)
        logging.info("Shell optimization complete!")
        return bool

    def run_center(self):
        "A method to run the reaction center optimization with the shell frozen"
        logging.info("Running center optimization with shell frozen...")
        self.reaction, bool = self.calculate(self.reaction, self.center_calc)
        logging.info("Center optimization complete!")
        return bool

    def run_overall(self):
        "A method to run the optimization of the entire TS"
        logging.info("Running overall optimization...")
        self.reaction, bool = self.calculate(self.reaction, self.overall_calc)
        logging.info("Overall optimization complete!")
        return bool

    def run_irc(self):
        "A method to run the IRC calculation"
        logging.info("Running IRC calculation")
        current_dir = os.getcwd
        try:
            scratch_path = os.path.expanduser(
                calc.scratch).replace(".", os.getcwd())
            os.chdir(scratch_path)
            self.irc_calc.calculate(self.reaction.ts.ase_ts)
        except:
            # This normally fails because of an issue with ase's `read_results` method.
            os.chdir(current_dir)
            pass
        logging.info("IRC calc complete!")

    def validate_irc(self):  # TODO: need to add more verification here
        logging.info("Validating IRC file...")
        irc_path = os.path.join(self.irc_calc.scratch,
                                self.irc_calc.label + ".log")
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
            if " Normal termination" in file_line:
                logging.info("IRC successfully ran")
                completed = True
        if completed == False:
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
            # cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            ircParse.logger.setLevel(logging.ERROR)

            atomcoords = ircParse.atomcoords
            atomnos = ircParse.atomnos
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = Molecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = Molecule()
            mol2.fromXYZ(atomnos, atomcoords[-1])

            testReaction = Reaction(
                reactants=mol1.split(),
                products=mol2.split(),
            )

            if isinstance(self.reaction.rmg_reaction.reactants[0], rmgpy.molecule.Molecule):
                targetReaction = Reaction(
                    reactants=[reactant.toSingleBonds()
                               for reactant in self.reaction.rmg_reaction.reactants],
                    products=[product.toSingleBonds()
                              for product in self.reaction.rmg_reaction.products],
                )
            elif isinstance(self.reaction.rmg_reaction.reactants[0], rmgpy.species.Species):
                targetReaction = Reaction(
                    reactants=[reactant.molecule[0].toSingleBonds()
                               for reactant in self.reaction.rmg_reaction.reactants],
                    products=[product.molecule[0].toSingleBonds()
                              for product in self.reaction.rmg_reaction.products],
                )

            if targetReaction.isIsomorphic(testReaction):
                return True
            else:
                return False

    def run_all(self, vibrational_analysis=True):
        """
        A method that is designed to run all of the automated quantum
        calculations for AutoTST. These can be run independently as well.

        :params:
        vibrational_analysis: (bool) A bool to tell AutoTST if you want to use
        vibrational analysis instead of IRC calcs to speed up calculations

        :returns:
        result: (bool) A bool to tell you if an AutoTST run successfully
        converged on a verified TS.
        """
        result = False
        r_and_p = self.run_reactants_and_products()
        if not r_and_p:
            return result
        shell = self.run_shell()
        self.fix_io_file(self.shell_calc)
        if not shell:
            return result
        center = self.run_center()
        self.fix_io_file(self.center_calc)
        if not center:
            return result
        overall = self.run_overall()
        self.fix_io_file(self.overall_calc)
        if not overall:
            return result

        vib = Vibrational_Analysis(
            reaction=self.reaction, scratch=self.scratch)
        logging.info("Performing Vibrational Analysis...")
        if vibrational_analysis and vib.validate_ts():
            logging.info(
                "Vibrational analysis successful! Successfully arrived at a TS.")
            result = True
        elif vibrational_analysis and not vib.validate_ts():
            logging.info(
                "Could not validate via vibrational analysis... \nRunning IRC instead...")
            self.run_irc()
            result = self.validate_irc()
        else:
            logging.info(
                "Running without vibrational analysis... \nRunning IRC instead...")
            self.run_irc()
            result = self.validate_irc()

        self.fix_io_file(self.irc_calc)

        if result:
            logging.info("Arrived at a TS!")
            return result

        else:
            logging.info("Could not arrive at a TS!")
            return result

    def fix_io_file(self, calc):
        """
        A method that removes the `left` and `right` text from a log, ase, and
        com files and turns it back into a smiles structure
        """
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

    def fix_io_files(self):
        """
        A method that removes the `left` and `right` text from a log, ase, and
        com files and turns it back into a smiles structure for ALL files.
        """
        for calc in self.reactant_calcs.values():
            self.fix_io_file(calc)

        for calc in self.product_calcs.values():
            self.fix_io_file(calc)

        self.fix_io_file(self.shell_calc)
        self.fix_io_file(self.center_calc)
        self.fix_io_file(self.overall_calc)
