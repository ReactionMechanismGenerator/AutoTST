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
from cclib.io import ccread

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.reaction import Reaction

import autotst
from autotst.reaction import Reaction, TS
from autotst.species import Species
from autotst.calculators.vibrational_analysis import VibrationalAnalysis
from autotst.calculators.calculator import Calculator

from rdkit import Chem
from cclib.io import ccread

from ase.io.gaussian import read_gaussian, read_gaussian_out
from ase.calculators.gaussian import Gaussian as ASEGaussian



class Gaussian(Calculator):

    def __init__(self,
                 reaction=None,
                 mem="5GB",
                 nprocshared=20,
                 scratch=".",
                 method="m062x",
                 basis="6-311+g(2df,2p)",
                 save_directory="."):
        """
        A method to create all of the calculators needed for AutoTST

        :params:
        autotst_reaction: (Reaction) The reaction of interest
        scratch: (str) The directory that you would like to use for calculations
        """

        self.reaction = reaction
        self.mem = mem
        self.nprocshared = nprocshared
        self.scratch = scratch
        self.method = method
        self.basis = basis
        self.save_directory = save_directory
        

        if reaction:
            self.label = reaction.label

            self.get_species_calcs(
                self.reaction, self.mem, self.nprocshared, self.scratch, self.method, self.basis)

            self.shell_calc = self.get_shell_calc(
                self.reaction, self.reaction.ts["forward"][0], self.mem, self.nprocshared, self.scratch, self.method, self.basis)
            self.center_calc = self.get_center_calc(
                self.reaction, self.reaction.ts["forward"][0], self.mem, self.nprocshared, self.scratch, self.method, self.basis)
            self.overall_calc = self.get_overall_calc(
                self.reaction, self.reaction.ts["forward"][0], self.mem, self.nprocshared, self.scratch, self.method, self.basis)
            self.irc_calc = self.get_irc_calc(
                self.reaction, self.reaction.ts["forward"][0], self.mem, self.nprocshared, self.scratch, self.method, self.basis)

            self.completed_irc = False

    def __repr__(self):
        return '<Gaussian Calculator "{0}">'.format(self.reaction.label)

    def get_species_calc(self, conformer, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method that creates a calculator for a reactant or product"

        conformer.rmg_molecule.updateMultiplicity()

        # using this round about way of doing stuff because rmg's `toAugumentedInChIKey` method doesn't work on our cluster

        smiles = conformer.rmg_molecule.toSMILES()
        label = Chem.rdinchi.InchiToInchiKey(
            Chem.MolToInchi(Chem.MolFromSmiles(smiles))).strip("-N")

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="opt=(verytight,gdiis,maxcycle=1000) freq IOP(2/16=3)",
                        multiplicity=conformer.rmg_molecule.multiplicity
                        )
        del calc.parameters['force']

        return calc

    def get_rotor_calc(self, conformer, torsion, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        """ A method to create all of the calculators needed to perform hindered rotor calculations"""
        
        string = ""
        for bond in conformer.bonds:
            i, j = bond.indices
            string += "B {} {}\n".format(i+1, j+1)

        i, j, k, l = torsion.atom_indices
        string += "D {} {} {} {} S 36 10.0".format(i+1, j+1, k+1, l+1)
        
        if isinstance(conformer, TS):
            label = self.label + "_tor_{}_{}".format(j,k)
        elif isisntance(conformer, Conformer):
            label = conformer.smiles + "_tor_{}_{}".format(j,k)
        
        label = conformer.smiles + "_tor_{}_{}".format(j,k)
        conformer.rmg_molecule.updateMultiplicity()
        mult = conformer.rmg_molecule.multiplicity

        if isinstance(autotst_object, Reaction):
            label = autotst_object.label + "_tor_{}_{}".format(j, k)
            mult = autotst_object.ts.rmg_ts.multiplicity
        elif isinstance(autotst_object, TS):
            label = autotst_object.label + "_tor_{}_{}".format(j, k)
            mult = autotst_object.rmg_ts.multiplicity
        elif isinstance(autotst_object, Species):
            smiles = autotst_object.rmg_molecule.toSMILES()
            label = Chem.rdinchi.InchiToInchiKey(
                Chem.MolToInchi(Chem.MolFromSmiles(smiles))).strip("-N")
            label += "_tor{}{}".format(j, k)
            mult = autotst_object.rmg_molecule.multiplicity

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="Opt=(CalcFC,ModRedun)",
                        multiplicity=mult,
                        addsec=[string])

        del calc.parameters['force']

        return calc

    def get_species_calcs(self, reaction=None, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method that collects all of the calculators for reactants and prods"

        if not reaction:
            reaction = self.reaction

        reactant_calcs = {}
        product_calcs = {}

        for reactant in self.reaction.reactants:
            r_calcs = {}
            for smiles, confs in reactant.conformers.iteritems():

                conf = confs[0]

                calc = self.get_species_calc(
                    conf, mem, nprocshared, scratch, method, basis)


                r_calcs[smiles] = calc
            reactant_calcs[reactant] = r_calcs

        self.reactant_calcs = reactant_calcs

        for product in self.reaction.products:
            p_calcs = {}
            for smiles, confs in product.conformers.iteritems():

                conf = confs[0]

                calc = self.get_species_calc(
                    conf, mem, nprocshared, scratch, method, basis)


                p_calcs[smiles] = calc
            product_calcs[product] = p_calcs

        self.product_calcs = product_calcs

        return reactant_calcs, product_calcs

    def get_shell_calc(self, reaction=None, conformer=None, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create a calculator that optimizes the reaction shell"

        if reaction is None:
            reaction = self.reaction

        if conformer is None:
            logging.info("No TS was provided, selecting the first TS geometry")
            conformer = reaction.ts["forward"][0]

        indicies = []
        for i, atom in enumerate(conformer.rmg_molecule.atoms):
            if not (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a+1, b+1)

        conformer.rmg_molecule.updateMultiplicity()

        label = reaction.label.replace(
            "(", "left").replace(")", "right") + "_shell"

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="Opt=(ModRedun,Loose,maxcycle=1000) Int(Grid=SG1)",
                        multiplicity=conformer.rmg_molecule.multiplicity,
                        addsec=[combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_center_calc(self, reaction=None, conformer=None, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the calculator to perform the reaction center opt"


        if reaction is None:
            reaction = self.reaction

        if conformer is None:
            logging.info("No TS was provided, selecting the first TS geometry")
            conformer = reaction.ts["forward"][0]

        indicies = []
        for i, atom in enumerate(conformer.rmg_molecule.atoms):
            if not (atom.label != ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            combos += "{0} {1} F\n".format(a+1, b+1)

        conformer.rmg_molecule.updateMultiplicity()

        label = reaction.label.replace(
            "(", "left").replace(")", "right") + "_center"

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="Opt=(ModRedun,Loose,maxcycle=1000) Int(Grid=SG1)",
                        multiplicity=conformer.rmg_molecule.multiplicity,
                        addsec=[combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_overall_calc(self, reaction=None, conformer=None, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the calculator to perform the full TS optimization"

        if reaction is None:
            reaction = self.reaction

        if conformer is None:
            logging.info("No TS was provided, selecting the first TS geometry")
            conformer = reaction.ts["forward"][0]


        conformer.rmg_molecule.updateMultiplicity()

        label = self.reaction.label.replace("(", "left").replace(")", "right")

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="opt=(ts,calcfc,noeigentest,maxcycle=1000) freq",
                        multiplicity=conformer.rmg_molecule.multiplicity)

        del calc.parameters['force']
        return calc

    def get_irc_calc(self, reaction=None, conformer=None, mem="5GB", nprocshared=20, scratch=".", method="m062x", basis="6-311+g(2df,2p)"):
        "A method to create the IRC calculator object"

        if reaction is None:
            reaction = self.reaction

        if conformer is None:
            logging.info("No TS was provided, selecting the first TS geometry")
            conformer = reaction.ts["forward"][0]

        conformer.rmg_molecule.updateMultiplicity()
        label = reaction.label.replace(
            "(", "left").replace(")", "right") + "_irc"

        calc = ASEGaussian(mem=mem,
                        nprocshared=nprocshared,
                        label=label,
                        scratch=scratch,
                        method=method,
                        basis=basis,
                        extra="irc=(calcall)",
                        multiplicity=conformer.rmg_molecule.multiplicity)

        del calc.parameters['force']
        return calc

    def calculate(self, conformer, calc):
        """
        A method to perform a calculation given a calculator and an AutoTST
        object. If the corresponding log file already exists, we will skip it

        :params:
        autotst_object: (Molecule, TS, Reaction) an
        AutoTST object that you want to run calculations on
        calc: (ase.calculators.calculator) the calculator that you want to run
        """

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
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True

                except:  # TODO: add error for seg fault
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
                for line in lines:
                    if "Entering Link" in line:
                        num = line.split()[-1][:-1]
                scratch_file = "Gau-" + num + ".int"
                while os.path.exists(scratch_file):
                    sleep(60)
                logging.info("Job complete, reading in results now by running calculate again...")

                sleep(30) # waiting a lil while to make sure that the file is fixed... just in case...
                try:
                    conformer.ase_molecule = read_gaussian_out(
                        old_file_name)
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True
                except IndexError:
                    logging.info("It appears that the previous log file wasn't finished... removing the files and rerunning")
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
                    logging.info("Something is wrong... it seems this run was interupted...")
                    logging.info("deleting {} and recalculating...".format(old_file_name))
                    os.remove(old_file_name)
                    return self.calculate(conformer, calc)
                
                scratch_file = "Gau-" + num + ".int"
                while os.path.exists(scratch_file):
                    sleep(60)
                logging.info("Job complete, reading in results now by running calculate again...")

                sleep(30) # waiting a lil while to make sure that the file is fixed... just in case...

                return self.calculate(conformer, calc)
            
            else:
                logging.info(
                    "Found previous file for {}, verifying it...".format(old_file_name))
                if success:
                    logging.info("Old output file verified, reading it in...")
                    conformer.ase_molecule = read_gaussian_out(old_file_name)
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
                conformer.update_coords()
                os.chdir(current_path)
                return conformer, True
            except:  # TODO: add error for seg fault
                # first calc failed, trying it once more
                logging.info(
                    "Failed first attempt for {}. Trying it once more...".format(new_file_name))
                try:
                    calc.calculate(conformer.ase_molecule)
                    conformer.ase_molecule = read_gaussian_out(old_file_name)
                    conformer.update_coords()
                    os.chdir(current_path)
                    return conformer, True
                except:  # TODO: add error for seg fault
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
            print "Not a valid path, cannot be verified..."
            return False

        f = open(path, "r")
        file_lines = f.readlines()[-5:]
        verified = (False, False)
        for file_line in file_lines:
            if " Normal termination" in file_line:
                verified = (True, True)
            if " Error termination" in file_line:
                verified = (True, False)

        return verified

    def calculate_rotor(self, conformer, calculator):
        """
        A method to run hindered rotor calculations
        """

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
            return False

    def calculate_species(self, reaction=None, reactant_calcs=None, product_calcs=None):
        "A method to run the calculations for all reactants and products"

        if reaction is None:
            reaction = self.reaction 
        if reactant_calcs is None:
            reactant_calcs = self.reactant_calcs
        if product_calcs is None:
            product_calcs = self.reactant_calcs

        bools = []
        for mol in reaction.reactants:
            for smiles, confs in mol.conformers.iteritems():
                calc = reactant_calcs[mol][smiles]

                origin_label = calc.label[:]

                for i, conf in enumerate(confs):
                    calc.label = origin_label + "_{}".format(i)
                    conf, b = self.calculate(conf, calc)
                    self.fix_io_file(calc)
                    bools.append(b)

        for mol in reaction.products:
            for smiles, confs in mol.conformers.iteritems():
                calc = product_calcs[mol][smiles]

                origin_label = calc.label[:]

                for i, conf in enumerate(confs):
                    calc.label = origin_label + "_{}".format(i)
                    conf, b = self.calculate(conf, calc)
                    self.fix_io_file(calc)
                    bools.append(b)

        return np.array(bools).all()

    def run_shell(self, conformer=None, shell_calc=None):
        "A method to run the shell optimization with the reaction center frozen"

        if conformer is None:
            conformer = self.reaction.ts["forward"][0]
        if shell_calc is None:
            shell_calc = self.shell_calc
        logging.info("Running shell optimization with center frozen...")
        conformer, bool = self.calculate(conformer, shell_calc)
        logging.info("Shell optimization complete!")
        return conformer, bool

    def run_center(self, conformer=None, center_calc=None):
        "A method to run the reaction center optimization with the shell frozen"
        if conformer is None:
            conformer = self.reaction.ts["forward"][0]
        if center_calc is None:
            center_calc = self.center_calc
        logging.info("Running center optimization with shell frozen...")
        conformer, bool = self.calculate(conformer, center_calc)
        logging.info("Center optimization complete!")
        return conformer, bool

    def run_overall(self, conformer=None, overall_calc=None):
        "A method to run the optimization of the entire TS"
        if conformer is None:
            conformer = self.reaction.ts["forward"][0]
        if overall_calc is None:
            overall_calc = self.overall_calc
        logging.info("Running overall optimization...")
        conformer, bool = self.calculate(conformer, overall_calc)
        logging.info("Overall optimization complete!")
        return conformer, bool

    def run_irc(self):
        "A method to run the IRC calculation"
        logging.info("Running IRC calculation")

        current_path = os.getcwd()
        scratch_path = os.path.expanduser(
            self.irc_calc.scratch)

        new_file_name = self.irc_calc.label.replace(
            "left", "(").replace("right", ")") + ".log"
        old_file_name = self.irc_calc.label + ".log"

        os.chdir(scratch_path)
        if os.path.exists(new_file_name):
            logging.info("It seems that an old IRC has been run, seeing if it's complete...")
            if self.verify_output_file(new_file_name):
                logging.info("Previous IRC complete and resulted in Normal Termination, verifying it...")
                os.chdir(current_path)

            else:
                logging.info("Previous IRC was not successful or incomplete... Rerunning it...")
                try:
                    self.irc_calc.calculate(self.reaction.ts.ase_ts)
                except:
                    # This normally fails because of an issue with ase's `read_results` method.
                    os.chdir(current_path)
                    pass
                logging.info("IRC calc complete!")
        else:
            logging.info("No previous IRC clac has been run, starting a new one...")
            try:
                self.irc_calc.calculate(self.reaction.ts.ase_ts)
            except:
                # This normally fails because of an issue with ase's `read_results` method.
                os.chdir(current_path)
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

            atomcoords = ircParse.atomcoords
            atomnos = ircParse.atomnos
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = RMGMolecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = RMGMolecule()
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
        r_and_p = self.calculate_species()
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
