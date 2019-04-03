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

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.reaction import Reaction as RMGReaction

import autotst
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.calculators.calculator import Calculator
from autotst.geometry import Torsion

from cclib.io import ccread

from ase import Atom, Atoms
from ase.io.gaussian import read_gaussian, read_gaussian_out
from ase.calculators.gaussian import Gaussian as ASEGaussian


class Gaussian(Calculator):

    def __init__(self,
                 conformer=None,
                 mem="5GB",
                 nprocshared=20,
                 scratch=".",
                 method="m062x",
                 basis="cc-pVTZ",
                 save_directory="."):

        self.command = "g16"

        assert isinstance(conformer, (type(None), Conformer)), "Please provide a Conformer object"
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
            label = self.conformer.smiles
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
                       basis="cc-pVTZ",
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
            label = d = conformer.smiles
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

    def get_conformer_calc(self,
                         conformer=None,
                         mem="5GB",
                         nprocshared=20,
                         scratch=".",
                         method="m062x",
                         basis="cc-pVTZ"):
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

        short_label = conformer.smiles
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
            extra="opt=(calcfc,verytight,gdiis,maxcycles=900) freq IOP(2/16=3) scf=(maxcycle=900)",
            multiplicity=conformer.rmg_molecule.multiplicity)
        calc.atoms = conformer.ase_molecule
        del calc.parameters['force']

        return calc

    def get_shell_calc(self,
                       ts=None,
                       direction="forward",
                       mem="5GB",
                       nprocshared=20,
                       scratch=".",
                       method="m062x",
                       basis="cc-pVTZ"):
        "A method to create a calculator that optimizes the reaction shell"

        assert direction.lower() in ["forward", "reverse"]

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        assert isinstance(ts, TS), "A TS object was not provided..."

        ts.rmg_molecule.updateMultiplicity()

        label = ts.reaction_label + "_" + direction.lower() + "_shell_" + str(ts.index)

        new_scratch = os.path.join(
                scratch,
                "ts",
                ts.reaction_label,
                "conformers"
            )

        if not os.path.isdir(new_scratch):
            os.makedirs(new_scratch)


        ind1 = ts.rmg_molecule.getLabeledAtom("*1").sortingLabel
        ind2 = ts.rmg_molecule.getLabeledAtom("*2").sortingLabel
        ind3 = ts.rmg_molecule.getLabeledAtom("*3").sortingLabel

        combos = ""
        combos += "{0} {1} F\n".format(ind1+1, ind2+1)
        combos += "{0} {1} F\n".format(ind2+1, ind3+1)
        combos += "{0} {1} {2} F".format(ind1+1, ind2+1, ind3+1)


        calc = ASEGaussian(mem=mem,
                           nprocshared=nprocshared,
                           label=label,
                           scratch=new_scratch,
                           method=method,
                           basis=basis,
                           extra="Opt=(ModRedun,Loose,maxcycles=900) Int(Grid=SG1) scf=(maxcycle=900)",
                           multiplicity=ts.rmg_molecule.multiplicity,
                           addsec=[combos])
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def get_center_calc(self,
                        ts=None,
                        direction="forward",
                        mem="5GB",
                        nprocshared=20,
                        scratch=".",
                        method="m062x",
                        basis="cc-pVTZ"):
        "A method to create a calculator that optimizes the reaction shell"

        assert direction.lower() in ["forward", "reverse"]

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

        label = ts.reaction_label + "_" + direction.lower() + "_center_" + str(ts.index)

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
                           extra="Opt=(ts,calcfc,noeigentest,ModRedun,maxcycles=900) scf=(maxcycle=900)",
                           multiplicity=ts.rmg_molecule.multiplicity,
                           addsec=[combos[:-1]])
        calc.atoms = ts.ase_molecule
        del calc.parameters['force']

        return calc

    def get_overall_calc(self,
                         ts=None,
                         direction="forward",
                         mem="5GB",
                         nprocshared=20,
                         scratch=".",
                         method="m062x",
                         basis="cc-pVTZ"):
        "A method to create a calculator that optimizes the reaction shell"

        assert direction.lower() in ["forward", "reverse"]

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        assert isinstance(ts, TS), "A TS object was not provided..."

        ts.rmg_molecule.updateMultiplicity()

        label = ts.reaction_label + "_" + direction.lower() + "_" + str(ts.index)

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
            extra="opt=(ts,calcfc,noeigentest,maxcycles=900) freq scf=(maxcycle=900)",
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
                     basis="cc-pVTZ"):
        "A method to create the IRC calculator object"

        if ts is None:
            if self.ts is None:
                return None
            elif not isinstance(self.conformer, TS):
                return None
            else:
                ts = self.conformer

        ts.rmg_molecule.updateMultiplicity()
        label = ts.reaction_label + "_irc_" + str(ts.index)

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
        irc_path = os.path.join(
            calc.scratch,
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
                reactants.append(react)

            for prod in p.split("+"):
                prod = RMGMolecule(SMILES=prod)
                products.append(prod)

            possible_reactants = []
            possible_products = []
            for reactant in reactants:
                possible_reactants.append(reactant.generate_resonance_structures())
                
            for product in products:
                possible_products.append(product.generate_resonance_structures())
                
            possible_reactants = list(itertools.product(*possible_reactants))
            possible_products = list(itertools.product(*possible_products))

            for possible_reactant in possible_reactants:
                reactant_list = []
                for react in possible_reactant:
                    reactant_list.append(react.toSingleBonds())
                    
                for possible_product in possible_products:
                    product_list = []
                    for prod in possible_product:
                        product_list.append(prod.toSingleBonds())
                    
                    targetReaction = RMGReaction(
                        reactants = list(reactant_list),
                        products = list(product_list)
                    )

                    if targetReaction.isIsomorphic(testReaction):
                        logging.info("IRC calculation was successful!")
                        return True
            logging.info("IRC calculation failed :(")
            return False


