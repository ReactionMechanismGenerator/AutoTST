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


import arkane
from arkane import Arkane as RMGArkane, KineticsJob, StatMechJob
from rdkit import Chem

from autotst.reaction import Reaction, TS
from autotst.species import Species
from autotst.calculators.calculator import Calculator


class Arkane(Calculator):

    def __init__(self, reaction, scratch=".", output_directory=".", model_chemistry="M06-2X/cc-pVTZ", freq_scale_factor=0.982):
        """
        A class to perform Arkane calculations:
        :param: reaction: (Reaction) The reaction of interest
        :param: output_directory: (str) The directory where you would like output files written to
        :param: model_chemistry: (str) The supported model_chemistry described by http://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html#model-chemistry
        :param: freq_scale_factor: (float) The scaling factor corresponding to the model chemistry - source:https://comp.chem.umn.edu/freqscale/version3b1.htm
        """

        self.reaction = reaction
        self.scratch = scratch

        self.arkane_job = RMGArkane()
        self.output_directory = output_directory
        self.arkane_job.outputDirectory = self.output_directory
        self.model_chemistry = model_chemistry
        self.freq_scale_factor = freq_scale_factor

    def get_atoms(self, rmg_mol):
        """
        A method to create an atom dictionary for an rmg molecule
        """
        atom_dict = {}

        for atom in rmg_mol.atoms:
            if atom.isCarbon():
                atom_type = "C"
            if atom.isHydrogen():
                atom_type = "H"
            if atom.isOxygen():
                atom_type = "O"

            try:
                atom_dict[atom_type] += 1
            except KeyError:
                atom_dict[atom_type] = 1

        return atom_dict

    def get_bonds(self, rmg_mol):
        bondList = []
        for atom in rmg_mol.atoms:
            for bond in atom.bonds.values():
                bondList.append(bond)
        bonds = list(set(bondList))
        bondDict = {}
        for bond in bonds:
            if bond.isSingle():
                if bond.atom1.symbol == 'C' and bond.atom2.symbol == 'C':
                    bondType = 'C-C'
                elif (bond.atom1.symbol == 'H' and bond.atom2.symbol == 'H'):
                    bondType = 'H-H'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'H') or (bond.atom1.symbol == 'H' and bond.atom2.symbol == 'C'):
                    bondType = 'C-H'
                elif (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'O'):
                    bondType = 'O-O'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'O') or (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'C'):
                    bondType = 'C-O'
                elif (bond.atom1.symbol == 'H' and bond.atom2.symbol == 'O') or (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'H'):
                    bondType = 'O-H'
                elif bond.atom1.symbol == 'N' and bond.atom2.symbol == 'N':
                    bondType = 'N-N'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'C'):
                    bondType = 'N-C'
                elif (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'O'):
                    bondType = 'N-O'
                elif (bond.atom1.symbol == 'H' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'H'):
                    bondType = 'N-H'
                elif bond.atom1.symbol == 'S' and bond.atom2.symbol == 'S':
                    bondType = 'S-S'
                elif (bond.atom1.symbol == 'H' and bond.atom2.symbol == 'S') or (bond.atom1.symbol == 'S' and bond.atom2.symbol == 'H'):
                    bondType = 'S-H'
            elif bond.isDouble:
                if bond.atom1.symbol == 'C' and bond.atom2.symbol == 'C':
                    bondType = 'C=C'
                elif (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'O'):
                    bondType = 'O=O'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'O') or (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'C'):
                    bondType = 'C=O'
                elif bond.atom1.symbol == 'N' and bond.atom2.symbol == 'N':
                    bondType = 'N=N'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'C'):
                    bondType = 'N=C'
                elif (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'O'):
                    bondType = 'N=O'
                elif (bond.atom1.symbol == 'O' and bond.atom2.symbol == 'S') or (bond.atom1.symbol == 'S' and bond.atom2.symbol == 'O'):
                    bondType = 'S=O'
            elif bond.isTriple:
                if bond.atom1.symbol == 'C' and bond.atom2.symbol == 'C':
                    bondType = 'C#C'
                elif bond.atom1.symbol == 'N' and bond.atom2.symbol == 'N':
                    bondType = 'N#N'
                elif (bond.atom1.symbol == 'C' and bond.atom2.symbol == 'N') or (bond.atom1.symbol == 'N' and bond.atom2.symbol == 'C'):
                    bondType = 'N#C'
            try:
                bondDict[bondType] += 1
            except KeyError:
                bondDict[bondType] = 1

        return bondDict

    def write_arkane_for_reacts_and_prods(self, mol):
        """
        a method to write species to an arkane input file. Mol is an RMGMolecule
        """

        output = ['#!/usr/bin/env python',
                  '# -*- coding: utf-8 -*-', '', 'atoms = {']

        atom_dict = self.get_atoms(mol)

        for atom, count in atom_dict.iteritems():
            output.append("    '{0}': {1},".format(atom, count))
        output = output + ['}', '']

        bond_dict = self.get_bonds(mol)
        if bond_dict != {}:
            output.append('bonds = {')
            for bond_type, num in bond_dict.iteritems():
                output.append("    '{0}': {1},".format(bond_type, num))
            output.append("}")
        else:
            output.append('bonds = {}')

        label = Chem.rdinchi.InchiToInchiKey(
            Chem.MolToInchi(Chem.MolFromSmiles(mol.getSMILES()))).strip("-N")

        external_symmetry = mol.getSymmetryNumber()

        output += ["", "linear = False", "", "externalSymmetry = {}".format(external_symmetry), "",
                   "spinMultiplicity = {}".format(mol.multiplicity), "", "opticalIsomers = 1", ""]

        output += ["energy = {", "    '{0}': Log('{1}.log'),".format(
            self.model_chemistry, label), "}", ""]

        output += ["geometry = Log('{0}.log')".format(label), ""]

        output += [
            "frequencies = Log('{0}.log')".format(label), ""]

        output += ["rotors = []"]

        input_string = ""

        for t in output:
            input_string += t + "\n"


        with open(os.path.join(self.scratch, label +".py"), "w") as f:
            f.write(input_string)

    def write_statmech_ts(self, rxn):
        output = ['#!/usr/bin/env python',
                  '# -*- coding: utf-8 -*-', '', 'atoms = {']

        atom_dict = self.get_atoms(rxn)

        for atom, count in atom_dict.iteritems():
            output.append("    '{0}': {1},".format(atom, count))
        output = output + ['}', '']

        bond_dict = self.get_bonds(rxn)
        if bond_dict != {}:
            output.append('bonds = {')
            for bond_type, num in bond_dict.iteritems():
                output.append("    '{0}': {1},".format(bond_type, num))

            output.append("}")
        else:
            output.append('bonds = {}')

        external_symmetry = ts.getSymmetryNumber()

        output += ["", "linear = False", "", "externalSymmetry = {}".format(external_symmetry), "",
                   "spinMultiplicity = {}".format(rxn.ts.rmg_ts.multiplicity), "", "opticalIsomers = 1", ""]

        output += ["energy = {", "    '{0}': Log('{1}.log'),".format(
            self.model_chemistry, rxn.label), "}", ""]

        output += ["geometry = Log('{0}.log')".format(rxn.label), ""]

        output += [
            "frequencies = Log('{0}.log')".format(rxn.label), ""]

        output += ["rotors = []", ""]

        input_string = ""

        for t in output:
            input_string += t + "\n"

        with open(os.path.join(self.scratch, rxn.label + ".py"), "w") as f:
            f.write(input_string)

    def write_arkane_ts(self, rxn):
        top = ["#!/usr/bin/env python", "# -*- coding: utf-8 -*-", "", 'modelChemistry = "{0}"'.format(
            self.model_chemistry), "frequencyScaleFactor = {0}".format(self.freq_scale_factor), "useHinderedRotors = False", "useBondCorrections = False", ""]

        labels = []

        for react in rxn.reactant_mols:
            label = Chem.rdinchi.InchiToInchiKey(
                Chem.MolToInchi(Chem.MolFromSmiles(react.smiles))).strip("-N")
            if label in labels:
                continue
            else:
                labels.append(label)
            line = "species('{0}', '{1}')".format(
                react.smiles, label + ".py")
            top.append(line)

        for prod in rxn.product_mols:
            label = Chem.rdinchi.InchiToInchiKey(
                Chem.MolToInchi(Chem.MolFromSmiles(prod.smiles))).strip("-N")
            if label in labels:
                continue
            else:
                labels.append(label)
            line = "species('{0}', '{1}')".format(
                prod.smiles, label + ".py")
            top.append(line)

        line = "transitionState('TS', '{0}')".format(rxn.label + ".py")
        top.append(line)


        line = ["",
                "reaction(",
                "    label = '{0}',".format(rxn.label),
                "    reactants = ['{0}', '{1}'],".format(
                    rxn.reactant_mols[0].smiles, rxn.reactant_mols[1].smiles),
                "    products = ['{0}', '{1}'],".format(
                    rxn.product_mols[0].smiles, rxn.product_mols[1].smiles),
                "    transitionState = 'TS',",
                "    tunneling = 'Eckart',",
                ")",
                "",
                "statmech('TS')",
                "kinetics('{0}')".format(rxn.label)]
        top += line

        input_string = ""

        for t in top:
            input_string += t + "\n"

        with open(os.path.join(self.scratch, rxn.label + ".canth.py"), "w") as f:
            f.write(input_string)

    def write_files(self):
        for mol in self.reaction.reactant_mols:

            if len(mol.conformers) < 1:
                logging.info("This molecule object has too")
            rmg_mol = None
            self.write_arkane_for_reacts_and_prods(mol)

        for mol in self.reaction.product_mols:
            self.write_arkane_for_reacts_and_prods(mol)

        self.write_statmech_ts(self.reaction)

        self.write_arkane_ts(self.reaction)

    def run(self):

        self.arkane_job.inputFile = os.path.join(
            self.scratch, self.reaction.label + ".canth.py")
        self.arkane_job.plot = False
        try:
            self.arkane_job.execute()
        except IOError:
            print "There was an issue with Cairo..."

        for job in self.arkane_job.jobList:
            if isinstance(job, KineticsJob):
                self.kinetics_job = job

    def set_reactants_and_products(self):

        for reactant in self.reaction.rmg_reaction.reactants:
            for r in self.kinetics_job.reaction.reactants:
                if reactant.toSMILES() == r.label:
                    r.molecule = [reactant]

        for product in self.reaction.rmg_reaction.products:
            for p in self.kinetics_job.reaction.products:
                if product.toSMILES() == p.label:
                    p.molecule = [product]

        self.reaction.rmg_reaction = self.kinetics_job.reaction

        return self.reaction
