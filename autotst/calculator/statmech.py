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
import logging
import rmgpy
from arkane.main import Arkane as RMGArkane, KineticsJob, StatMechJob, ThermoJob
from rdkit import Chem
from ase import Atom, Atoms
from cclib.io import ccread

from autotst.reaction import Reaction, TS
from autotst.species import Species
from autotst.calculators.calculator import Calculator

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)



class StatMech(Calculator):

    def __init__(
            self,
            reaction=None,
            scratch=".",
            output_directory=".",
            model_chemistry="M06-2X/cc-pVTZ",
            freq_scale_factor=0.982):
        """
        A class to perform Arkane calculations:
        :param: reaction: (Reaction) The reaction of interest
        :param: output_directory: (str) The directory where you would like output files written to
        :param: model_chemistry: (str) The supported model_chemistry described by http://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html#model-chemistry
        :param: freq_scale_factor: (float) The scaling factor corresponding to the model chemistry - source:https://comp.chem.umn.edu/freqscale/version3b1.htm
        """
  
        self.reaction = reaction
        self.scratch = scratch

        self.kinetics_job = RMGArkane()
        self.thermo_job = RMGArkane()
        self.output_directory = output_directory
        self.model_chemistry = model_chemistry
        self.freq_scale_factor = freq_scale_factor

    def get_atoms(self, conformer=None):
        """
        A method to create an atom dictionary for an rmg molecule
        """
        atom_dict = {}

        conf = conformer

        rmg_mol = conf.rmg_molecule

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

    def get_bonds(self, conformer=None):

        conf = conformer

        rmg_mol = conf.rmg_molecule

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

    def write_species_files(self, species=None, scratch="."):

        for smiles, confs in species.conformers.items():
            if os.path.exists(os.path.join(scratch, "species", smiles, smiles +".log")):
                logging.info("Lowest energy conformer log file exists for {}".format(smiles))
                write_arkane_conformer(conformer=confs[0], scratch=scratch)
            else:
                logging.info("Lowest energy conformer log file DOES NOT exist for {}".format(smiles))
        
    def write_conformer_file(self, conformer=None, scratch="."):
        
        conf = conformer

        
        label = conf.smiles
        
        if not os.path.exists(os.path.join(scratch,"species", label, label + ".log")):
            logging.info("There is no lowest energy conformer file...")
            return False
        
        parser = ccread(os.path.join(scratch,"species", label, label + ".log"))
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
            atoms.append(Atom(symbol=symbol_dict[atom_num], position=coords))
            
        conf.ase_molecule = Atoms(atoms)
        conf.update_coords_from("ase")
        mol = conf.rmg_molecule

        mol = conf.rmg_molecule
        output = ['#!/usr/bin/env python',
                '# -*- coding: utf-8 -*-', '', 'atoms = {']

        atom_dict = self.get_atoms(conformer=conf) ### Fix this

        for atom, count in atom_dict.iteritems():
            output.append("    '{0}': {1},".format(atom, count))
        output = output + ['}', '']

        bond_dict = self.get_bonds(conformer=conf) ### fix this
        if bond_dict != {}:
            output.append('bonds = {')
            for bond_type, num in bond_dict.iteritems():
                output.append("    '{0}': {1},".format(bond_type, num))
            output.append("}")
        else:
            output.append('bonds = {}')

        external_symmetry = conformer.calculate_symmetry_number()

        output += ["",
                "linear = {}".format(mol.isLinear()),
                "",
                "externalSymmetry = {}".format(external_symmetry),
                "",
                "spinMultiplicity = {}".format(mol.multiplicity),
                "",
                "opticalIsomers = 1",
                ""]
        

        output += ["energy = {", "    '{0}': Log('{1}.log'),".format(
            self.model_chemistry, label), "}", ""] ###fix this

        output += ["geometry = Log('{0}.log')".format(label), ""]

        output += [
            "frequencies = Log('{0}.log')".format(label), ""]

        """
        TODO: add rotor information @carl
        output += ["rotors = ["]
        for torsion in conf.torsions:
            output += [self.get_rotor_info(conf, torsion)]
        output += ["]"]
        """

        input_string = ""

        for t in output:
            input_string += t + "\n"


        with open(os.path.join(scratch,"species", label, label + '.py'), "w") as f:
            f.write(input_string)
        return True

    def get_rotor_info(self, conformer, torsion):
        """
        Formats and returns info about torsion as it should appear in an Arkane species.py

        conformer :: autotst conformer object
        torsion   :: autotst torsion object 

        Needed for Arkane species file:
        scanLog :: Gaussian output log of freq calculation on optimized geometry
        pivots  :: torsion center: j,k of i,j,k,l (Note Arkane begins indexing with 1)
        top     :: ID of all atoms in one top (Note Arkane begins indexing with 1)

        """
        i,j,k,l = torsion.atom_indices

        tor_center = [j,k]
        tor_center_adj = [j+1, k+1] # Adjusted since mol's IDs start from 0 while Arkane's start from 1

        if isinstance(conformer, TS):
            tor_log = os.path.join(
                self.scratch,
                comformer.reaction_label + "_tor{0}{1}.log".format(j,k)
            )
        else:
            tor_log = os.path.join(
                self.scratch,
                conformer.rmg_molecule.toAugmentedInChIKey().strip("-N") + '_tor{0}{1}.log'.format(j,k)
            )

        top_IDs = []
        for num, tf in enumerate(torsion.mask):
            if tf: top_IDs.append(num)

        top_IDs_adj = [ID+1 for ID in top_IDs] # Adjusted to start from 1 instead of 0

        info = "     HinderedRotor(scanLog=Log('{0}'), pivots={1}, top={2}, fit='fourier'),".format(tor_log, tor_center_adj, top_IDs_adj)

        return info


    def write_ts_input(self, transitionstate=None, scratch=None):
        
        label = transitionstate.reaction_label
        
        if not os.path.exists(os.path.join(scratch,"ts", label, label + ".log")):
            logging.info("There is no lowest energy conformer file...")
            return False
        
        parser = ccread(os.path.join(scratch,"ts", label, label + ".log"))
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
            atoms.append(Atom(symbol=symbol_dict[atom_num], position=coords))
            
        transitionstate.ase_molecule = Atoms(atoms)
        transitionstate.update_coords_from("ase")
        mol = transitionstate.rmg_molecule

        output = ['#!/usr/bin/env python',
                '# -*- coding: utf-8 -*-', '', 'atoms = {']

        atom_dict = self.get_atoms(conformer=transitionstate) # need to fix

        for atom, count in atom_dict.iteritems():
            output.append("    '{0}': {1},".format(atom, count))
        output = output + ['}', '']

        bond_dict = self.get_bonds(conformer=transitionstate) # need to fix
        if bond_dict != {}:
            output.append('bonds = {')
            for bond_type, num in bond_dict.iteritems():
                output.append("    '{0}': {1},".format(bond_type, num))

            output.append("}")
        else:
            output.append('bonds = {}')
        transitionstate.rmg_molecule.updateMultiplicity()
        
        external_symmetry = transitionstate.calculate_symmetry_number()

        output += ["",
                "linear = False",
                "",
                "externalSymmetry = {}".format(external_symmetry),
                "",
                "spinMultiplicity = {}".format(transitionstate.rmg_molecule.multiplicity),
                "",
                "opticalIsomers = 1",
                ""]

        output += ["energy = {", "    '{0}': Log('{1}.log'),".format(
            self.model_chemistry, label), "}", ""] #fix this 

        output += ["geometry = Log('{0}.log')".format(label), ""]

        output += [
            "frequencies = Log('{0}.log')".format(label), ""]

        output += ["rotors = []", ""] ### TODO: Fix this

        input_string = ""

        for t in output:
            input_string += t + "\n"

        with open(os.path.join(scratch,"ts", label, label + '.py'), "w") as f:
            f.write(input_string)
        return True

    def write_kinetics_input(self, reaction=None, scratch="."):
        
        top = [
            "#!/usr/bin/env python",
            "# -*- coding: utf-8 -*-",
            "",
            'modelChemistry = "{0}"'.format(
                self.model_chemistry), #fix this
            "frequencyScaleFactor = {0}".format(
                self.freq_scale_factor), #fix this
            "useHinderedRotors = False", #fix this @carl
            "useBondCorrections = False",
            ""]

        labels = []
        r_smiles = []
        p_smiles = []
        for i, react in enumerate(reaction.reactants):
            lowest_energy = 1e5
            lowest_energy_conf = None
            
            if len(react.conformers.keys()) > 1:
                for smiles in react.conformers.keys():
                    path = os.path.join(scratch, "species", smiles, smiles + ".log")
                    if not os.path.exists(path):
                        logging.info("It looks like {} doesn't have any optimized geometries".format(smiles))
                        continue
                    
                    parser = ccread(path)
                    energy = parser.scfenergies[-1]
                    if energy < lowest_energy:
                        lowest_energy = energy
                        lowest_energy_conf = react.conformers[smiles][0]
                            
            else:
                smiles = react.conformers.keys()[0]
                path = os.path.join(scratch, "species", smiles, smiles + ".log")
                if not os.path.exists(path):
                    logging.info("It looks like {} doesn't have any optimized geometries".format(smiles))
                    continue

                parser = ccread(path)
                lowest_energy = parser.scfenergies[-1]
                lowest_energy_conf = react.conformers.values()[0][0]

            #r_smiles.append(lowest_energy_conf.smiles)
            r_smiles.append("react_{}".format(i))
            label = lowest_energy_conf.smiles
            if label in labels:
                continue
            else:
                labels.append(label)
            line = "species('{0}', '{1}', structure=SMILES('{2}'))".format(
                "react_{}".format(i), os.path.join(scratch,"species", label, label + ".py"), label)
            top.append(line)
                
        for i, prod in enumerate(reaction.products):

            lowest_energy = 1e5
            lowest_energy_conf = None
            
            
            if len(prod.conformers.keys()) > 1:
                for smiles in prod.conformers.keys():
                    path = os.path.join(scratch, "species", smiles, smiles + ".log")
                    if not os.path.exists(path):
                        logging.info("It looks like {} doesn't have any optimized geometries".format(smiles))
                        continue
                    
                    parser = ccread(path)
                    energy = parser.scfenergies[-1]
                    if energy < lowest_energy:
                        lowest_energy = energy
                        lowest_energy_conf = prod.conformers[smiles][0]
                            
            else:
                smiles = prod.conformers.keys()[0]
                path = os.path.join(scratch, "species", smiles, smiles + ".log")
                if not os.path.exists(path):
                    logging.info("It looks like {} doesn't have any optimized geometries".format(smiles))
                    continue

                parser = ccread(path)
                lowest_energy = parser.scfenergies[-1]
                lowest_energy_conf = prod.conformers.values()[0][0]

            #p_smiles.append(lowest_energy_conf.smiles)
            p_smiles.append("prod_{}".format(i))
            label = lowest_energy_conf.smiles
            if label in labels:
                continue
            else:
                labels.append(label)
            line = "species('{0}', '{1}', structure=SMILES('{2}'))".format(
                "prod_{}".format(i), os.path.join(scratch, "species", label, label + ".py"), label)
            top.append(line)

        line = "transitionState('TS', '{0}')".format(os.path.join(scratch, "ts", reaction.label, reaction.label + ".py"))
        top.append(line)

        line = ["",
                "reaction(",
                "    label = '{0}',".format(reaction.label),
                "    reactants = {},".format(
                    r_smiles),
                "    products = {},".format(
                    p_smiles),
                "    transitionState = 'TS',",
                "    tunneling = 'Eckart',",
                ")",
                "",
                "statmech('TS')",
                "kinetics('{0}')".format(reaction.label)]

        top += line
        
        input_string = ""

        for t in top:
            input_string += t + "\n"

        with open(os.path.join(scratch, "ts", reaction.label, reaction.label + ".kinetics.py"), "w") as f:
            f.write(input_string)

    def write_thermo_input(self, conformer=None, scratch="."): 
        model_chemistry="M06-2X/cc-pVTZ"
        freq_scale_factor=0.982
        
        top = [
            "#!/usr/bin/env python",
            "# -*- coding: utf-8 -*-",
            "",
            'modelChemistry = "{0}"'.format(
                model_chemistry), #fix this
            "frequencyScaleFactor = {0}".format(
                freq_scale_factor), #fix this
            "useHinderedRotors = False", #fix this @carl
            "useBondCorrections = False",
            ""]


        line = "species('species', '{1}', structure=SMILES('{2}'))".format(
            "species", os.path.join(conformer.smiles + ".py"), conformer.smiles)
        top.append(line)
            
        top.append("statmech('species')")
        top.append("thermo('species', 'NASA')")

        input_string = ""

        for t in top:
            input_string += t + "\n"

        with open(os.path.join(scratch, "species", conformer.smiles, conformer.smiles + ".thermo.py"), "w") as f:
            f.write(input_string)
########################################

    def write_files(self):

        for mol in self.reaction.reactants:
            for smiles, confs in mol.conformers.items():
                conf =  confs[0]
                self.write_conformer_file(conf, scratch=self.scratch)

        for mol in self.reaction.products:
            for smiles, confs in mol.conformers.items():
                conf =  confs[0]
                self.write_conformer_file(conf, scratch=self.scratch)

        self.write_ts_input(self.reaction.ts["forward"][0], scratch=self.scratch)

        self.write_kinetics_input(self.reaction, scratch=self.scratch)

    def run(self):
        
        self.kinetics_job.inputFile = os.path.join(
            self.scratch, "ts", self.reaction.label, self.reaction.label + ".kinetics.py")
        self.kinetics_job.plot = False
        self.kinetics_job.outputDirectory = os.path.join(self.scratch, "ts")
        # try:
        self.kinetics_job.execute()
        # except IOError:

        for job in self.kinetics_job.jobList:
            if isinstance(job, KineticsJob):
                self.kinetics_job = job
            elif isinstance(job, ThermoJob):
                self.thermo_job = job

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
