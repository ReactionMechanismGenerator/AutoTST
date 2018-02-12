import os
import sys
import cPickle as pkl
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import re
import imp
import itertools
import random
import numpy as np
from numpy import array
import pandas as pd

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdChemReactions import ChemicalReaction

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel

from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDepository, KineticsRules
from rmgpy.qm.main import QMCalculator, QMSettings
from rmgpy.qm.qmdata import QMData
from rmgpy.qm.reaction import QMReaction
from rmgpy.qm.molecule import QMMolecule

from multi_molecule import *
from geometry import *


rmg_database = RMGDatabase()
database_path = os.path.abspath(os.path.join(os.getenv('RMGpy', '..'), '..', 'RMG-database', 'input'))
rmg_database.load(database_path,
                 kineticsFamilies=['H_Abstraction'],
                 transportLibraries=[],
                 reactionLibraries=[],
                 seedMechanisms=[],
                 thermoLibraries=['primaryThermoLibrary', 'KlippensteinH2O2', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR' ],
                 solvation=False,
                 )
home = os.path.expandvars("$HOME")
ts_file = home + "/Code/ga_conformer/python_code/ts_database.pkl"
f = open(ts_file, "r")
ts_database = pkl.load(f)

settings = QMSettings(
    software='gaussian',
    method='m062x',
    fileStore=os.path.expandvars('.'),
    scratchDirectory=os.path.expandvars('.'),
    )

class Multi_Reaction():

    def __init__(self, label, reaction_family, rmg_reaction = None):
        self.label = label
        self.reaction_family = reaction_family

        if rmg_reaction:

            reactant_mols = []
            product_mols = []

            reactants = rmg_reaction.reactants
            products = rmg_reaction.products

            for reactant in reactants:
                reactant_mols.append(Multi_Molecule(reactant.SMILES))

            for product in products:
                product_mols.append(Multi_Molecule(product))

            self.reactant_mols = reactant_mols
            self.product_mols = product_mols
        else:
            self.get_reactants_and_products()

        self.get_rmg_reactions()
        self.create_ts_geometries()

    def get_reactants_and_products(self):

        """
        This uses the reaction label to creat multi_molecule objects of

        :return:
        """
        reactants, products = self.label.split("_")

        if "+" in reactants:
            reactants = reactants.split("+")

        if "+" in products:
            products = products.split("+")

        reactant_mols = []
        product_mols = []

        for reactant in reactants:
            reactant_mols.append(Multi_Molecule(reactant))

        for product in products:
            product_mols.append(Multi_Molecule(product))

        self.reactant_mols = reactant_mols
        self.product_mols = product_mols


    def get_rmg_reactions(self):
        """
        This method creates a labeled rmg_reaction from the reaction string

        the reaction string should look as follows: r1+r2_p1+p2
        """

        rmg_reactants = []
        rmg_products = []

        for reactant_mol in self.reactant_mols:
            rmg_reactants.append(reactant_mol.rmg_molecule)

        for product_mol in self.product_mols:
            rmg_products.append(product_mol.rmg_molecule)

        test_reaction = Reaction(reactants=rmg_reactants, products=rmg_products, reversible=True)

        reaction_list = rmg_database.kinetics.generateReactionsFromFamilies(
            rmg_reactants,
            rmg_products,
            only_families=self.reaction_family)

        for reaction in reaction_list:
            # Check if any of the RMG proposed reactions matches the reaction in the mechanism
            if test_reaction.isIsomorphic(reaction):
                atom_labels_reactants = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
                atom_labels_products = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

                for reactant in reaction.reactants:
                    reactant.clearLabeledAtoms()
                    for atom in reactant.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_reactants[atom_label[0]] = True

                for product in reaction.products:
                    product.clearLabeledAtoms()
                    for atom in product.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_products[atom_label[0]] = True

                if all(atom_labels_reactants.values()) and all(atom_labels_products.values()):
                    # We successfully labeled all of the atoms
                    break
        self.rmg_reaction = reaction
        self.rmg_qm_reaction = QMReaction(reaction=reaction, settings=settings, tsDatabase=ts_database)

    def create_ts_geometries(self):
        """
        A method to use the tools in rmg / autotst to create a reasonable TS geometry
        This will create the geometry in both rdkit and ase

        :return:
        self.multi_ts: a multi_ts object that contains geometries of a ts in
                        rdkit, ase, and rmg molecules
        """
        self.multi_ts = Multi_TS(self)


class Multi_TS():
    def __init__(self, Multi_Reaction):

        self.multi_reaction = Multi_Reaction
        self.label = Multi_Reaction.label  # make sure that both the reaction and TS have same label

        self.create_rdkit_ts_geometry()
        self.create_ase_ts_geometry()
        self.create_rmg_ts_geometry()
        self.get_ts_torsion_list()
        self.get_ts_torsions()
        self.get_ts_angle_indices()
        self.angle_indices
        self.get_ts_angle()
        i, j, k =self.angle_indices
        r_mask = self.angle.right_mask

        self.ase_ts.set_angle(a1=i, a2=j, a3=k, angle=float(180), mask=r_mask)
        labels = []
        for label, atom in self.rmg_ts.getLabeledAtoms().iteritems():
            labels.append(atom.sortingLabel)

        for torsion in self.torsions:
            if set(labels).issubset(torsion.indices[:3]):
                i, j, k, l = torsion.indices
                r_mask = torsion.right_mask

                self.ase_ts.set_dihedral(a1=i,
                                         a2=j,
                                         a3=k,
                                         a4=l,
                                         mask=r_mask,
                                         angle=float(180))
        self.update_ts_from_ase_ts()

    def create_rdkit_ts_geometry(self):

        self.rmg_ts, product = self.multi_reaction.rmg_qm_reaction.setupMolecules()

        labels, atom_match = self.multi_reaction.rmg_qm_reaction.getLabels(self.rmg_ts)

        self.rdkit_ts, bm, self.multi_reaction.rmg_qm_reaction.reactantGeom = self.multi_reaction.rmg_qm_reaction.generateBoundsMatrix(self.rmg_ts)

        bm = self.multi_reaction.rmg_qm_reaction.editMatrix(self.rmg_ts, bm, labels)

        self.rdkit_ts = self.multi_reaction.rmg_qm_reaction.reactantGeom.rd_embed(self.rdkit_ts, 1000, bm=bm, match=atom_match)[0]

    def create_ase_ts_geometry(self):

        mol_list = AllChem.MolToMolBlock(self.rdkit_ts).split('\n')
        ase_atoms = []
        for i, line in enumerate(mol_list):

            if i > 3:

                try:
                    atom0, atom1, bond, rest = line
                    atom0 = int(atom0)
                    atom0 = int(atom1)
                    bond = float(bond)

                except ValueError:
                    try:
                        x, y, z, symbol = line.split()[0:4]
                        x = float(x)
                        y = float(y)
                        z = float(z)
                        # print symbol

                        ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

                    except:
                        continue

        self.ase_ts = Atoms(ase_atoms)

    def create_rmg_ts_geometry(self):

        for i, position in enumerate(self.ase_ts.get_positions()):
            self.rmg_ts.atoms[i].coords = position

    def view_ts(self):
        """
        A method designed to create a 3D figure of the Multi_Molecule with py3Dmol
        """
        mb = Chem.MolToMolBlock(self.rdkit_ts)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mb, "sdf")
        p.setStyle({'stick':{}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

    def create_pseudo_geometry(self):

        rdmol_copy = self.rdkit_ts.__copy__()
        rdmol_copy = Chem.RWMol(rdmol_copy)
        for atom in rdmol_copy.GetAtoms():
            idx = atom.GetIdx()
            rmg_atom = self.rmg_ts.atoms[idx]

            if rmg_atom.label:
                if rmg_atom.label == "*1":
                    atom1_star = atom
                if rmg_atom.label == "*2":
                    atom2_star = atom
                if rmg_atom.label == "*3":
                    atom3_star = atom

        try:
            rdmol_copy.AddBond(atom1_star.GetIdx(), atom2_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)
        except RuntimeError:
            rdmol_copy.AddBond(atom2_star.GetIdx(), atom3_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)

        return rdmol_copy

    def get_ts_angle_indices(self):
        rdmol_copy = self.rdkit_ts.__copy__()
        rdmol_copy = Chem.RWMol(rdmol_copy)
        for atom in rdmol_copy.GetAtoms():
            idx = atom.GetIdx()
            rmg_atom = self.rmg_ts.atoms[idx]

            if rmg_atom.label:
                if rmg_atom.label == "*1":
                    atom1_star = atom
                if rmg_atom.label == "*2":
                    atom2_star = atom
                if rmg_atom.label == "*3":
                    atom3_star = atom

        if rdmol_copy.GetBondBetweenAtoms(atom1_star.GetIdx(), atom2_star.GetIdx()):
            bond_list = list(atom3_star.GetBonds())

            for bond in bond_list:
                atomX = bond.GetOtherAtom(atom3_star)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2_star.GetIdx():
                    other_atom = atomX
            angle_indices = [atom2_star.GetIdx(), atom3_star.GetIdx(), other_atom.GetIdx()]

        else:
            bond_list = list(atom1_star.GetBonds())

            for bond in bond_list:
                atomX = bond.GetOtherAtom(atom1_star)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2_star.GetIdx():
                    other_atom = atomX

            angle_indices = [atom2_star.GetIdx(), atom1_star.GetIdx(), other_atom.GetIdx()]

        self.angle_indices = angle_indices

    def get_ts_angle(self):

        i, j, k = self.angle_indices
        ang = self.ase_ts.get_angle(i, j, k)
        angle = Angle(indices=self.angle_indices, degree=ang, left_mask=[], right_mask=[])
        left_mask = self.get_ts_left_mask(angle)
        right_mask = self.get_ts_right_mask(angle)

        self.angle = Angle(self.angle_indices, ang, left_mask, right_mask)

        return self.angle

    def get_ts_torsion_list(self):

        rdmol_copy = self.create_pseudo_geometry()

        torsion_list = []
        for bond1 in rdmol_copy.GetBonds():
            atom1 = bond1.GetBeginAtom()
            atom2 = bond1.GetEndAtom()
            if atom1.IsInRing() or atom2.IsInRing():
                # Making sure that bond1 we're looking at are in a ring
                continue

            bond_list1 = list(atom1.GetBonds())
            bond_list2 = list(atom2.GetBonds())

            if not len(bond_list1) > 1 and not len(bond_list2) > 1:
                # Making sure that there are more than one bond attached to
                # the atoms we're looking at
                continue

            # Getting the 0th and 3rd atom and insuring that atoms
            # attached to the 1st and 2nd atom are not terminal hydrogens
            # We also make sure that all of the atoms are properly bound together

            # If the above are satisfied, we append a tuple of the torsion our torsion_list
            got_atom0 = False
            got_atom3 = False

            for bond0 in bond_list1:
                atomX = bond0.GetOtherAtom(atom1)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2.GetIdx():
                    got_atom0 = True
                    atom0 = atomX

            for bond2 in bond_list2:
                atomY = bond2.GetOtherAtom(atom2)
                if atomY.GetAtomicNum() == 1 and len(atomY.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    continue
                if atomY.GetIdx() != atom1.GetIdx():
                    got_atom3 = True
                    atom3 = atomY

            if not (got_atom0 and got_atom3):
                # Making sure atom0 and atom3 were not found
                continue

            # Looking to make sure that all of the atoms are properly bonded to eached
            if (
                            rdmol_copy.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
                            rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
                        rdmol_copy.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())):
                torsion_tup = (atom0.GetIdx(), atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx())
                torsion_list.append(torsion_tup)

        self.torsion_list = torsion_list

    def get_ts_torsions(self):
        torsions = []
        for indices in self.torsion_list:
            i, j, k, l = indices

            dihedral = self.ase_ts.get_dihedral(i, j, k, l)
            tor = Torsion(indices=indices, dihedral=dihedral, left_mask=[], right_mask=[])
            left_mask = self.get_ts_left_mask(tor)
            right_mask = self.get_ts_right_mask(tor)

            torsions.append(Torsion(indices, dihedral, left_mask, right_mask))
        self.torsions = torsions
        return self.torsions

    def get_ts_right_mask(self, torsion_or_angle):

        rdmol_copy = self.create_pseudo_geometry()

        rdkit_atoms = rdmol_copy.GetAtoms()

        if "Torsion" in str(torsion_or_angle.__class__):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif "Angle" in str(torsion_or_angle.__class__):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]

        complete_RHS = False
        i = 0
        atom_index = RHS_atoms_index[0]
        while complete_RHS == False:
            try:
                RHS_atom = rdkit_atoms[atom_index]
                for neighbor in RHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in RHS_atoms_index) or (neighbor.GetIdx() in LHS_atoms_index):
                        continue
                    else:
                        RHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = RHS_atoms_index[i]

            except IndexError:
                complete_RHS = True

        right_mask = [index in RHS_atoms_index for index in range(len(self.ase_ts))]

        return right_mask

    def get_ts_left_mask(self, torsion_or_angle):

        rdmol_copy = self.create_pseudo_geometry()

        rdkit_atoms = rdmol_copy.GetAtoms()

        if "Torsion" in str(torsion_or_angle.__class__):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif "Angle" in str(torsion_or_angle.__class__):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]


        complete_LHS = False
        i = 0
        atom_index = LHS_atoms_index[0]
        while complete_LHS == False:
            try:
                LHS_atom = rdkit_atoms[atom_index]
                for neighbor in LHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in LHS_atoms_index) or (neighbor.GetIdx() in RHS_atoms_index):
                        continue
                    else:
                        LHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = LHS_atoms_index[i]

            except IndexError:
                complete_LHS = True

        left_mask = [index in LHS_atoms_index for index in range(len(self.ase_ts))]

        return left_mask

    def set_rmg_ts_coords(self, molecule_base):

        if molecule_base == "RDKit":
            mol_list = AllChem.MolToMolBlock(self.rdkit_ts).split('\n')
            for i, atom in enumerate(self.rmg_ts.atoms):

                j = i + 4
                coords = mol_list[j].split()[:3]

                for k, coord in enumerate(coords):

                    coords[k] = float(coord)
                atom.coords = np.array(coords)

        elif molecule_base == "ASE":
            for i, position in enumerate(self.ase_ts.get_positions()):
                self.rmg_ts.atoms[i].coords = position

    def update_ts_from_rdkit_ts(self):
        # In order to update the ase molecule you simply need to rerun the get_ase_molecule method
        self.create_ase_ts_geometry()
        self.set_rmg_ts_coords("RDKit")

        # Getting the new torsion angles
        self.get_ts_torsions()

    def update_ts_from_ase_ts(self):

        self.set_rmg_ts_coords("ASE")

        # setting the geometries of the rdkit molecule

        positions = self.ase_ts.get_positions()

        conf = self.rdkit_ts.GetConformers()[0]

        for i, atom in enumerate(self.rdkit_ts.GetAtoms()):
            conf.SetAtomPosition(i, positions[i])

        # Getting the new torsion angles
        self.get_ts_torsions()


    def update_ts_from_rmg_ts(self):

        conf = self.rdkit_ts.GetConformers()[0]
        ase_atoms = []
        for i, atom in enumerate(self.rmg_ts.atoms):
            x, y, z = atom.coords
            symbol = atom.symbol

            conf.SetAtomPosition(i, [x, y, z])

            ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

        self.ase_ts = Atoms(ase_atoms)

        # Getting the new torsion angles
        self.get_ts_torsion_list()
