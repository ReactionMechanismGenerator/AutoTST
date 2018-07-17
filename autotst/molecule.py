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
import logging

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
import autotst
import ase
from ase import Atom, Atoms
import rmgpy
from rmgpy.molecule import Molecule

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

try:
    import py3Dmol
except ImportError:
    logging.info("Error importing py3Dmol")

import numpy as np

from autotst.geometry import CisTrans, Torsion, Angle, Bond


class AutoTST_Molecule():
    """
    A class that allows for one to create RMG, RDKit and ASE
    molecules from a single string with identical atom indicies

    Inputs:
    * smiles (str): a SMILES string that describes the molecule of interest
    * rmg_molecule (RMG Molecule object): an rmg molecule that we can extract information from
    """

    def __init__(self, smiles=None, rmg_molecule=None):

        assert (
            smiles or rmg_molecule), "Please provide a SMILES string and / or an RMG Molecule object."

        if smiles and rmg_molecule:
            assert rmg_molecule.isIsomorphic(
                Molecule(SMILES=smiles)), "SMILES string did not match RMG Molecule object"
            self.smiles = smiles
            self.rmg_molecule = rmg_molecule

        elif rmg_molecule:
            self.rmg_molecule = rmg_molecule
            self.smiles = rmg_molecule.toSMILES()

        else:
            self.smiles = smiles
            self.rmg_molecule = Molecule(SMILES=smiles)

        self.get_rdkit_molecule()
        self.set_rmg_coords("RDKit")
        self.get_ase_molecule()
        self.get_torsions()
        self.get_angles()
        self.get_bonds()

    def __repr__(self):
        return '<AutoTST Molecule "{0}">'.format(self.smiles)

    def get_rdkit_molecule(self):
        """
        A method to create an RDKit Molecule from the rmg_molecule.
        Indicies will be the same as in the RMG Molecule
        """

        RDMol = self.rmg_molecule.toRDKitMol(removeHs=False)

        rdkit.Chem.AllChem.EmbedMolecule(RDMol)

        self.rdkit_molecule = RDMol

    def get_ase_molecule(self):
        """
        A method to create an ASE Molecule from the rdkit_molecule.
        Indicies will be the same as in the RMG and RDKit Molecule.
        """
        mol_list = AllChem.MolToMolBlock(self.rdkit_molecule).split('\n')
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

                        ase_atoms.append(
                            Atom(symbol=symbol, position=(x, y, z)))
                    except:
                        continue

        self.ase_molecule = Atoms(ase_atoms)
        return self.ase_molecule

    def view_mol(self):
        """
        A method designed to create a 3D figure of the AutoTST_Molecule with py3Dmol from the rdkit_molecule
        """
        mb = Chem.MolToMolBlock(self.rdkit_molecule)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mb, "sdf")
        p.setStyle({'stick': {}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

#############################################################################

    def get_bonds(self):

        rdmol_copy = self.rdkit_molecule
        bond_list = []
        for bond in rdmol_copy.GetBonds():
            bond_list.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

        bonds = []
        for indices in bond_list:
            i, j = indices

            length = self.ase_molecule.get_distance(i, j)

            reaction_center = "No"

            bond = Bond(indices=indices, length=length,
                        reaction_center=reaction_center)

            bonds.append(bond)
        self.bonds = bonds
        return self.bonds

    def get_angles(self):

        rdmol_copy = self.rdkit_molecule

        angle_list = []
        for atom1 in rdmol_copy.GetAtoms():
            for atom2 in atom1.GetNeighbors():
                for atom3 in atom2.GetNeighbors():
                    if atom1.GetIdx() == atom3.GetIdx():
                        continue

                    to_add = (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx())
                    if (to_add in angle_list) or (tuple(reversed(to_add)) in angle_list):
                        continue
                    angle_list.append(to_add)

        angles = []
        for indices in angle_list:
            i, j, k = indices

            degree = self.ase_molecule.get_angle(i, j, k)
            ang = Angle(indices=indices, degree=degree,
                        left_mask=[], right_mask=[])
            left_mask = self.get_left_mask(ang)
            right_mask = self.get_right_mask(ang)

            reaction_center = "No"

            angles.append(Angle(indices, degree, left_mask,
                                right_mask, reaction_center))
        self.angles = angles
        return self.angles

    def get_torsions(self):
        rdmol_copy = self.rdkit_molecule

        torsion_list = []
        cistrans_list = []
        for bond1 in rdmol_copy.GetBonds():
            atom1 = bond1.GetBeginAtom()
            atom2 = bond1.GetEndAtom()
            if atom1.IsInRing() or atom2.IsInRing():
                # Making sure that bond1 we're looking at are not in a ring
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
                # if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                # This means that we have a terminal hydrogen, skip this
                # NOTE: for H_abstraction TSs, a non teminal H should exist
                #    continue
                if atomX.GetIdx() != atom2.GetIdx():
                    got_atom0 = True
                    atom0 = atomX

            for bond2 in bond_list2:
                atomY = bond2.GetOtherAtom(atom2)
                # if atomY.GetAtomicNum() == 1 and len(atomY.GetBonds()) == 1:
                # This means that we have a terminal hydrogen, skip this
                #    continue
                if atomY.GetIdx() != atom1.GetIdx():
                    got_atom3 = True
                    atom3 = atomY

            if not (got_atom0 and got_atom3):
                # Making sure atom0 and atom3 were not found
                continue

            # Looking to make sure that all of the atoms are properly bonded to eached
            if ("SINGLE" in str(rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()).GetBondType()) and
                rdmol_copy.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
                rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
                    rdmol_copy.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())):

                torsion_tup = (atom0.GetIdx(), atom1.GetIdx(),
                               atom2.GetIdx(), atom3.GetIdx())

                already_in_list = False
                for torsion_entry in torsion_list:
                    a, b, c, d = torsion_entry
                    e, f, g, h = torsion_tup

                    if (b, c) == (f, g) or (b, c) == (g, f):
                        already_in_list = True

                if not already_in_list:
                    torsion_list.append(torsion_tup)

            if ("DOUBLE" in str(rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()).GetBondType()) and
                rdmol_copy.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
                rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
                    rdmol_copy.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())):

                torsion_tup = (atom0.GetIdx(), atom1.GetIdx(),
                               atom2.GetIdx(), atom3.GetIdx())

                already_in_list = False
                for torsion_entry in torsion_list:
                    a, b, c, d = torsion_entry
                    e, f, g, h = torsion_tup

                    if (b, c) == (f, g) or (b, c) == (g, f):
                        already_in_list = True

                if not already_in_list:
                    cistrans_list.append(torsion_tup)

        torsions = []
        cistrans = []
        for indices in torsion_list:
            i, j, k, l = indices

            dihedral = self.ase_molecule.get_dihedral(i, j, k, l)
            tor = Torsion(indices=indices, dihedral=dihedral,
                          left_mask=[], right_mask=[])
            left_mask = self.get_left_mask(tor)
            right_mask = self.get_right_mask(tor)
            reaction_center = "No"

            torsions.append(Torsion(indices, dihedral,
                                    left_mask, right_mask, reaction_center))

        for indices in cistrans_list:
            i, j, k, l = indices

            dihedral = self.ase_molecule.get_dihedral(i, j, k, l)
            tor = CisTrans(indices=indices, dihedral=dihedral,
                           left_mask=[], right_mask=[])
            left_mask = self.get_left_mask(tor)
            right_mask = self.get_right_mask(tor)
            reaction_center = "No"

            cistrans.append(CisTrans(indices, dihedral,
                                     left_mask, right_mask, reaction_center))
        self.torsions = torsions
        self.cistrans = cistrans
        return self.torsions

    def get_right_mask(self, torsion_or_angle):

        rdmol_copy = self.rdkit_molecule

        rdkit_atoms = rdmol_copy.GetAtoms()

        if (isinstance(torsion_or_angle, autotst.geometry.Torsion) or
                isinstance(torsion_or_angle, autotst.geometry.CisTrans)):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif isinstance(torsion_or_angle, autotst.geometry.Angle):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]

        complete_RHS = False
        i = 0
        atom_index = RHS_atoms_index[0]
        while complete_RHS is False:
            try:
                RHS_atom = rdkit_atoms[atom_index]
                for neighbor in RHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in RHS_atoms_index) or (neighbor.GetIdx() in LHS_atoms_index):
                        continue
                    else:
                        RHS_atoms_index.append(neighbor.GetIdx())
                i += 1
                atom_index = RHS_atoms_index[i]

            except IndexError:
                complete_RHS = True

        right_mask = [index in RHS_atoms_index for index in range(
            len(self.ase_molecule))]

        return right_mask

    def get_left_mask(self, torsion_or_angle):

        rdmol_copy = self.rdkit_molecule

        rdkit_atoms = rdmol_copy.GetAtoms()

        if (isinstance(torsion_or_angle, autotst.geometry.Torsion) or
                isinstance(torsion_or_angle, autotst.geometry.CisTrans)):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif isinstance(torsion_or_angle, autotst.geometry.Angle):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]

        complete_LHS = False
        i = 0
        atom_index = LHS_atoms_index[0]
        while complete_LHS is False:
            try:
                LHS_atom = rdkit_atoms[atom_index]
                for neighbor in LHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in LHS_atoms_index) or (neighbor.GetIdx() in RHS_atoms_index):
                        continue
                    else:
                        LHS_atoms_index.append(neighbor.GetIdx())
                i += 1
                atom_index = LHS_atoms_index[i]

            except IndexError:
                complete_LHS = True

        left_mask = [index in LHS_atoms_index for index in range(
            len(self.ase_molecule))]

        return left_mask

    def set_rmg_coords(self, molecule_base):

        if molecule_base == "RDKit":
            mol_list = AllChem.MolToMolBlock(self.rdkit_molecule).split('\n')
            for i, atom in enumerate(self.rmg_molecule.atoms):
                j = i + 4
                coords = mol_list[j].split()[:3]
                for k, coord in enumerate(coords):
                    coords[k] = float(coord)
                atom.coords = np.array(coords)

        elif molecule_base == "ASE":
            for i, position in enumerate(self.ase_molecule.get_positions()):
                self.rmg_molecule.atoms[i].coords = position

    def update_from_rdkit_mol(self):

        # In order to update the ase molecule you simply need to rerun the get_ase_molecule method
        self.get_ase_molecule()
        self.set_rmg_coords("RDKit")
        # Getting the new torsion angles
        self.get_torsions()

    def update_from_ase_mol(self):

        self.set_rmg_coords("ASE")
        # setting the geometries of the rdkit molecule
        positions = self.ase_molecule.get_positions()
        conf = self.rdkit_molecule.GetConformers()[0]
        for i, atom in enumerate(self.rdkit_molecule.GetAtoms()):
            conf.SetAtomPosition(i, positions[i])

        # Getting the new torsion angles
        self.get_torsions()

    def update_from_rmg_mol(self):

        conf = self.rdkit_molecule.GetConformers()[0]
        ase_atoms = []
        for i, atom in enumerate(self.rmg_molecule.atoms):
            x, y, z = atom.coords
            symbol = atom.symbol

            conf.SetAtomPosition(i, [x, y, z])

            ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

        self.ase_molecule = Atoms(ase_atoms)

        # Getting the new torsion angles
        self.get_torsions()
