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

import unittest

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
import autotst
import ase
import os
from ase import Atom, Atoms
import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.species import Species as RMGSpecies
from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
from autotst.species import Species, Conformer
import numpy as np

class TestConformer(unittest.TestCase):
    def setUp(self):
        self.conformer = Conformer(smiles='CC')
    def test_rmg_molecules(self):
        self.assertIsInstance(self.conformer.rmg_molecule,RMGMolecule)
    def test_rdkit_mol(self):
        autotst_rdkit = self.conformer.get_rdkit_mol()
        self.assertIsInstance(autotst_rdkit,Mol)
    def test_ase_mol(self):
        autotst_ase_mol = self.conformer.get_ase_mol()
        self.assertIsInstance(autotst_ase_mol,Atoms)
    def test_get_molecules(self):
        autotst_rdkit, autotst_ase_mol = self.conformer.get_molecules()
        self.assertIsInstance(autotst_rdkit,Mol)
        self.assertIsInstance(autotst_ase_mol,Atoms)
    def test_get_bonds(self):
        bonds = self.conformer.get_bonds()
        self.assertIsInstance(bonds,list)
        self.assertIsInstance(bonds[0],Bond)
        self.assertEquals(len(bonds),7)
    def test_get_angles(self):
        angles = self.conformer.get_angles()
        self.assertIsInstance(angles,list)
        self.assertIsInstance(angles[0],Angle)
        self.assertEquals(len(angles),12)
    def test_get_torsions(self):
        torsions = self.conformer.get_torsions()
        self.assertIsInstance(torsions,list)
        self.assertIsInstance(torsions[0],Torsion)
        self.assertEquals(len(torsions),1)
    def test_get_cistrans(self):
        cistrans = self.conformer.get_cistrans()
        self.assertIsInstance(cistrans,list)
        self.assertEquals(len(cistrans),0)
    def test_get_chiralcenters(self):
        chiralcenters = self.conformer.get_chiral_centers()
        self.assertIsInstance(chiralcenters,list)
        self.assertEquals(len(chiralcenters),0)
    def test_get_geometries(self):
        geometries = self.conformer.get_geometries()
        self.assertIsInstance(geometries,tuple)
        self.assertIsInstance(geometries[0],list)
        self.assertIsInstance(geometries[0][0],Bond)
        self.assertIsInstance(geometries[1],list)
        self.assertIsInstance(geometries[1][0],Angle)
        self.assertIsInstance(geometries[2],list)
        self.assertIsInstance(geometries[2][0],Torsion)
        self.assertIsInstance(geometries[3],list)
        self.assertIsInstance(geometries[4],list)
    def test_calculate_symmetry_number(self):
        self.assertEquals(self.conformer.calculate_symmetry_number(),1)
        os.remove("./CC.symm")
    def test_get_xyz_block(self):
        xyz_block = self.conformer.get_xyz_block()
        positions = self.conformer.ase_molecule.arrays["positions"]
        for n in range(len(positions)):
            self.assertTrue((np.array([float(x) for x in xyz_block.split('\n')[n].split()[1:]]) == positions[n]).all())

class TestSpecies(unittest.TestCase):
    def setUp(self):
        self.species = Species(smiles=["CC"])
    def test_generate_structures(self):
        self.assertIsInstance(self.species.generate_structures(),dict)
        self.assertIsInstance(list(self.species.generate_structures().values())[0][0],Conformer)

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
