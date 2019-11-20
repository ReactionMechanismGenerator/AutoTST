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

import os, sys
import unittest
from rdkit.Chem.rdchem import Mol, RWMol
from ase import Atoms
from autotst.reaction import Reaction, TS
from autotst.data.base import TransitionStates
from rmgpy.reaction import Reaction as RMGReaction
from rmgpy.species import Species as RMGSpecies
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.data.rmg import RMGDatabase



class TestReaction(unittest.TestCase):
    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.reaction2 = Reaction(rmg_reaction=RMGReaction(
            reactants = [RMGMolecule(smiles="CC"), RMGMolecule(smiles="[O]O")],
            products = [RMGMolecule(smiles="[CH2]C"), RMGMolecule(smiles="OO")]
        ))

    def test_label(self):
        self.assertEqual(self.reaction.label, "CC+[O]O_[CH2]C+OO")
        self.assertEqual(self.reaction2.get_label(), "CC+[O]O_C[CH2]+OO")

    def test_rmg_reaction(self):
        test_reaction = RMGReaction(
            reactants=[
                RMGMolecule(smiles="CC"),
                RMGMolecule(smiles="[O]O")
            ],
            products = [
                RMGMolecule(smiles="[CH2]C"),
                RMGMolecule(smiles="OO")
            ]
        )

        self.assertTrue(test_reaction.is_isomorphic(self.reaction.get_rmg_reaction()))
        self.assertTrue(test_reaction.is_isomorphic(self.reaction2.get_rmg_reaction()))


    def test_databases(self):
        rmg_database, ts_databases = self.reaction.load_databases()

        self.assertIsInstance(rmg_database, RMGDatabase)
        self.assertIsInstance(ts_databases, dict)
        self.assertIsInstance(ts_databases["H_Abstraction"], TransitionStates)

        rmg_database, ts_databases = self.reaction2.load_databases()

        self.assertIsInstance(rmg_database, RMGDatabase)
        self.assertIsInstance(ts_databases, dict)
        self.assertIsInstance(ts_databases["H_Abstraction"], TransitionStates)

    def test_disance_data(self):
        d12 = 1.38
        d13 = 2.43 
        # d13 in reactionTest is smaller than the distance in baseTest.py because
        # d13 is edited to be smaller in reaction.py. The value returned from the database
        # is ~2.53 but is reduced to 2.43 when called from the reaction object itself
        d23 = 1.16

        self.assertAlmostEquals(d12, self.reaction.distance_data.distances["d12"], places=1)
        self.assertAlmostEquals(d13, self.reaction.distance_data.distances["d13"], places=1)
        self.assertAlmostEquals(d23, self.reaction.distance_data.distances["d23"], places=1)

        self.assertAlmostEquals(d12, self.reaction2.distance_data.distances["d12"], places=1)
        self.assertAlmostEquals(d13, self.reaction2.distance_data.distances["d13"], places=1)
        self.assertAlmostEquals(d23, self.reaction2.distance_data.distances["d23"], places=1)

    def test_generate_reactants_and_products(self):
        reactants, products = self.reaction.generate_reactants_and_products()

        self.assertIsInstance(reactants, list)
        self.assertIsInstance(products, list)
        self.assertEquals(len(reactants), 2)
        self.assertAlmostEquals(len(products), 2)

        reactants, products = self.reaction2.generate_reactants_and_products()

        self.assertIsInstance(reactants, list)
        self.assertIsInstance(products, list)
        self.assertEquals(len(reactants), 2)
        self.assertAlmostEquals(len(products), 2)

    def test_labeled_reaction(self):
        test_reaction = RMGReaction(
            reactants=[
                RMGMolecule(smiles="CC"),
                RMGMolecule(smiles="[O]O")
            ],
            products = [
                RMGMolecule(smiles="[CH2]C"),
                RMGMolecule(smiles="OO")
            ]
        )
        labeled_reaction, reaction_family = self.reaction.get_labeled_reaction()

        self.assertEquals(reaction_family.lower(), "h_abstraction")
        self.assertTrue(test_reaction.is_isomorphic(labeled_reaction))

        merged = labeled_reaction.reactants[0].merge(labeled_reaction.reactants[1])
        self.assertTrue(merged.get_labeled_atoms("*1")[0].is_carbon())
        self.assertTrue(merged.get_labeled_atoms("*2")[0].is_hydrogen())
        self.assertTrue(merged.get_labeled_atoms("*3")[0].is_oxygen)

        merged = labeled_reaction.products[0].merge(labeled_reaction.products[1])
        self.assertTrue(merged.get_labeled_atoms("*3")[0].is_carbon())
        self.assertTrue(merged.get_labeled_atoms("*2")[0].is_hydrogen())
        self.assertTrue(merged.get_labeled_atoms("*1")[0].is_oxygen)


        labeled_reaction, reaction_family = self.reaction2.get_labeled_reaction()

        self.assertEquals(reaction_family.lower(), "h_abstraction")
        self.assertTrue(test_reaction.is_isomorphic(labeled_reaction))

        merged = labeled_reaction.reactants[0].merge(labeled_reaction.reactants[1])
        self.assertTrue(merged.get_labeled_atoms("*1")[0].is_carbon())
        self.assertTrue(merged.get_labeled_atoms("*2")[0].is_hydrogen())
        self.assertTrue(merged.get_labeled_atoms("*3")[0].is_oxygen)

        merged = labeled_reaction.products[0].merge(labeled_reaction.products[1])
        self.assertTrue(merged.get_labeled_atoms("*3")[0].is_carbon())
        self.assertTrue(merged.get_labeled_atoms("*2")[0].is_hydrogen())
        self.assertTrue(merged.get_labeled_atoms("*1")[0].is_oxygen)

    def test_rmg_complexes(self):
        self.reaction.get_labeled_reaction()
        self.reaction.get_rmg_complexes()

        self.assertEquals(len(self.reaction.complexes), 2)
        self.assertEquals(len(self.reaction.complexes["forward"].get_all_labeled_atoms()), 3)
        self.assertEquals(len(self.reaction.complexes["reverse"].get_all_labeled_atoms()), 3)


        self.reaction2.get_labeled_reaction()
        self.reaction2.get_rmg_complexes()

        self.assertEquals(len(self.reaction2.complexes), 2)
        self.assertEquals(len(self.reaction2.complexes["forward"].get_all_labeled_atoms()), 3)
        self.assertEquals(len(self.reaction2.complexes["reverse"].get_all_labeled_atoms()), 3)

    def test_ts(self):
        self.assertEquals(len(self.reaction.ts), 2)
        self.assertEquals(len(self.reaction.ts["forward"]), 1)
        self.assertEquals(len(self.reaction.ts["reverse"]), 1)
        self.assertIsInstance(self.reaction.ts["forward"][0], TS)
        self.assertIsInstance(self.reaction.ts["reverse"][0], TS)

        self.assertEquals(len(self.reaction2.ts), 2)
        self.assertEquals(len(self.reaction2.ts["forward"]), 1)
        self.assertEquals(len(self.reaction2.ts["reverse"]), 1)
        self.assertIsInstance(self.reaction2.ts["forward"][0], TS)
        self.assertIsInstance(self.reaction2.ts["reverse"][0], TS)

    def test_reaction_families(self):
        # R_Addition_MultipleBond
        reaction = Reaction("C#C+[OH]_[CH]=CO")
        _, family = reaction.get_labeled_reaction()
        self.assertEqual(family.lower(), "R_Addition_MultipleBond".lower())
        
        # intra_H_migration
        reaction = Reaction("C[CH]O_CC[O]")
        _, family = reaction.get_labeled_reaction()
        self.assertEqual(family.lower(), "intra_H_migration".lower())

class TestTS(unittest.TestCase):

    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.reaction2 = Reaction(rmg_reaction=RMGReaction(
            reactants = [RMGMolecule(smiles="CC"), RMGMolecule(smiles="[O]O")],
            products = [RMGMolecule(smiles="[CH2]C"), RMGMolecule(smiles="OO")]
        ))
        self.ts = self.reaction.ts["forward"][0]
        self.ts2 = self.reaction.ts["forward"][0]
        self.ts.get_molecules()
        self.ts2.get_molecules()
        
    def test_reaction_label(self):
        self.assertEquals(self.ts.reaction_label, "CC+[O]O_[CH2]C+OO")
        self.assertEquals(self.ts2.reaction_label, "CC+[O]O_[CH2]C+OO")
    
    def test_smiles(self):
        self.assertEquals(self.ts.smiles, "CC.[O]O")
        self.assertEquals(self.ts2.smiles, "CC.[O]O")
    
    def test_rdkit_molecule(self):
        rdkit_molecule = self.ts.rdkit_molecule
        self.assertIsInstance(rdkit_molecule, Mol)
        self.assertEquals(len(rdkit_molecule.GetAtoms()), 11)
        self.assertEquals(len(rdkit_molecule.GetBonds()), 9)

        rdkit_molecule = self.ts2.rdkit_molecule
        self.assertIsInstance(rdkit_molecule, Mol)
        self.assertEquals(len(rdkit_molecule.GetAtoms()), 11)
        self.assertEquals(len(rdkit_molecule.GetBonds()), 9)
    
    def test_pseudo_geometry(self):
        self.assertIsInstance(self.ts._pseudo_geometry, RWMol)
        self.assertEquals(len(self.ts._pseudo_geometry.GetBonds()), 10)

        self.assertIsInstance(self.ts2._pseudo_geometry, RWMol)
        self.assertEquals(len(self.ts2._pseudo_geometry.GetBonds()), 10)

    def test_ase_molecule(self):
        ase_molecule = self.ts.ase_molecule
        self.assertIsInstance(ase_molecule, Atoms)
        self.assertEquals(len(ase_molecule.get_atomic_numbers()), 11)

        ase_molecule = self.ts2.ase_molecule
        self.assertIsInstance(ase_molecule, Atoms)
        self.assertEquals(len(ase_molecule.get_atomic_numbers()), 11)
    
    def test_symmetry_number(self):
        self.assertEquals(self.ts.symmetry_number, 1)
        self.assertEquals(self.ts2.symmetry_number, 1)

    def test_bounds_matrix(self):

        lbl1 = self.ts.rmg_molecule.get_labeled_atoms("*1")[0].sorting_label
        lbl2 = self.ts.rmg_molecule.get_labeled_atoms("*2")[0].sorting_label
        lbl3 = self.ts.rmg_molecule.get_labeled_atoms("*3")[0].sorting_label

        d12 = self.ts.distance_data.distances["d12"]
        u12 = self.ts.distance_data.uncertainties["d12"]
        d13 = self.ts.distance_data.distances["d13"]
        u13 = self.ts.distance_data.uncertainties["d13"]
        d23 = self.ts.distance_data.distances["d23"]
        u23 = self.ts.distance_data.uncertainties["d23"]

        bm = self.ts.bm
        low, high = sorted([lbl1, lbl2])
        self.assertAlmostEqual(d12, bm[low, high], delta=u12/2 )
        low, high = sorted([lbl1, lbl3])
        self.assertAlmostEqual(d13, bm[low, high], delta=u13/2 )
        low, high = sorted([lbl2, lbl3])
        self.assertAlmostEqual(d23, bm[low, high], delta=u23/2 )


        lbl1 = self.ts2.rmg_molecule.get_labeled_atoms("*1")[0].sorting_label
        lbl2 = self.ts2.rmg_molecule.get_labeled_atoms("*2")[0].sorting_label
        lbl3 = self.ts2.rmg_molecule.get_labeled_atoms("*3")[0].sorting_label

        d12 = self.ts2.distance_data.distances["d12"]
        u12 = self.ts2.distance_data.uncertainties["d12"]
        d13 = self.ts2.distance_data.distances["d13"]
        u13 = self.ts2.distance_data.uncertainties["d13"]
        d23 = self.ts2.distance_data.distances["d23"]
        u23 = self.ts2.distance_data.uncertainties["d23"]

        bm = self.ts2.bm
        low, high = sorted([lbl1, lbl2])
        self.assertAlmostEqual(d12, bm[low, high], delta=u12/2 )
        low, high = sorted([lbl1, lbl3])
        self.assertAlmostEqual(d13, bm[low, high], delta=u13/2 )
        low, high = sorted([lbl2, lbl3])
        self.assertAlmostEqual(d23, bm[low, high], delta=u23/2 )
        
    def test_get_bonds(self):
        self.assertEqual(len(self.ts.get_bonds()), 10)
        true_count = 0
        for bond in self.ts.get_bonds():
            if bond.reaction_center: true_count += 1
        self.assertEqual(true_count, 2)

        self.assertEqual(len(self.ts2.get_bonds()), 10)
        true_count = 0
        for bond in self.ts2.get_bonds():
            if bond.reaction_center: true_count += 1
        self.assertEqual(true_count, 2)

    def test_get_angles(self):
        self.assertEqual(len(self.ts.get_angles()), 15)
        self.assertEqual(len(self.ts2.get_angles()), 15)

    def test_get_torsions(self):
        self.assertEqual(len(self.ts.get_torsions()), 4)
        self.assertEqual(len(self.ts2.get_torsions()), 4)


if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

