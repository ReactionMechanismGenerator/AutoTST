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

from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
import unittest

class TestBond(unittest.TestCase):
    def setUp(self):
        self.bond = Bond(index=1,
                        atom_indices=[1,2],
                        length=1.8,
                        reaction_center=True,
                        mask=[True,True,True])
    def test_index(self):
        self.assertIsInstance(self.bond.index,int)
    def test_atom_indices(self):
        self.assertIsInstance(self.bond.atom_indices,list)
        self.assertEqual(len(self.bond.atom_indices),2)
    def test_length(self):
        self.assertIsInstance(self.bond.length,float)
    def test_reaction_center(self):
        self.assertIsInstance(self.bond.reaction_center,bool)
    def test_mask(self):
        self.assertIsInstance(self.bond.mask,list)

class TestAngle(unittest.TestCase):
    def setUp(self):
        self.angle = Angle(index=1,
                        atom_indices=[1,2],
                        degree=90.0,
                        reaction_center=True,
                        mask=[True,True,True])
    def test_index(self):
        self.assertIsInstance(self.angle.index,int)
    def test_atom_indices(self):
        self.assertIsInstance(self.angle.atom_indices,list)
        self.assertEqual(len(self.angle.atom_indices),2)
    def test_length(self):
        self.assertIsInstance(self.angle.degree,float)
    def test_reaction_center(self):
        self.assertIsInstance(self.angle.reaction_center,bool)
    def test_mask(self):
        self.assertIsInstance(self.angle.mask,list)

class TestTorsion(unittest.TestCase):
    def setUp(self):
        self.torsion = Torsion(index=1,
                        atom_indices=[1,2,3],
                        dihedral=60.0,
                        reaction_center=True,
                        mask=[True,True,True])
    def test_index(self):
        self.assertIsInstance(self.torsion.index,int)
    def test_atom_indices(self):
        self.assertIsInstance(self.torsion.atom_indices,list)
        self.assertEqual(len(self.torsion.atom_indices),3)
    def test_length(self):
        self.assertIsInstance(self.torsion.dihedral,float)
    def test_reaction_center(self):
        self.assertIsInstance(self.torsion.reaction_center,bool)
    def test_mask(self):
        self.assertIsInstance(self.torsion.mask,list)

class TestCisTrans(unittest.TestCase):
    def setUp(self):
        self.cistrans = CisTrans(index=1,
                        atom_indices=[1,2,3],
                        dihedral=60.0,
                        reaction_center=True,
                        stero='str',
                        mask=[True,True,True])
    def test_index(self):
        self.assertIsInstance(self.cistrans.index,int)
    def test_atom_indices(self):
        self.assertIsInstance(self.cistrans.atom_indices,list)
        self.assertEqual(len(self.cistrans.atom_indices),3)
    def test_length(self):
        self.assertIsInstance(self.cistrans.dihedral,float)
    def test_reaction_center(self):
        self.assertIsInstance(self.cistrans.reaction_center,bool)
    def test_mask(self):
        self.assertIsInstance(self.cistrans.mask,list)
    def test_stereo(self):
        self.assertIsInstance(self.cistrans.stero,str)

class TestChiralCenter(unittest.TestCase):
    def setUp(self):
        self.chiralcenter = ChiralCenter(index=1,
                        atom_index=1,
                        chirality='chirality')
    def test_index(self):
        self.assertIsInstance(self.chiralcenter.index,int)
    def test_atom_index(self):
        self.assertIsInstance(self.chiralcenter.atom_index,int)
    def test_chirality(self):
        self.assertIsInstance(self.chiralcenter.chirality,str)
        

if __name__ == "__main__":
  unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

