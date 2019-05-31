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

import unittest, os, sys, shutil
from autotst.species import Conformer
from autotst.conformer.utilities import get_energy, find_terminal_torsions
from ase.calculators.emt import EMT

class TestUtilities(unittest.TestCase):

    def setUp(self):
        self.conformer = Conformer("NC(O)C=CC")
        self.conformer.get_molecules()

    def test_get_potential_energy(self):

        self.conformer.ase_molecule.set_calculator(EMT())
        test_energy = get_energy(self.conformer)
        energy = self.conformer.ase_molecule.get_potential_energy()
        self.assertAlmostEqual(test_energy, energy, places=2)

    def test_find_torsions(self):

        terminal_torsions, non_terminal_torsions = find_terminal_torsions(self.conformer)
        self.assertEqual(len(self.conformer.torsions), len(terminal_torsions + non_terminal_torsions))
        self.assertEqual(len(terminal_torsions), 1) # only one terminal methyl group
        self.assertEqual(len(non_terminal_torsions), 3) # N-C, O-C, C-C

    
if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))