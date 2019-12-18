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
from ..species import Conformer
from .utilities import get_energy, find_terminal_torsions
from .systematic import find_all_combos, systematic_search
import ase.calculators.emt

class TestSystematic(unittest.TestCase):

    def setUp(self):
        self.conformer = Conformer("CC(C)C=CC")
        self.conformer_3rad = Conformer("[C]=[CH]")
        self.conformer_4rad = Conformer("[C]=[C]")

    def test_find_all_combos(self):
        "Test that we can identify all possble conformers in a systematic search"
        combos = find_all_combos(self.conformer, delta=120, cistrans=True, chiral_centers=True)
        self.assertTrue(len(combos), 108)
        self.assertTrue(len(combos[0]), 3)
        
    def test_systematic_search(self):
        "Test that the systematic search can find more than 1 conformer"
        self.conformer.ase_molecule.set_calculator(ase.calculators.emt.EMT())
        confs = systematic_search(self.conformer, delta=180.0)
        self.assertTrue(1 < len(confs) <= 3)

    def test_systematic_search_multiplicity(self):
        self.conformer_3rad.ase_molecule.set_calculator(ase.calculators.emt.EMT())
        self.conformer_4rad.ase_molecule.set_calculator(ase.calculators.emt.EMT())
        confs_3rad = systematic_search(self.conformer_3rad, delta=180.0, energy_cutoff = "default",
                                      rmsd_cutoff = "default", multiplicity = True)
        confs_4rad = systematic_search(self.conformer_4rad, delta=180.0, energy_cutoff = "high",
                                      rmsd_cutoff = "loose", multiplicity = True)
        self.assertTrue(confs_3rad[0].rmg_molecule.multiplicity == 2)
        self.assertTrue(confs_3rad[1].rmg_molecule.multiplicity == 4)
        self.assertTrue(confs_4rad[0].rmg_molecule.multiplicity == 1)
        self.assertTrue(confs_4rad[1].rmg_molecule.multiplicity == 3)
        self.assertTrue(confs_4rad[2].rmg_molecule.multiplicity == 5)
        

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))