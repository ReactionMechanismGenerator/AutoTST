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

import os, shutil
import itertools
import logging
import numpy as np
from cclib.io import ccread

import autotst
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.geometry import Torsion
from autotst.calculator.gaussian import Gaussian

from cclib.io import ccread

from ase import Atom, Atoms
from ase.io.gaussian import read_gaussian, read_gaussian_out
from ase.calculators.gaussian import Gaussian as ASEGaussian

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.reaction import Reaction as RMGReaction

class TestGaussian(unittest.TestCase):
    def setUp(self):
        os.environ["PATH"] = os.path.expandvars("$AUTOTST/test/bin:") + os.environ["PATH"]
        rxn = Reaction(label='C+[O]O_[CH3]+OO')
        ts = rxn.ts["forward"][0]
        ts.get_molecules()
        self.gaussian = Gaussian(conformer=ts)
    def test_rotor_calc(self):
        autotst_gaussian_rotor = self.gaussian.get_rotor_calc()
        calc_dict = autotst_gaussian_rotor.todict()
        self.assertIsInstance(autotst_gaussian_rotor,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)
        """
        default_settings = {
            "method": "m062x",
            "basis": "cc-pVTZ",
            "mem": "5GB",
            "nprocshared": 20,
        }
        default_settings['multiplicity'] = 
        self.assertDictEqual(calc_dict,default_settings)
        """
    def test_conformer_calc(self):
        autotst_gaussian_confomer = Gaussian(conformer=Conformer(smiles='CCC')).get_conformer_calc()
        calc_dict = autotst_gaussian_confomer.todict()
        self.assertIsInstance(autotst_gaussian_confomer,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)

        #self.assertEqual(calc_dict['extra'],"opt=(calcfc,maxcycles=900) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900)")

    def test_shell_calc(self):
        autotst_gaussian_shell = self.gaussian.get_shell_calc()
        calc_dict = autotst_gaussian_shell.todict()
        self.assertIsInstance(autotst_gaussian_shell,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)
    
    def test_center_calc(self):
        autotst_gaussian_center = self.gaussian.get_center_calc()
        calc_dict = autotst_gaussian_center.todict()
        self.assertIsInstance(autotst_gaussian_center,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)
    
    def test_overall_calc(self):
        autotst_gaussian_overall = self.gaussian.get_overall_calc()
        calc_dict = autotst_gaussian_overall.todict()
        self.assertIsInstance(autotst_gaussian_overall,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)

    def test_irc_calc(self):
        autotst_gaussian_irc = self.gaussian.get_irc_calc()
        calc_dict = autotst_gaussian_irc.todict()
        self.assertIsInstance(autotst_gaussian_irc,ASEGaussian)
        self.assertIsInstance(calc_dict,dict)

    def tearDown(self):

        if os.path.exists(os.path.expandvars("$AUTOTST/autotst/calculator/ts")):
            shutil.rmtree(os.path.expandvars("$AUTOTST/autotst/calculator/ts"))
        if os.path.exists(os.path.expandvars("$AUTOTST/autotst/calculator/species")):
            shutil.rmtree(os.path.expandvars("$AUTOTST/autotst/calculator/species"))

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    
    