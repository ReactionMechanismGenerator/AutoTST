#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2020 Richard H. West (r.west@northeastern.edu)
#   and the AutoTST Team
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
import os, shutil, logging
import pandas as pd
import numpy as np
import ase
import cclib.io
import rmgpy
import rmgpy.molecule
from ..reaction import Reaction, TS
from ..species import Species, Conformer
from .vibrational_analysis import percent_change, VibrationalAnalysis

class VibrationalAnalysisTest(unittest.TestCase):

    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.reaction.get_labeled_reaction()
        self.ts = self.reaction.ts["forward"][0]
        self.ts.get_molecules()

        directory = os.path.expandvars("$AUTOTST/test")
        if not os.path.exists(os.path.join(directory, "ts", self.reaction.label, "conformers")):
            os.makedirs(os.path.join(directory, "ts", self.reaction.label, "conformers"))
        if not os.path.exists(os.path.join(directory, "ts", self.reaction.label, self.reaction.label + ".log")):
            shutil.copy(
                os.path.join(directory, "bin", "log-files", self.reaction.label + "_forward_0.log"),
                os.path.join(directory, "ts", self.reaction.label, "conformers", self.reaction.label + "_forward_0.log")
            )

        self.directory = directory
        self.vibrational_analysis = VibrationalAnalysis(
            transitionstate = self.ts,
            directory = self.directory
        )

        self.reaction2 = Reaction("[CH3]+CC(F)(F)F_C+[CH2]C(F)(F)F")
        self.reaction2.get_labeled_reaction()
        self.ts2 = self.reaction2.ts["forward"][0]
        log_file = os.path.join(
            directory, "bin", "log-files", "[CH3]+CC(F)(F)F_C+[CH2]C(F)(F)F.log")
        self.vibrational_analysis2 = VibrationalAnalysis(
            transitionstate = self.ts2,
            log_file = log_file
        )

    def test_get_log_file(self):
        log_file = self.vibrational_analysis.get_log_file()

        actual_path = os.path.join(
            self.directory,
            "ts",
            self.ts.reaction_label,
            "conformers",
            f"{self.ts.reaction_label}_{self.ts.direction}_{self.ts.index}.log"
        )

        self.assertEqual(log_file, actual_path)

    def test_parse_vibrations(self):
        print(self.directory)
        print(self.vibrational_analysis.directory)
        vibrations = self.vibrational_analysis.parse_vibrations()
        self.assertIsInstance(vibrations, list)
        self.assertEqual(len(vibrations), 27)

    def test_obtain_geometries(self):

        vibrations = self.vibrational_analysis.parse_vibrations()

        symbol_dict = {
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }
        atoms = []

        parser = cclib.io.ccread(self.vibrational_analysis.log_file, loglevel=logging.ERROR)

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
            atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))

        test_pre_geometry = ase.Atoms(atoms)
        test_post_geometry = test_pre_geometry.copy()

        for vib, displacement in vibrations:
            if vib < 0:
                test_post_geometry.arrays["positions"] -= displacement

        pre_geometry, post_geometry = self.vibrational_analysis.obtain_geometries()

        for i, positions in enumerate(test_pre_geometry.arrays["positions"]):
            for j, x in enumerate(positions):
                self.assertEqual(x, pre_geometry.arrays["positions"][i][j])
        for i, positions in enumerate(test_post_geometry.arrays["positions"]):
            for j, x in enumerate(positions):
                self.assertEqual(x, post_geometry.arrays["positions"][i][j])

    def test_validate_ts(self):
        self.assertTrue(self.vibrational_analysis.validate_ts())
    
    def test_validate_by_connecting_the_dots(self):
        self.assertTrue(self.vibrational_analysis2.validate_by_connecting_the_dots())

    def test_validate(self):
        self.assertTrue(self.vibrational_analysis2.validate())

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    try:
        directory = os.path.expandvars("$AUTOTST/test")
        if os.path.exists(os.path.join(directory, "ts")):
            shutil.rmtree(os.path.join(directory, "ts"))

        for head, _, files in os.walk(os.path.expandvars("$AUTOTST")):
            for fi in files:
                if fi.endswith(".symm"):
                    os.remove(os.path.join(head, fi))
    except:
        None

