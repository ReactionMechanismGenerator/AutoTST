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
import os
import logging
import shutil
import numpy as np
import autotst
from ..reaction import Reaction
from .base import QMData, DistanceData
from .inputoutput import InputOutput
import rmgpy.kinetics

class TestInputOutput(unittest.TestCase):

    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.io = InputOutput(
            reaction=self.reaction, 
            directory=os.path.expandvars("$AUTOTST/test/")
        )
        try:
            os.makedirs(os.path.join(
                os.path.expandvars("$AUTOTST/test/"),
                "ts",
                self.reaction.label
            ))
        except OSError:
            try:
                shutil.copy(
                    os.path.join(
                        os.path.expandvars("$AUTOTST/test/bin/log-files"), 
                        self.reaction.label + "_forward_0.log"
                        ),
                    os.path.join(
                        os.path.expandvars("$AUTOTST/test/"),
                        "ts",
                        self.reaction.label,
                        self.reaction.label + ".log"
                    )
                )
            except:
                pass
            

    def test_ts_file_path(self):
        path = self.io.get_ts_file_path()
        self.assertEqual(
            path,
            os.path.join(
                os.path.expandvars("$AUTOTST/test/"),
                "ts",
                self.reaction.label,
                self.reaction.label + ".ts"
            )
        )

    def test_kinetics_file_path(self):
        path = self.io.get_kinetics_file_path()
        self.assertEqual(
            path,
            os.path.join(
                os.path.expandvars("$AUTOTST/test/"),
                "ts",
                self.reaction.label,
                self.reaction.label + ".kinetics"
            )
        )
    
    def test_get_qmdata(self):

        self.qmdata = self.io.get_qmdata()
        self.assertEqual(self.qmdata.ground_state_degeneracy, 2)
        self.assertAlmostEqual(self.qmdata.molecular_mass[0], 126.1, places=1)
        self.assertAlmostEqual(self.qmdata.energy[0], -6277.0, places=1)
        self.assertEqual(len(self.qmdata.atom_numbers), 11)
        self.assertEqual(self.qmdata.number_of_atoms, 11)
        self.assertEqual(len(self.qmdata.atom_coords[0]), 11)
        self.assertEqual(len(self.qmdata.frequencies[0]), 27)
        self.assertEqual(self.qmdata.method.lower(), "m062x")
        self.assertEqual(self.qmdata.source, "AutoTST")

    def test_ts_io(self):
        self.assertTrue(self.io.save_ts())
        self.assertIsInstance(self.io.read_ts_file(), dict)

    def test_kinetics_io(self):

        self.io.reaction.get_rmg_reaction()
        self.io.reaction.rmg_reaction.kinetics = rmgpy.kinetics.Arrhenius()

        self.assertTrue(self.io.save_kinetics())
        self.assertIsInstance(self.io.read_kinetics_file(), dict)


if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))