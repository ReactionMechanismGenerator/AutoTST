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
import os
import logging
import numpy as np
import rmgpy
from rmgpy.data.rmg import RMGDatabase
import autotst
from autotst.data.base import QMData, DistanceData, TransitionStates, TransitionStateDepository, TSGroups

class TestQMData(unittest.TestCase):
    def setUp(self):
        self.qmdata = QMData()
    def test_init(self):
        self.assertEqual(self.qmdata.groundStateDegeneracy, 0)
        self.assertEqual(self.qmdata.numberOfAtoms, 0)
        self.assertEqual(self.qmdata.stericEnergy, None)
        self.assertEqual(self.qmdata.molecularMass[0], 0)
        self.assertEqual(self.qmdata.energy[0], 0)
        self.assertEqual(len(self.qmdata.atomicNumbers), 0)
        self.assertEqual(len(self.qmdata.rotationalConstants[0]), 0)
        self.assertEqual(len(self.qmdata.atomCoords[0]), 1)
        self.assertEqual(len(self.qmdata.frequencies[0]), 0)
        self.assertEqual(self.qmdata.source, None)
        self.assertEqual(self.qmdata.method, None)
        self.assertEqual(self.qmdata.label, "")

    def test_get_qmdata(self):
        """
        A method that is designed to obtain the QM data for a transitionstate or molecule
        Returns a qmdata object
        """
        self.qmdata.get_qmdata(os.path.expandvars("$AUTOTST/test_scratch/ts/CC+[O]O_[CH2]C+OO/CC+[O]O_[CH2]C+OO.log"))

        self.assertEqual(self.qmdata.groundStateDegeneracy, 2)
        self.assertAlmostEqual(self.qmdata.molecularMass[0], 126.1, places=1)
        self.assertAlmostEqual(self.qmdata.energy[0], -6277.0, places=1)
        self.assertEqual(len(self.qmdata.atomNumbers), 11)
        self.assertEqual(self.qmdata.numberOfAtoms, 11)
        self.assertEqual(len(self.qmdata.atomCoords[0]), 11)
        self.assertEqual(len(self.qmdata.frequencies[0]), 27)
        self.assertEqual(self.qmdata.method.lower(), "m062x")
        self.assertEqual(self.qmdata.source, "AutoTST")


class TestDistanceData(unittest.TestCase):

    def setUp(self):
        self.distancedata = DistanceData(
            distances = {"d12": 1.0, "d13": 1.0, "d23": 1.0},
            uncertainties = {"d12": 0.1, "d13": 0.1, "d23": 0.1}
        )

    def test_add(self):
        self.distancedata.add(self.distancedata)
        for key in self.distancedata.distances.keys():
            self.assertEqual(self.distancedata.distances[key], 2.0)
            self.assertEqual(self.distancedata.uncertainties[key], 0.2)

class TestTransitionStates(unittest.TestCase):
    def setUp(self):
        rmg_database = RMGDatabase()
        rmg_database.load(
            rmgpy.settings['database.directory'],
            kineticsFamilies=[
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration"
        ],
            transportLibraries=[],
            reactionLibraries=[],
            seedMechanisms=[],
            thermoLibraries=[
                'primaryThermoLibrary',
                'thermo_DFT_CCSDTF12_BAC',
                'CBS_QB3_1dHR'],
            solvation=False,
        )
        self.rmg_database = rmg_database

    def test_load(self):

        ts_database = TransitionStates()
        path = os.path.join(autotst.settings["tst_database_path"], "H_Abstraction")
        global_context = {'__builtins__': None}
        local_context = {'DistanceData': DistanceData}
        family = self.rmg_database.kinetics.families["H_Abstraction"]
        ts_database.family = family
        ts_database.load(path, local_context, global_context)
        


if __name__ == "__main__":
    unittest.main()