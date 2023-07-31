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
import numpy as np
import autotst
from autotst.common import AUTOTST_PATH, AUTOTST_DATABASE_PATH
from autotst.reaction import Reaction
from autotst.data.base import QMData, DistanceData, TransitionStates, TransitionStateDepository, TSGroups
import rmgpy
import rmgpy.data.rmg

class TestQMData(unittest.TestCase):
    def setUp(self):
        self.qmdata = QMData()
    def test_init(self):
        self.assertEqual(self.qmdata.ground_state_degeneracy, 0)
        self.assertEqual(self.qmdata.number_of_atoms, 0)
        self.assertEqual(self.qmdata.steric_energy, None)
        self.assertEqual(self.qmdata.molecular_mass[0], 0)
        self.assertEqual(self.qmdata.energy[0], 0)
        self.assertEqual(len(self.qmdata.atomic_numbers), 0)
        self.assertEqual(len(self.qmdata.rotational_constants[0]), 0)
        self.assertEqual(len(self.qmdata.atom_coords[0]), 1)
        self.assertEqual(len(self.qmdata.frequencies[0]), 0)
        self.assertEqual(self.qmdata.source, None)
        self.assertEqual(self.qmdata.method, None)
        self.assertEqual(self.qmdata.label, "")

    def test_get_qmdata(self):
        """
        A method that is designed to obtain the QM data for a transitionstate or molecule
        Returns a qmdata object
        """
        self.qmdata.get_qmdata(os.path.join(AUTOTST_PATH, "test", "bin", "log-files", "CC+[O]O_[CH2]C+OO_forward_0.log"))

        self.assertEqual(self.qmdata.ground_state_degeneracy, 2)
        self.assertAlmostEqual(self.qmdata.molecular_mass[0], 126.1, places=1)
        self.assertAlmostEqual(self.qmdata.energy[0], -6277.0, places=1)
        self.assertEqual(len(self.qmdata.atom_numbers), 11)
        self.assertEqual(self.qmdata.number_of_atoms, 11)
        self.assertEqual(len(self.qmdata.atom_coords[0]), 11)
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
        rmg_database = rmgpy.data.rmg.RMGDatabase()
        rmg_database.load(
            rmgpy.settings['database.directory'],
            kinetics_families=[
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration"
        ],
            transport_libraries=[],
            reaction_libraries=[],
            seed_mechanisms=[],
            thermo_libraries=[
                'primaryThermoLibrary',
                'thermo_DFT_CCSDTF12_BAC',
                'CBS_QB3_1dHR'],
            solvation=False,
        )
        self.rmg_database = rmg_database
        ts_database = TransitionStates()
        path = os.path.join(autotst.settings["tst_database_path"], "H_Abstraction")
        global_context = {'__builtins__': None}
        local_context = {'DistanceData': DistanceData}
        family = self.rmg_database.kinetics.families["H_Abstraction"]
        ts_database.family = family
        ts_database.load(path, local_context, global_context)
        self.ts_database = ts_database

    def test_load(self):

        ts_database = TransitionStates()
        path = os.path.join(autotst.settings["tst_database_path"], "H_Abstraction")
        global_context = {'__builtins__': None}
        local_context = {'DistanceData': DistanceData}
        family = self.rmg_database.kinetics.families["H_Abstraction"]
        ts_database.family = family
        ts_database.load(path, local_context, global_context)

        self.assertIsInstance(ts_database.depository, TransitionStateDepository)
        self.assertIsInstance(ts_database.groups, TSGroups)

    def test_estimate_distances(self):

        reaction = Reaction("CC+[O]O_[CH2]C+OO")
        labeled_reaction, _ = reaction.get_labeled_reaction()

        distance_data = self.ts_database.estimate_distances(labeled_reaction)

        d12 = 1.38
        d13 = 2.53
        # d13 in reactionTest is smaller than the distance in baseTest.py because
        # d13 is edited to be smaller in reaction.py. The value returned from the database
        # is ~2.53 but is reduced to ~2.43 when called from the reaction object itself
        d23 = 1.16

        self.assertAlmostEquals(d12, distance_data.distances["d12"], places=1)
        self.assertAlmostEquals(d13, distance_data.distances["d13"], places=1)
        self.assertAlmostEquals(d23, distance_data.distances["d23"], places=1)

class TestTransitionStateDepository(unittest.TestCase):

    def setUp(self):
        self.ts_depository = TransitionStateDepository(label="test")

        self.settings = {
            "file_path": os.path.join(
                os.path.join(AUTOTST_DATABASE_PATH, "database", "H_Abstraction", "TS_training", "reactions.py")
            ),
            "local_context": {"DistanceData":DistanceData},
            "global_context": {'__builtins__': None}
        }

    def test_load(self):

        self.ts_depository.load(
            self.settings["file_path"],
            self.settings["local_context"],
            self.settings["global_context"]
        )

class TestTSGroups(unittest.TestCase):

    def setUp(self):

        self.ts_groups = TSGroups(label="test")

        self.settings = {
            "file_path": os.path.join(
                os.path.join(AUTOTST_DATABASE_PATH, "database", "H_Abstraction", "TS_groups.py")
            ),
            "local_context": {"DistanceData":DistanceData},
            "global_context": {'__builtins__': None}
        }
    def test_load(self):

        self.ts_groups.load(
            self.settings["file_path"],
            self.settings["local_context"],
            self.settings["global_context"]
        )

    def test_estimate_distances_using_group_additivity(self):

        self.test_load()

        reaction = Reaction("CC+[O]O_[CH2]C+OO")
        labeled_reaction, _ = reaction.get_labeled_reaction()

        distance_data = self.ts_groups.estimate_distances_using_group_additivity(labeled_reaction)

        d12 = 1.38
        d13 = 2.53
        # d13 in reactionTest is smaller than the distance in baseTest.py because
        # d13 is edited to be smaller in reaction.py. The value returned from the database
        # is ~2.53 but is reduced to ~2.43 when called from the reaction object itself
        d23 = 1.16

        self.assertAlmostEquals(d12, distance_data.distances["d12"], places=1)
        self.assertAlmostEquals(d13, distance_data.distances["d13"], places=1)
        self.assertAlmostEquals(d23, distance_data.distances["d23"], places=1)

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))