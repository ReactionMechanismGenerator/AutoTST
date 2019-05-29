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

from autotst.data.base import * #QMData #, DistanceData, TransitionStates, TransitionStateDepository, TSGroups

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
        self.assertEqual(len(self.qmdata.atomCoords[0]), 0)
        self.assertEqual(len(self.qmdata.frequencies[0]), 0)
        self.assertEqual(self.qmdata.source, None)
        self.assertEqual(self.qmdata.method, None)
        self.assertEqual(self.qmdata.label, "")

    


if __name__ == "__main__":
  unittest.main()