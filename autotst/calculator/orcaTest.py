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
import shutil

from autotst.calculator.orca import Orca
from autotst.species import Conformer


class TestOrca(unittest.TestCase):

    def setUp(self):
        conf = Conformer(smiles='C')
        self.orca = Orca(conformer=conf, directory=os.path.expandvars(
            "$AUTOTST/autotst/calculator/fod"))

    def test_load_conformer_attributes(self):
        charge = 0
        mult = 1
        label = 'CC(C)C'
        base = 'CC{C}C'
        self.orca.conformer = Conformer(smiles='CC(C)C')
        self.orca.load_conformer_attributes()
        self.assertTrue([charge,mult,label,base] == 
                        [self.orca.charge,self.orca.mult,self.orca.label,self.orca.base])

    def test_write_fod_input(self):
        if os.path.exists(self.orca.directory):
            shutil.rmtree(self.orca.directory)
        os.makedirs(self.orca.directory)
        self.orca.write_fod_input()
        self.assertTrue(os.path.exists(os.path.join(self.orca.directory,'C_fod.inp')))

    def test_check_normal_termination(self):
        path = os.path.expandvars(
            "$AUTOTST/test/bin/log-files/C_fod.log")
        self.assertTrue(self.orca.check_normal_termination(path))

    def test_read_fod_log(self):
        path = os.path.expandvars(
            "$AUTOTST/test/bin/log-files/C_fod.log")
        fod = self.orca.read_fod_log(path)
        self.assertEquals(float(0.000025),fod)

    def tearDown(self):

        if os.path.exists(os.path.expandvars("$AUTOTST/autotst/calculator/fod")):
            shutil.rmtree(os.path.expandvars(
                "$AUTOTST/autotst/calculator/fod"))

if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
