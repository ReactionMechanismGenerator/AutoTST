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

import unittest, os, sys, shutil
from ..reaction import Reaction
from .update import get_unknown_species


class TestUpdateMethods(unittest.TestCase):

    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.reaction.get_labeled_reaction()
        self.family = "H_Abstraction"
        self.method = "m062x/cc-pVTZ"
        self.short_description = "test case"
    """
    def test_update(self):
        UpdateMethods.update_al(self.reaction, self.family, method=self.method, short_desc=self.short_description)
        self.assertTrue(True)"""

    def test_get_unknown_species(self):
        "A test designed to count the number of unknown species"
        unknowns = get_unknown_species([self.reaction], [])
        self.assertEqual(len(unknowns), 4)

        
if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))