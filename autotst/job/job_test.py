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

import os, sys, subprocess, shutil
import unittest
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.data.base import TransitionStates
from autotst.job.job import Job
from autotst.calculator.gaussian import Gaussian
import multiprocessing
import subprocess
import time

class JobTest(unittest.TestCase):

    def setUp(self):
        os.environ["PATH"] = os.path.expandvars("$AUTOTST/test/bin:") + os.environ["PATH"]
        os.environ["TEST_STATUS"] = "None"
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.calculator = Gaussian(directory=os.path.expandvars("$AUTOTST/test"))
        self.job = Job(
            reaction=self.reaction,
            calculator=self.calculator,
            partition="test",
            username="test",
            exclude="test",
            account="test"
        )
        self.job2 = Job(
            reaction=self.reaction,
            calculator=self.calculator,
            partition="test",
            username="test",
            exclude="test",
            account=["test"]
        )

    def test_setup(self):
        self.assertEqual(self.job.username, "test")
        self.assertEqual(self.job.exclude, "test")
        self.assertEqual(self.job.partition, "test")
        self.assertEqual(self.job.account, "test")
        self.assertEqual(self.job.label, self.reaction.label)

    def test_setup2(self):
        job = Job(directory=".")
        self.assertEqual(job.directory, ".")
        self.assertEqual(job.scratch, ".")

    def test_read_log(self):

        path = os.path.expandvars("$AUTOTST/test/bin/log-files/CC_0.log")

        atoms = self.job.read_log(path)

        self.assertEqual(len(atoms), 8)

        carbon_count = 0
        hydrogen_count = 0
        for atom in atoms:
            if atom.symbol == "H":
                hydrogen_count += 1
            elif atom.symbol == "C":
                carbon_count += 1

        self.assertEqual(hydrogen_count, 6)
        self.assertEqual(carbon_count, 2)

    def test_write_input(self):
        self.assertTrue(True)

    def test_check_complete(self):
        ### I don't know how to create alaises in a python script
        os.environ["TEST_STATUS"] = "None"
        self.assertFalse(self.job.check_complete("test1"))
        self.assertTrue(self.job.check_complete("test2"))

    def test_submit(self):

        result = self.job.submit("echo testing")
        self.assertTrue(result)

    ### For conformers
    def test_submit_conformer(self):
        self.reaction.generate_reactants_and_products()
        conformer = list(self.reaction.reactants[0].conformers.values())[0][0]
        label = self.job.submit_conformer(conformer)
        self.assertEqual(label, f"{conformer.smiles}_{conformer.index}")

    def test_submit_conformer2(self):
        self.reaction.generate_reactants_and_products()
        conformer = list(self.reaction.reactants[0].conformers.values())[0][0]
        label = self.job2.submit_conformer(conformer)
        self.assertEqual(label, f"{conformer.smiles}_{conformer.index}")

    def test_calculate_conformer(self):
        conformer = Conformer(smiles='CC',index=0)
        result = self.job.calculate_conformer(conformer=conformer)
        self.assertTrue(result)

    def test_calculate_species(self):
        self.reaction.generate_reactants_and_products()

        for species in self.reaction.reactants + self.reaction.products:
            self.job.calculate_species(species)
            for smiles in species.conformers.keys():
                self.assertTrue(os.path.exists(os.path.join(
                    os.path.expandvars("$AUTOTST/test/species/"),
                    smiles,
                    smiles + ".log"
                )))
    
    def test_submit_transitionstate(self):
        ts = self.reaction.ts["forward"][0]
        ts.get_molecules()
        for opt_type in ["shell", "center", "overall"]:
            label = self.job.submit_transitionstate(ts, opt_type=opt_type)
            if opt_type == "overall":
                self.assertEqual(label, f"{ts.reaction_label}_{ts.direction}_{ts.index}")
            else:
                self.assertEqual(label, f"{ts.reaction_label}_{ts.direction}_{opt_type}_{ts.index}")

    def test_submit_transitionstate2(self):
        ts = self.reaction.ts["forward"][0]
        ts.get_molecules()
        for opt_type in ["shell", "center", "overall"]:
            label = self.job2.submit_transitionstate(ts, opt_type=opt_type)
            if opt_type == "overall":
                self.assertEqual(label, f"{ts.reaction_label}_{ts.direction}_{ts.index}")
            else:
                self.assertEqual(label, f"{ts.reaction_label}_{ts.direction}_{opt_type}_{ts.index}")

    def test_calculate_transitionstate(self):
        ts = self.reaction.ts["forward"][0]
        ts.get_molecules()
        result = self.job.calculate_transitionstate(ts)
        self.assertTrue(result)

    def test_calculate_reaction(self):

        del self.reaction.ts["reverse"]
        result = self.job.calculate_reaction()
        self.assertTrue(result)
    
if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))


