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
from autotst.reaction import Reaction
from autotst.calculator.statmech import StatMech
from rmgpy.reaction import Reaction as RMGReaction
from rmgpy.kinetics import Arrhenius

class TestStatMech(unittest.TestCase):
    def setUp(self):
        self.reaction = Reaction("CC+[O]O_[CH2]C+OO")
        self.reaction.get_labeled_reaction()
        self.reaction.generate_reactants_and_products()

        directory = os.path.expandvars("$AUTOTST/test")
        if not os.path.exists(os.path.join(directory, "ts", self.reaction.label)):
            os.makedirs(os.path.join(directory, "ts", self.reaction.label))
        if not os.path.exists(os.path.join(directory, "ts", self.reaction.label, self.reaction.label + ".log")):
            shutil.copy(
                os.path.join(directory, "bin", "log-files", self.reaction.label + "_forward_0.log"),
                os.path.join(directory, "ts", self.reaction.label, self.reaction.label + ".log")
            )

        for sp in self.reaction.reactants + self.reaction.products:
            for smiles in sp.conformers.keys():
                if not os.path.exists(os.path.join(directory, "species", smiles)):
                    os.makedirs(os.path.join(directory, "species", smiles))
                if not os.path.exists(os.path.join(directory, "species", smiles, smiles+".log")):
                    shutil.copy(
                        os.path.join(directory, "bin", "log-files", smiles + "_0.log"),
                        os.path.join(directory, "species", smiles, smiles+".log")
                    )

        self.statmech = StatMech(
            reaction = self.reaction,
            directory = directory
        )

    def tearDown(self):
        directory = os.path.expandvars("$AUTOTST/test")
        if os.path.exists(os.path.join(directory, "ts")):
            shutil.rmtree(os.path.join(directory, "ts"))
        if os.path.exists(os.path.join(directory, "species")):
            shutil.rmtree(os.path.join(directory, "species"))

        for head, _, files in os.walk(os.path.expandvars("$AUTOTST")):
            for fi in files:
                if fi.endswith(".symm"):
                    os.remove(os.path.join(head, fi))

    def test_get_atoms(self):

        ts = self.reaction.ts["forward"][0]
        atom_dict = self.statmech.get_atoms(ts)
        self.assertEqual(atom_dict["H"], 7)
        self.assertEqual(atom_dict["C"], 2)
        self.assertEqual(atom_dict["O"], 2)
    
    def test_get_bonds(self):
        ts = self.reaction.ts["forward"][0]
        bond_dict = self.statmech.get_bonds(ts)
        self.assertEqual(bond_dict["C-C"], 1)
        self.assertTrue(bond_dict["C-H"] <= 6)

    def test_write_conformer_file(self):
        species = self.reaction.reactants[0]
        conformer = species.conformers.values()[0][0]
        self.assertTrue(self.statmech.write_conformer_file(conformer))

        self.assertTrue(
            os.path.exists(os.path.join(
                self.statmech.directory,
                "species",
                conformer.smiles,
                conformer.smiles + ".py"
            ))
        )
        # Running it again to see if it recognizes that a .py file was already written
        self.assertTrue(self.statmech.write_conformer_file(conformer))

    def test_write_species_file(self):
        species = self.reaction.reactants[0]

        self.statmech.write_species_files(species)
        for smiles in species.conformers.keys():
            self.assertTrue(
                os.path.exists(os.path.join(
                    self.statmech.directory,
                    "species",
                    smiles,
                    smiles + ".py"
                ))
            )
    def test_write_ts_input(self):
        ts = self.reaction.ts["forward"][0]
        self.assertTrue(self.statmech.write_ts_input(ts))
        self.assertTrue(os.path.exists(os.path.join(
            self.statmech.directory,
            "ts",
            self.reaction.label,
            self.reaction.label + ".py"
        )))
        self.assertTrue(self.statmech.write_ts_input(ts))

    def test_write_kinetics_input(self):
        self.statmech.write_kinetics_input()

        self.assertTrue(os.path.exists(os.path.join(
            self.statmech.directory,
            "ts",
            self.reaction.label,
            self.reaction.label + ".kinetics.py"
        )))

    def test_write_files(self):

        self.statmech.write_files()
        for mol in self.reaction.reactants + self.reaction.products:
            for confs in mol.conformers.values():
                conf = confs[0]
                self.assertTrue(os.path.exists(os.path.join(
                    self.statmech.directory,
                    "species", 
                    conf.smiles,
                    conf.smiles + ".py"
                )))
        self.assertTrue(os.path.exists(os.path.join(
            self.statmech.directory,
            "ts",
            self.reaction.label,
            self.reaction.label + ".py"
        )))
        self.assertTrue(os.path.exists(os.path.join(
            self.statmech.directory,
            "ts",
            self.reaction.label,
            self.reaction.label + ".kinetics.py"
        )))

    def test_run(self):
        self.statmech.write_files()
        self.statmech.run()
        self.assertIsInstance(self.statmech.kinetics_job.reaction, RMGReaction)
        self.assertIsInstance(self.statmech.kinetics_job.reaction.kinetics, Arrhenius)

    def test_set_results(self):
        self.test_run()
        self.statmech.set_results()
        self.assertTrue(
            self.reaction.rmg_reaction.isIsomorphic(
                self.statmech.kinetics_job.reaction
            )
        )
        self.assertTrue(
            self.statmech.reaction.rmg_reaction.isIsomorphic(
                self.statmech.kinetics_job.reaction
            )
        )
if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

    