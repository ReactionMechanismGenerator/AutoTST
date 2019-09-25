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

    def test_check_NormalTermination(self):
        path = os.path.expandvars(
            "$AUTOTST/test/bin/log-files/C_fod.log")
        self.assertTrue(self.orca.check_NormalTermination(path))

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
