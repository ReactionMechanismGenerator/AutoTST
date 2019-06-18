import unittest

import os
import shutil

from autotst.calculator.orca import Orca
from autotst.species import Conformer


class TestOrca(unittest.TestCase):

    def setUp(self):
        conf = Conformer(smiles='C')
        self.orca = Orca(conformer=conf)

    def test_write_fod_input(self):
        path = os.path.expandvars(
            "$AUTOTST/autotst/calculator/fod")
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)
        self.orca.write_fod_input(path=path)
        self.assertTrue(os.path.exists(os.path.join(path,'C_fod.inp')))

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
