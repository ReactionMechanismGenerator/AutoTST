# To create TS geometry guesses

import numpy as np
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem.Pharm3D import EmbedLib

import ase

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

import os
from autotst.reaction import AutoTST_Reaction, AutoTST_TS
from autotst.molecule import AutoTST_Molecule
from autotst.geometry import Bond, Angle, Torsion
from ase.io.gaussian import read_gaussian_out

# To perform TS search
from ase.calculators.gaussian import Gaussian
from autotst.calculators.gaussian import AutoTST_Gaussian
from autotst.calculators.vibrational_analysis import Vibrational_Analysis, percent_change
from autotst.calculators.cantherm import AutoTST_CanTherm

# For conformer analysis
from ase.calculators.lj import LennardJones
from autotst.conformer.utilities  import create_initial_population, select_top_population
from autotst.conformer.ga import perform_ga
from autotst.conformer.simple_es import perform_simple_es

test_reaction = AutoTST_Reaction("[CH]=CC=C+[O]O_[CH]=C[C]=C+OO", "H_Abstraction")

scratch_dir = "./gaussian_example"

### Performing the partial optimizations
tst_calculators = AutoTST_Gaussian(test_reaction, scratch=scratch_dir, save_directory=".")

kinetics = tst_calculators.read_kinetics_file()

if kinetics:
    logging.info("We have previously loaded kinetics:")
    logging.info("{0!r}".format(kinetics['reaction']))

else:
    gaussian_results = tst_calculators.run_all(vibrational_analysis = False)
    if gaussian_results:
        ### Running CanTherm ###
        cantherm = AutoTST_CanTherm(tst_calculators.reaction, scratch=scratch_dir, output_directory=scratch_dir)
        cantherm.write_files()
        cantherm.run()
        cantherm.set_reactants_and_products()

        ### Printing Results ###
        logging.info("The kinetics of intrest are as follows:")
        logging.info("{0!r}".format(cantherm.kinetics_job.reaction))

        tst_calculators.save_kinetics(tst_calculators.method, cantherm.kinetics_job.reaction)

    else:
        logging.info("Failed gaussian... :(")
