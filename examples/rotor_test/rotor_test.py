#Import a variety of things
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
import sys
from autotst.reaction import AutoTST_Reaction, AutoTST_TS
from autotst.molecule import AutoTST_Molecule
from autotst.geometry import Bond, Angle, Torsion
from ase.io.gaussian import read_gaussian_out

# To perform TS search
from ase.calculators.gaussian import Gaussian
from autotst.calculators.gaussian import AutoTST_Gaussian
from autotst.calculators.vibrational_analysis import Vibrational_Analysis, percent_change
from autotst.calculators.cantherm import AutoTST_CanTherm

if len(sys.argv)>1:
    i = int(sys.argv[-1])
elif os.getenv('SLURM_ARRAY_TASK_ID'):
    i = int(os.getenv('SLURM_ARRAY_TASK_ID'))
elif os.getenv('LSB_JOBINDEX'):
    i = int(os.getenv('LSB_JOBINDEX'))
else:
    #raise Exception("Specify a TS number!")
    logging.warning("Number not specified as script argument or via environment variable, so using default")
    i = 1
logging.info("RUNNING WITH JOB NUMBER i = {}".format(i))

mol = AutoTST_Molecule("C(C)C(CC)CC")
calcs = AutoTST_Gaussian().get_rotor_calcs(mol)

tor = mol.torsions[i-1]

logging.info("Starting Calculations for torsion {0} of {1}".format(tor.indices, mol.smiles))

try:
    calc = calcs[tor.indices[1:3]]
    calc.calculate(mol.ase_molecule)
except:
    logging.info("Couldn't calculate the following torsion: {}".format(tor.indices))

