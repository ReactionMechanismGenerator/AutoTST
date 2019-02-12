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
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.species import Species as RMGSpecies
from rmgpy.reaction import Reaction as RMGReaction, ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

import os
from autotst.reaction import Reaction, TS
from autotst.molecule import Molecule
from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
from ase.io.gaussian import read_gaussian_out

# To perform TS search
from autotst.calculators.calculator import Calculator
from ase.calculators.gaussian import Gaussian as ASEGaussian
from autotst.calculators.gaussian import Gaussian
from autotst.calculators.vibrational_analysis import VibrationalAnalysis, percent_change
from autotst.calculators.statmech import StatMech

class 

class Job():
    """
    A class to deal with the input and output of calculations
    """

    def __init__(self, reaction=None):

        self.reaction = reaction
        if isinstance(reaction, Reaction):
            self.label = reaction.label

    def __repr__(self):
        return "< Job '{}'>".format(label)


    def submit(self, calculator=None):
        return None

    def 
    

    