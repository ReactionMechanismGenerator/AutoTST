import os
import sys
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import re
import imp
import itertools
import random
import numpy as np
from numpy import array
import pandas as pd


# do this before we have a chance to import openbabel!
import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

#import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction

import ase
from ase import Atom, Atoms

from autotst.molecule import *
from autotst.reaction import *
from autotst.conformer.utilities import *

from ase.calculators.morse import * #chosing this calculator for now because it's fast
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *

import cPickle as pickle


def perform_brute_force(multi_object,
                        delta=float(30),
                        store_directory="."):

    if "AutoTST_Molecule" in str(multi_object.__class__):
        ase_object = multi_object.ase_molecule
        torsions = multi_object.torsions
        file_name = "mol_brute_force.csv"

    elif "AutoTST_Reaction" in str(multi_object.__class__):
        ase_object = multi_object.ts.ase_ts
        torsions = multi_object.ts.torsions
        file_name = "rxn_brute_force.csv"

    elif "AutoTST_TS" in str(multi_object.__class__):
        ase_object = multi_object.ase_ts
        torsions = multi_object.torsions
        file_name = "ts_brute_force.csv"

    torsion_angles = np.arange(0, 360, delta)
    torsion_combos = list(itertools.combinations_with_replacement(torsion_angles, len(torsions)))
    if len(torsions) != 1:
        torsion_combos = list(
            set(
                torsion_combos +
                list(itertools.combinations_with_replacement(
                    torsion_angles[::-1], len(torsions)
                ))))

    results = []
    for index, combo in enumerate(torsion_combos):
        logging.info("Generating conformer {}".format(index))
        logging.info("Applying the torsion combo {0} to the molecule or TS.".format(combo))
        geo = zip(torsions, combo)
        for torsion in geo:
            tor = torsion[0]
            angle = torsion[1]

            i, j, k, l = tor.indices
            right_mask = tor.right_mask
            ase_object.set_dihedral(a1=i,
                                    a2=j,
                                    a3=k,
                                    a4=l,
                                    angle=float(angle),
                                    mask=right_mask)

        results.append([ase_object.get_potential_energy()] + list(combo))
    brute_force = pd.DataFrame(results)
    columns = ["Energy"]
    for i in range(len(torsions)):
        columns = columns + ["Torsion " + str(i)]

    brute_force.columns = columns

    brute_force = brute_force.sort_values("Energy")

    f = os.path.join(store_directory,filename)
    df.to_csv(f)

    return brute_force
