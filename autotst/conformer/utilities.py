#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
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
################################################################################

import os
import sys
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import itertools
import random
import numpy as np
from numpy import array
import pandas as pd
import cPickle as pickle


# do this before we have a chance to import openbabel!
import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction

from autotst.molecule import *
from autotst.reaction import *

from ase.calculators.morse import *
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *

def create_initial_population(autotst_object, delta=30, population_size=200):
    """
    A function designed to take a multi_molecule, multi_rxn or multi_ts object
    and create an initial population of conformers.

    :param:
     multi_object: a multi_molecule, multi_ts, or multi_rxn object
     calc: an ASE calculator. If none is chosen, an EMT() calculator will be used
     delta: the step size between possible dihedral angles in degrees.
     population_size: the number of individuals to be used for the population

    :return:
     df: a DataFrame of the results sorted by the lowest energy conformers
    """
    logging.info("Creating initial population of {} individuals from random guesses".format(population_size))

    df = None

    def update(autotst_obj):
        if isinstance(autotst_object, autotst.molecule.AutoTST_Molecule):
            autotst_obj.update_from_ase_mol()

        if isinstance(autotst_object, autotst.reaction.AutoTST_Reaction):
            autotst_obj.ts.update_from_ase_ts()

        if isinstance(autotst_object, autotst.reaction.AutoTST_TS):
            autotst_obj.update_from_ase_ts()



    possible_dihedrals = np.arange(0, 360, delta)
    population = []
    if isinstance(autotst_object, autotst.molecule.AutoTST_Molecule):
        logging.info("The object given is a `AutoTST_Molecule` object")

        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_molecule
        rdkit_object = autotst_object.rdkit_molecule

    if isinstance(autotst_object, autotst.reaction.AutoTST_Reaction):
        logging.info("The object given is a `AutoTST_Reaction` object")
        torsions = autotst_object.ts.torsions
        ase_object = autotst_object.ts.ase_ts
        ase_object = autotst_object.ts.rdkit_ts

    if isinstance(autotst_object, autotst.reaction.AutoTST_TS):
        logging.info("The object given is a `AutoTST_TS` object")
        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_ts
        ase_object = autotst_object.rdkit_ts


    for indivudual in range(population_size):
        dihedrals = []
        for torsion in torsions:
            dihedral = np.random.choice(possible_dihedrals)
            dihedrals.append(dihedral)
            i, j, k, l = torsion.indices
            right_mask = torsion.right_mask

            ase_object.set_dihedral(
                a1=i,
                a2=j,
                a3=k,
                a4=l,
                angle=float(dihedral),
                mask=right_mask
            )

        update(autotst_object)

        ff1 = UFFGetMoleculeForceField(rdkit_object)
        constrained_energy = ff1.CalcEnergy()

        rd_copy = rdkit_object.__copy__()
        UFFOptimizeMolecule(rd_copy)
        ff2 = UFFGetMoleculeForceField(rd_copy)
        relaxed_energy = ff2.CalcEnergy()

        relaxed_torsions = []

        for torsion in torsions:

            i, j, k, l = torsion.indices
            angle = round(rdkit.Chem.rdMolTransforms.GetDihedralDeg(rd_copy.GetConformer(0), i,j,k,l), -1)
            if angle < 0: angle += 360
            relaxed_torsions.append(angle)

        population.append([constrained_energy, relaxed_energy] + dihedrals + relaxed_torsions)


    if len(population) > 0:
        logging.info("Creating a dataframe of the initial population")
        df = pd.DataFrame(population)
        columns = ["constrained_energy", "relaxed_energy"]
        for i in range(len(torsions)):
            columns = columns + ["torsion_" + str(i)]

        for i in range(len(torsions)):
            columns = columns + ["relaxed_torsion_" + str(i)]
        df.columns = columns
        df = df.sort_values("constrained_energy")

    return df


def select_top_population(df=None, top_percent=0.30):
    """
    :param:
     df: a DataFrame of a population of torsions with columns of `Energy` and `Torsion N`
     top_percent: a float of the top percentage of the population requested

    :return:
     top: a DataFrame containing the top percentage of the population
    """

    logging.info("Selecting the top population")
    assert "constrained_energy" in df.columns
    df.sort_values("constrained_energy")
    population_size = df.shape[0]
    top_population = population_size * top_percent
    top = df.iloc[:int(top_population), :]

    return top
