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

from autotst.molecule import *
from autotst.reaction import *
from autotst.conformer.utilities import *

from ase.calculators.morse import *
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *


def perform_simple_es(multi_object,
                      df,
                      top_percent=0.3,
                      tolerance=1e-4,
                      max_generations=np.inf,
                      store_generations=False,
                      store_directory="."):
    """
    Performs a simple evolutionary strategy to determine the lowest energy conformer of a TS or molecule

    :param multi_object: a multi_ts, multi_rxn, or multi_molecule that you want to perform conformer analysis on
       * the ase_object of the multi_object must have a calculator attached to it.
    :param df: a DataFrame containing the initial population
    :param top_percent: float of the top percentage of conformers you want to select
    :param tolerance: float of one of the possible cut off points for the analysis
    :param max_generations: int of one of the possible cut off points for the analysis
    :param store_generations: do you want to store pickle files of each generation
    :param store_directory: the director where you want the pickle files stored
    :return df: a DataFrame containing the final generation
    """
    top = select_top_population(df,
                                top_percent=top_percent
                                )
    top_population = top.shape[0]

    population_size = df.shape[0]

    # Takes each of the molecule objects
    if "AutoTST_Molecule" in str(multi_object.__class__):
        ase_object = multi_object.ase_molecule
        torsions = multi_object.torsions

    elif "AutoTST_Reaction" in str(multi_object.__class__):
        ase_object = multi_object.ts.ase_ts
        torsions = multi_object.ts.torsions

    elif "AutoTST_TS" in str(multi_object.__class__):
        ase_object = multi_object.ase_ts
        torsions = multi_object.torsions

    gen_number = 0
    complete = False
    while complete == False:
        gen_number += 1
        logging.info("Performing Simple ES on generation {}".format(gen_number))

        results = []
        for individual in range(population_size):
            dihedrals = []
            for index, torsion in enumerate(torsions):
                i, j, k, l = torsion.indices
                right_mask = torsion.right_mask

                dihedral = random.gauss(top.mean()[index + 1], top.std()[index + 1])
                dihedrals.append(dihedral)
                ase_object.set_dihedral(a1=i,
                                        a2=j,
                                        a3=k,
                                        a4=l,
                                        angle=float(dihedral),
                                        mask=right_mask)

                # Updating the molecule
                if "AutoTST_Molecule" in str(multi_object.__class__):
                    multi_object.update_from_ase_mol()

                elif "AutoTST_Reaction" in str(multi_object.__class__):
                    multi_object.ts.update_from_ase_ts()

                elif "AutoTST_TS" in str(multi_object.__class__):
                    multi_object.update_from_ase_ts()

            e = ase_object.get_potential_energy()
            results.append([e] + dihedrals)

        df = pd.DataFrame(results)
        logging.info("Creating the DataFrame of results for the {}th generation".format(gen_number))

        columns = ["Energy"]
        for i in range(len(torsions)):
            columns = columns + ["Torsion " + str(i)]

        df.columns = columns
        df = df.sort_values("Energy")

        if store_generations == True:
            # This portion stores each generation if desired

            if "AutoTST_Reaction" in str(multi_object.__class__):
                generation_name = "rxn_simple_es_generation_{}.csv".format(gen_number)

            elif "AutoTST_Molecule" in str(multi_object.__class__):
                generation_name = "mol_simple_es_generation_{}.csv".format(gen_number)

            elif "AutoTST_TS" in str(multi_object.__class__):
                generation_name = "ts_simple_es_generation_{}.csv".format(gen_number)

            else:
                generation_name = "simple_es_generation_{}.csv".format(gen_number)

            f = os.path.join(store_directory, generation_name)
            df.to_csv(f)

        top = df.iloc[:int(top_population), :]

        if gen_number >= max_generations:
            complete = True
            logging.info("Max generations reached. Simple ES complete.")
        if top.std()[0] / top.mean()[0] < tolerance:
            complete = True
            logging.info("Cutoff criteria reached. Simple ES complete.")

    return df
