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

import autotst
from autotst.geometry import Bond, Angle, Torsion, CisTrans
from autotst.molecule import AutoTST_Molecule
from autotst.reaction import AutoTST_Reaction, AutoTST_TS
from autotst.conformer.utilities import update_from_ase, create_initial_population, \
                                        select_top_population, get_unique_conformers, get_energies


def perform_ga(autotst_object,
               initial_pop=None,
               top_percent=0.3,
               tolerance=0.0001,
               max_generations=500,
               store_generations=False,
               store_directory=".",
               mutation_probability=0.2,
               delta=30):

    """
    Performs a genetic algorithm to determine the lowest energy conformer of a TS or molecule. 

    :param autotst_object: am autotst_ts, autotst_rxn, or autotst_molecule that you want to perform conformer analysis on
       * the ase_object of the autotst_object must have a calculator attached to it.
    :param initial_pop: a DataFrame containing the initial population
    :param top_percent: float of the top percentage of conformers you want to select
    :param tolerance: float of one of the possible cut off points for the analysis
    :param max_generations: int of one of the possible cut off points for the analysis
    :param store_generations: do you want to store pickle files of each generation
    :param store_directory: the director where you want the pickle files stored
    :param mutation_probability: float of the chance of mutation
    :param delta: the degree change in dihedral angle between each possible dihedral angle

    :return results: a DataFrame containing the final generation
    :return unique_conformers: a dictionary with indicies of unique torsion combinations and entries of energy of those torsions
    """
    assert autotst_object, "No AutoTST object provided..."
    if initial_pop is None:
        logging.info("No initial population provided, creating one using base parameters...")
        initial_pop = create_initial_population(autotst_object)

    possible_dihedrals = np.arange(0, 360, delta)
    top = select_top_population(initial_pop,
                                top_percent=top_percent
                                )

    population_size = initial_pop.shape[0]

    results = initial_pop

    if isinstance(autotst_object, autotst.molecule.AutoTST_Molecule):
        logging.info("The object given is a `AutoTST_Molecule` object")
        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_molecule
        label = autotst_object.smiles

    if isinstance(autotst_object, autotst.reaction.AutoTST_Reaction):
        logging.info("The object given is a `AutoTST_Reaction` object")
        torsions = autotst_object.ts.torsions
        ase_object = autotst_object.ts.ase_ts
        label = autotst_object.label

    if isinstance(autotst_object, autotst.reaction.AutoTST_TS):
        logging.info("The object given is a `AutoTST_TS` object")
        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_ts
        label = autotst_object.label

    assert ase_object.get_calculator(), "To use GA, you must attach an ASE calculator to the `ase_molecule`."
    gen_number = 0
    complete = False
    unique_conformers = {}
    while complete == False:
        gen_number += 1
        logging.info("Performing GA on generation {}".format(gen_number))

        r = []
        for individual in range(population_size):
            parent_0, parent_1 = random.sample(top.index, 2)
            dihedrals = []
            for index, torsion in enumerate(torsions):

                if random.random() < mutation_probability:
                    dihedral = np.random.choice(possible_dihedrals)
                else:
                    if 0.5 > random.random():
                        dihedral = results["torsion_" + str(index)].loc[parent_0]
                    else:
                        dihedral = results["torsion_" + str(index)].loc[parent_1]

                i, j, k, l = torsion.indices
                right_mask = torsion.right_mask

                dihedrals.append(dihedral)
                ase_object.set_dihedral(a1=i,
                                        a2=j,
                                        a3=k,
                                        a4=l,
                                        angle=float(dihedral),
                                        mask=right_mask)

            # Updating the molecule
            update_from_ase(autotst_object)

            constrained_energy, relaxed_energy, ase_copy = get_energies(autotst_object)

            relaxed_torsions = []

            for torsion in torsions:

                i, j, k, l = torsion.indices

                angle = round(ase_copy.get_dihedral(i,j,k,l), -1)
                angle = int(30 * round(float(angle)/30))
                if angle < 0: angle += 360
                relaxed_torsions.append(angle)

            r.append([constrained_energy, relaxed_energy] + dihedrals + relaxed_torsions)


        results = pd.DataFrame(r)
        logging.info("Creating the DataFrame of results for the {}th generation".format(gen_number))

        columns = ["constrained_energy", "relaxed_energy"]
        for i in range(len(torsions)):
            columns = columns + ["torsion_" + str(i)]

        for i in range(len(torsions)):
            columns = columns + ["relaxed_torsion_" + str(i)]
        results.columns = columns
        results = results.sort_values("constrained_energy")

        unique_conformers = get_unique_conformers(results, unique_conformers)

        if store_generations == True:
            # This portion stores each generation if desired
            logging.info("Saving the results DataFrame")

            generation_name = "{0}_ga_generation_{1}.csv".format(label, gen_number)
            f = os.path.join(store_directory, generation_name)
            results.to_csv(f)

        top = select_top_population(results, top_percent)

        stats = top.describe()

        if gen_number >= max_generations:
            complete = True
            logging.info("Max generations reached. GA complete.")
        if abs((stats["constrained_energy"]["max"] - stats["constrained_energy"]["min"]) / stats["constrained_energy"]["min"]) < tolerance:
            complete = True
            logging.info("Cutoff criteria reached. GA complete.")

    return results, unique_conformers
