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

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction


from multi_molecule import *
from multi_reaction import *
from utilities import *

from ase.calculators.morse import * #chosing this calculator for now because it's fast
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *


def perform_ga(multi_object,
               df,
               top_percent=0.3,
               tolerance=1e-4,
               max_generations=np.inf,
               store_generations=False,
               store_directory=".",
               mutation_probability=0.2,
               delta=30):

    """
    Performs a genetic algorithm to determine the lowest energy conformer of a TS or molecule

    :param multi_object: a multi_ts, multi_rxn, or multi_molecule that you want to perform conformer analysis on
       * the ase_object of the multi_object must have a calculator attached to it.
    :param df: a DataFrame containing the initial population
    :param top_percent: float of the top percentage of conformers you want to select
    :param tolerance: float of one of the possible cut off points for the analysis
    :param max_generations: int of one of the possible cut off points for the analysis
    :param store_generations: do you want to store pickle files of each generation
    :param store_directory: the director where you want the pickle files stored
    :param mutation_probability: float of the chance of mutation
    :param delta: the degree change in dihedral angle between each possible dihedral angle
    :return df: a DataFrame containing the final generation
    """
    possible_dihedrals = np.arange(0, 360, delta)
    top = select_top_population(df,
                                top_percent=top_percent
                                )

    top_population = top.shape[0]

    population_size = df.shape[0]

    # Takes each of the molecule objects
    if "Multi_Molecule" in str(multi_object.__class__):
        ase_object = multi_object.ase_molecule
        torsions = multi_object.torsions

    elif "Multi_Reaction" in str(multi_object.__class__):
        ase_object = multi_object.multi_ts.ase_ts
        torsions = multi_object.multi_ts.torsions

    elif "Multi_TS" in str(multi_object.__class__):
        ase_object = multi_object.ase_ts
        torsions = multi_object.torsions

    gen_number = 0
    complete = False
    while complete == False:
        gen_number += 1
        logging.info("Performing GA on generation {}".format(gen_number))

        results = []
        for individual in range(population_size):
            parent_0, parent_1 = random.sample(np.arange(top.shape[0]), 2)
            dihedrals = []
            for index, torsion in enumerate(torsions):

                if random.random() < mutation_probability:
                    dihedral = np.random.choice(possible_dihedrals)
                    logging.info("Mutation on the {0}th torsion of the {1}th individual".format(index, individual))
                else:
                    logging.info("Crossover on the {0}th torsion of the {1}th individual".format(index, individual))
                    if 0.5 > random.random():
                        dihedral = df.iloc[parent_0, index + 1]
                    else:
                        dihedral = df.iloc[parent_1, index + 1]

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
                if "Multi_Molecule" in str(multi_object.__class__):
                    multi_object.update_geometry_from_ase_mol()

                elif "Multi_Reaction" in str(multi_object.__class__):
                    multi_object.multi_ts.update_ts_from_ase_ts()

                elif "Multi_TS" in str(multi_object.__class__):
                    multi_object.update_ts_from_ase_ts()

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
            logging.info("Pickling the DataFrame")

            if "Multi_Reaction" in str(multi_object.__class__):
                generation_name = "rxn_ga_generation_{}.pkl".format(gen_number)

            elif "Multi_Molecule" in str(multi_object.__class__):
                generation_name = "mol_ga_generation_{}.pkl".format(gen_number)

            elif "Multi_TS" in str(multi_object.__class__):
                generation_name = "ts_ga_generation_{}.pkl".format(gen_number)

            else:
                generation_name = "ga_generation_{}.pkl".format(gen_number)

            f = open(os.path.join(store_directory, generation_name), "w")
            pickle.dump(df, f)

        top = df.iloc[:int(top_population), :]

        if gen_number >= max_generations:
            complete = True
            logging.info("Max generations reached. GA complete.")
        if top.std()[0] / top.mean()[0] < tolerance:
            complete = True
            logging.info("Cutoff criteria reached. GA complete.")

    return df
