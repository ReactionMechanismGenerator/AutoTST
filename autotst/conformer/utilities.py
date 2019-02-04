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
import ase
from ase import units
from ase.constraints import FixBondLengths
from ase.optimize import BFGS

import autotst
from autotst.geometry import Bond, Angle, Torsion, CisTrans
from autotst.molecule import Molecule
from autotst.reaction import Reaction, TS


def update_from_ase(autotst_obj):
    """
    A function designed to update all objects based off of their ase objects
    """
    if isinstance(autotst_obj, autotst.molecule.Molecule):
        autotst_obj.update_from_ase_mol()

    if isinstance(autotst_obj, autotst.reaction.Reaction):
        autotst_obj.ts.update_from_ase_ts()

    if isinstance(autotst_obj, autotst.reaction.TS):
        autotst_obj.update_from_ase_ts()



def create_initial_population(autotst_object, delta=30, population_size=30):
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
    logging.info("Creating initial population of {} individuals from random guesses".format(
        population_size))

    df = None

    possible_dihedrals = np.arange(0, 360, delta)
    population = []
    if isinstance(autotst_object, autotst.molecule.Molecule):
        logging.info("The object given is a `Molecule` object")

        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_molecule

    if isinstance(autotst_object, autotst.reaction.Reaction):
        logging.info("The object given is a `Reaction` object")
        torsions = autotst_object.ts.torsions
        ase_object = autotst_object.ts.ase_ts

    if isinstance(autotst_object, autotst.reaction.TS):
        logging.info("The object given is a `TS` object")
        torsions = autotst_object.torsions
        ase_object = autotst_object.ase_ts
        
    terminal_torsions, non_terminal_torsions = find_terminal_torsions(autotst_object)

    for indivudual in range(population_size):
        dihedrals = []
        for torsion in non_terminal_torsions:
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

        energy = get_energy(autotst_object)

        population.append(
            [energy] + dihedrals)

    if len(population) > 0:
        logging.info("Creating a dataframe of the initial population")
        df = pd.DataFrame(population)
        columns = ["energy"]
        for i in range(len(non_terminal_torsions)):
            columns = columns + ["torsion_" + str(i)]

        df.columns = columns
        df = df.sort_values("energy")

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
    assert "energy" in df.columns
    df.sort_values("energy")
    population_size = df.shape[0]
    top_population = population_size * top_percent
    top = df.iloc[:int(top_population), :]

    return top


def get_unique_conformers(df, unique_torsions={}, min_rms=60):
    """
    A function designed to identify all low energy conformers within a standard deviation of the data given.

    :param:
     df: a DataFrame of a population of torsions with columns of `Energy` and `Torsion N`
     unique_torsions: a dict of unique torsions already present in the dataframe

    :return:
     unique_torsion: a dict containing all of the unique torsions from the dataframe appended
    """
    ### CURRENTLY BROKEN
    columns = []

    for c in df.columns:
        if "torsion" in c:
            columns.append(c)

    assert len(columns) > 0
    assert "energy" in df.columns

    mini = df[df.energy < (df.energy.min(
    ) + (units.kcal / units.mol) / units.eV)].sort_values("energy")
    for i, combo in enumerate(mini[columns].values):
        c = []
        for tor in combo:
            c.append(tor % 360)
        combo = np.array(c)
        energy = mini.energy.iloc[i]
        
        unique = []
        for key in unique_torsions.keys():
            key = np.array(key)
            
            if ((key - combo)**2).mean() > min_rms:
                unique.append(True)
                
            else:
                unique.append(False)
                
        if np.array(unique).all():
            combo = tuple(combo)
            unique_torsions[combo] = energy
            
            
    min_e = min(unique_torsions.values())

    to_delete = []
    for key, value  in unique_torsions.iteritems():

        if value > min_e + ((units.kcal / units.mol) / units.eV):
            to_delete.append(key)

    for delete in to_delete:
        del unique_torsions[delete]
                
        
    return unique_torsions

def partial_optimize_mol(autotst_object):
    """
    # CURRENTLY BROKEN
    """
    
    if isinstance(autotst_object, autotst.molecule.Molecule):
        ase_object = autotst_object.ase_molecule
        labels = []

    if isinstance(autotst_object, autotst.reaction.Reaction):
        ase_object = autotst_object.ts.ase_ts

        labels = []
        for atom in autotst_object.ts.rmg_ts.getLabeledAtoms().values():
            labels.append(atom.sortingLabel)

    if isinstance(autotst_object, autotst.reaction.TS):
        ase_object = autotst_object.ase_ts

        labels = []
        for atom in autotst_object.ts.rmg_ts.getLabeledAtoms().values():
            labels.append(atom.sortingLabel)
            
    ase_copy = ase_object.copy()
    ase_copy.set_calculator(ase_object.get_calculator())

    ase_copy.set_constraint(FixBondLengths(
        list(itertools.combinations(labels, 2))))

    opt = BFGS(ase_copy)
    complete = opt.run(fmax=0.01, steps=5)
    
    relaxed_energy = ase_copy.get_potential_energy()
    return relaxed_energy, ase_copy

def get_energy(conformer):
    """
    A function to find the potential energy of a conformer
    """

    energy = conformer.ase_molecule.get_potential_energy()

    return energy


def find_terminal_torsions(conformer):
    """
    A function that can find the terminal and non terminal torsions in a conformer object
    """
    terminal_torsions = []
    non_terminal_torsions = []
    for torsion in conformer.torsions:

        i,j,k,l = torsion.atom_indices
        rmg_mol = conformer.rmg_molecule
            
        assert rmg_mol

        atom_j = rmg_mol.atoms[j]
        atom_k = rmg_mol.atoms[k]
        
        terminal = False

        if (atom_j.isCarbon()) and (len(atom_j.bonds) == 4):
            num_hydrogens = 0
            for atom_other in atom_j.bonds.keys():
                if atom_other.isHydrogen():
                    num_hydrogens += 1

            if num_hydrogens == 3:
                terminal = True

        if (atom_k.isCarbon()) and (len(atom_k.bonds) == 4):
            num_hydrogens = 0
            for atom_other in atom_k.bonds.keys():
                if atom_other.isHydrogen():
                    num_hydrogens += 1

            if num_hydrogens == 3:
                terminal = True
                      
        if terminal:
            terminal_torsions.append(torsion)
        else:
            non_terminal_torsions.append(torsion)
            
    return terminal_torsions, non_terminal_torsions
