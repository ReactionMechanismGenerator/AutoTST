#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
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
##########################################################################
import itertools
import logging
import pandas as pd
import numpy as np
import os

import ase
from ase import Atoms
from ase import calculators
from ase.optimize import BFGS

import autotst
from autotst.conformer.utilities import get_energy, find_terminal_torsions


def find_all_combos(
        conformer,
        delta=float(30),
        cistrans=True,
        chiral_centers=True):
    """
    A function to find all possible conformer combinations for a given conformer
    """

    terminal_torsions, torsions = find_terminal_torsions(conformer)
    cistranss = conformer.cistrans
    chiral_centers = conformer.chiral_centers

    torsion_angles = np.arange(0, 360, delta)
    torsion_combos = list(itertools.combinations_with_replacement(
        torsion_angles, len(torsions)))
    if len(torsions) != 1:
        torsion_combos = list(
            set(
                torsion_combos +
                list(itertools.combinations_with_replacement(
                    torsion_angles[::-1], len(torsions)
                ))))

    if cistrans:
        cistrans_options = ["E", "Z"]
        cistrans_combos = list(itertools.combinations_with_replacement(
            cistrans_options, len(cistranss)))
        if len(cistranss) != 1:
            cistrans_combos = list(
                set(
                    cistrans_combos +
                    list(itertools.combinations_with_replacement(
                        cistrans_options[::-1], len(cistranss)
                    ))))

    else:
        cistrans_combos = [()]

    if chiral_centers:
        chiral_options = ["R", "S"]
        chiral_combos = list(itertools.combinations_with_replacement(
            chiral_options, len(chiral_centers)))
        if len(chiral_centers) != 1:
            chiral_combos = list(
                set(
                    chiral_combos +
                    list(itertools.combinations_with_replacement(
                        chiral_options[::-1], len(chiral_centers)
                    ))))
    else:
        chiral_combos = [()]

    all_combos = list(
        itertools.product(
            torsion_combos,
            cistrans_combos,
            chiral_combos))
    return all_combos


def systematic_search(conformer,
                      delta=float(30),
                      cistrans=True,
                      chiral_centers=True,
                      store_results=False,
                      store_directory="."):
    """
    Perfoms a brute force conformer analysis of a molecule or a transition state

    :param autotst_object: am autotst_ts, autotst_rxn, or autotst_molecule that you want to perform conformer analysis on
       * the ase_object of the autotst_object must have a calculator attached to it.
    :param store_generations: do you want to store pickle files of each generation
    :param store_directory: the director where you want the pickle files stored
    :param delta: the degree change in dihedral angle between each possible dihedral angle

    :return results: a DataFrame containing the final generation
    :return unique_conformers: a dictionary with indicies of unique torsion combinations and entries of energy of those torsions
    """
    # Takes each of the molecule objects

    combos = find_all_combos(
        conformer,
        delta=delta,
        cistrans=cistrans,
        chiral_centers=chiral_centers)

    terminal_torsions, torsions = find_terminal_torsions(conformer)
    file_name = conformer.smiles + "_brute_force.csv"

    calc = conformer.ase_molecule.get_calculator()

    results = []
    for combo in combos:

        torsions, cistrans, chiral_centers = combo

        for i, torsion in enumerate(torsions):

            tor = conformer.torsions[i]
            i, j, k, l = tor.atom_indices
            mask = tor.mask

            conformer.ase_molecule.set_dihedral(
                a1=i,
                a2=j,
                a3=k,
                a4=l,
                angle=torsion,
                mask=mask
            )
            conformer.update_coords()

        for i, e_z in enumerate(cistrans):
            ct = conformer.cistrans[i]
            conformer.set_cistrans(ct.index, e_z)

        for i, s_r in enumerate(chiral_centers):
            center = conformer.chiral_centers[i]
            conformer.set_chiral_center(center.index, s_r)

        conformer.ase_molecule.set_calculator(calc)
        # Something funky happening with this optimization
        # conformer.update_coords()
        #opt = BFGS(conformer.ase_molecule)
        # opt.run()
        conformer.update_coords()
        try:
            energy = get_energy(conformer)
        except BaseException:
            energy = 1e5

        sample = ["torsion_{}", "cistrans_{}", "chiral_center_{}"]
        columns = ["energy"]
        long_combo = [energy]

        for c, name in zip(combo, sample):
            for i, info in enumerate(c):
                columns.append(name.format(i))
                long_combo.append(info)

        results.append(long_combo)

    brute_force = pd.DataFrame(results, columns=columns)

    if store_results:
        f = os.path.join(store_directory, file_name)
        brute_force.to_csv(f)

    from ase import units

    unique_conformers = []
    for i, ind in enumerate(brute_force[brute_force.energy < (
            brute_force.energy.min() + units.kcal / units.mol / units.eV)].index):
        copy_conf = conformer.copy()
        copy_conf.index = i
        for col in brute_force.columns[1:]:
            geometry = col.split("_")[0]
            index = int(col.split("_")[-1])
            value = brute_force.loc[ind][col]

            if "torsion" in geometry:
                copy_conf.set_torsion(index, float(value))
            elif "cistrans" in geometry.lower():
                copy_conf.set_cistrans(index, value)
            elif "chiral" in geometry.lower():
                copy_conf.set_chirality(index, value)

        unique_conformers.append(copy_conf)

    return brute_force, unique_conformers
