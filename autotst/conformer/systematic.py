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
import multiprocessing
from multiprocessing import Process, Manager
import time

import ase
from ase import Atoms
from ase import calculators
from ase.calculators.calculator import FileIOCalculator
from ase.optimize import BFGS

import autotst
from autotst.species import Conformer
from autotst.reaction import TS
from autotst.conformer.utilities import get_energy, find_terminal_torsions


def find_all_combos(
        conformer,
        delta=float(60),
        cistrans=True,
        chiral_centers=True):
    """
    A function to find all possible conformer combinations for a given conformer

    Params:
    - conformer (`Conformer`) an AutoTST `Conformer` object of interest
    - delta (int or float): a number between 0 and 180 or how many conformers to generate per dihedral
    - cistrans (bool): indication of if one wants to consider cistrans bonds
    - chiral_centers (bool): indication of if one wants to consider chiral centers bonds

    Returns:
    - all_combos (list): a list corresponding to the number of unique conformers created.
    """

    conformer.get_geometries()

    _, torsions = find_terminal_torsions(conformer)
    cistranss = conformer.cistrans
    chiral_centers = conformer.chiral_centers

    torsion_angles = np.arange(0, 360, delta)
    torsion_combos = list(itertools.product(
        torsion_angles, repeat=len(torsions)))

    if cistrans:
        cistrans_options = ["E", "Z"]
        cistrans_combos = list(itertools.product(
            cistrans_options, repeat=len(cistranss)))



    else:
        cistrans_combos = [()]

    if chiral_centers:
        chiral_options = ["R", "S"]
        chiral_combos = list(itertools.product(
            chiral_options, repeat=len(chiral_centers)))

    else:
        chiral_combos = [()]

    all_combos = list(
        itertools.product(
            torsion_combos,
            cistrans_combos,
            chiral_combos))
    return all_combos


def systematic_search(conformer,
                      delta=float(60),
                      cistrans=True,
                      chiral_centers=True,
                      tmp="."):
    """
    Perfoms a systematic conformer analysis of a `Conformer` or a `TS` object

    Variables:
    - conformer (`Conformer` or `TS`): a `Conformer` or `TS` object of interest
    - delta (int or float): a number between 0 and 180 or how many conformers to generate per dihedral
    - cistrans (bool): indication of if one wants to consider cistrans bonds
    - chiral_centers (bool): indication of if one wants to consider chiral centers bonds

    Returns:
    - confs (list): a list of unique `Conformer` objects within 1 kcal/mol of the lowest energy conformer determined
    """
    # Takes each of the molecule objects

    combos = find_all_combos(
        conformer,
        delta=delta,
        cistrans=cistrans,
        chiral_centers=chiral_centers)

    if len(combos) == 0:
        logging.info(
            "This species has no torsions, cistrans bonds, or chiral centers")
        logging.info("Returning origional conformer")
        return [conformer]

    _, torsions = find_terminal_torsions(conformer)

    calc = conformer.ase_molecule.get_calculator()
    if isinstance(calc, FileIOCalculator):
        logging.info("The calculator generates input and output files.")

    results = []
    conformers = {}
    combinations = {}
    for index, combo in enumerate(combos):

        combinations[index] = combo

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
            conformer.set_chirality(center.index, s_r)

        conformer.update_coords_from("ase")

        conformers[index] = conformer.copy()

    logging.info(
        "There are {} unique conformers generated".format(len(conformers)))


    def opt_conf(conformer, calculator, i):
        """
        A helper function to optimize the geometry of a conformer.
        Only for use within this parent function
        """

        labels = []
        for bond in conformer.bonds:
            labels.append(bond.atom_indices)

        if isinstance(conformer, TS):
            label = conformer.reaction_label
            ind1 = conformer.rmg_molecule.getLabeledAtom("*1").sortingLabel
            ind2 = conformer.rmg_molecule.getLabeledAtom("*3").sortingLabel
            labels.append([ind1, ind2])
            type = 'ts'
        else:
            label = conformer.smiles
            type = 'species'

        if isinstance(calc, FileIOCalculator):
            calculator.directory(tmp, label, '{}_{}'.format(conformer.smiles, i))
            if not os.path.exists(calculator.directory):
                try:
                    os.mkdirs(calculator.directory)
                except:
                    logging.info("An error occured when creating {}".format(calculator.directory))

            calculator.atoms = conformer.ase_molecule

        from ase.constraints import FixBondLengths
        c = FixBondLengths(labels)
        conformer.ase_molecule.set_constraint(c)

        conformer.ase_molecule.set_calculator(calculator)

        opt = BFGS(conformer.ase_molecule, logfile=None)
        opt.run()

        conformer.update_coords_from("ase")
        energy = get_energy(conformer)
        return_dict[i] = (energy, conformer.ase_molecule.arrays,
                          conformer.ase_molecule.get_all_distances())

    manager = Manager()
    return_dict = manager.dict()

    processes = []
    for i, conf in list(conformers.items()):
        p = Process(target=opt_conf, args=(conf, calc, i))
        processes.append(p)

    active_processes = []
    for process in processes:
        if len(active_processes) < multiprocessing.cpu_count():
            process.start()
            active_processes.append(process)
            continue

        else:
            one_done = False
            while not one_done:
                for i, p in enumerate(active_processes):
                    if not p.is_alive():
                        one_done = True
                        break

            process.start()
            active_processes[i] = process
    complete = np.zeros_like(active_processes, dtype=bool)
    while not np.all(complete):
        for i, p in enumerate(active_processes):
            if not p.is_alive():
                complete[i] = True

    from ase import units
    results = []
    for _, values in list(return_dict.items()):
        results.append(values)

    df = pd.DataFrame(results, columns=["energy", "arrays", 'distances'])
    df = df[df.energy < df.energy.min() + units.kcal / units.mol /
            units.eV].sort_values("energy")

    tolerance = 0.1
    scratch_index = []
    unique_index = []
    for index in df.index:
        if index in scratch_index:
            continue
        unique_index.append(index)
        scratch_index.append(index)
        distances = df.distances[index]
        for other_index in df.index:
            if other_index in scratch_index:
                continue

            other_distances = df.distances[other_index]

            if tolerance > np.sqrt((distances - other_distances)**2).mean():
                scratch_index.append(other_index)

    logging.info("We have identified {} unique conformers for {}".format(
        len(unique_index), conformer))
    confs = []
    i = 0
    for info in df[["energy", "arrays"]].loc[unique_index].values:
        copy_conf = conformer.copy()

        energy, array = info
        copy_conf.energy = energy
        copy_conf.ase_molecule.set_positions(array["positions"])
        copy_conf.update_coords_from("ase")
        c = copy_conf.copy()
        c.index = i
        confs.append(c)
        i += 1

    return confs
