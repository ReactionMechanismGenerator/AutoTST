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
import itertools
import logging
import pandas as pd
import numpy as np
import os

import ase
from ase import Atoms
from ase import calculators

import autotst
from autotst.conformer.utilities import update_from_ase, get_unique_conformers, get_energies

def perform_brute_force(autotst_object,
                        delta=float(30),
                        store_results=True,
                        store_directory="."):

    # Takes each of the molecule objects
    if isinstance(autotst_object, autotst.molecule.AutoTST_Molecule):
        ase_object = autotst_object.ase_molecule
        torsions = autotst_object.torsions
        file_name = autotst_object.smiles + "_brute_force.csv"


    elif isinstance(autotst_object, autotst.reaction.AutoTST_Reaction):
        ase_object = autotst_object.ts.ase_ts
        torsions = autotst_object.ts.torsions
        file_name = autotst_object.label + "_brute_force.csv"

    elif isinstance(autotst_object, autotst.reaction.AutoTST_TS):
        ase_object = autotst_object.ase_ts
        torsions = autotst_object.torsions
        file_name = autotst_object.label + "_brute_force.csv"


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

        update_from_ase(autotst_object)

        constrained_energy, relaxed_energy, ase_copy = get_energies(autotst_object)

        relaxed_torsions = []

        for torsion in torsions:

            i, j, k, l = torsion.indices

            angle = round(ase_copy.get_dihedral(i,j,k,l), -1)
            angle = int(30 * round(float(angle)/30))
            if angle < 0: angle += 360
            relaxed_torsions.append(angle)

        results.append([constrained_energy, relaxed_energy] + list(combo) + relaxed_torsions)

    brute_force = pd.DataFrame(results)
    columns = ["constrained_energy", "relaxed_energy"]
    for i in range(len(torsions)):
        columns = columns + ["torsion_" + str(i)]

    for i in range(len(torsions)):
        columns = columns + ["relaxed_torsion_" + str(i)]

    brute_force.columns = columns

    if store_results:
        f = os.path.join(store_directory, file_name)
        brute_force.to_csv(f)

    return brute_force
