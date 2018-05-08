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

    f = os.path.join(store_directory, file_name)
    df.to_csv(f)

    return brute_force
