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

import os, sys, logging, itertools
import pandas as pd
import numpy as np
import autotst
from ..reaction import Reaction, TS
from ..species import Species
from ..geometry import Bond, Angle, Torsion, CisTrans
import ase
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def get_energy(conformer):
    """
    A function to find the potential energy of a conformer

    Variables:
    - conformer (Conformer): the conformer object of interest

    Returns:
    - energy (float): the corresponding energy of that conformer. 
        Will result in an error if there is no ASECalculator object attached to the ase_molecule
    """

    energy = conformer.ase_molecule.get_potential_energy()

    return energy


def find_terminal_torsions(conformer):
    """
    A function that can find the terminal and non terminal torsions in a conformer object

    Variables:
    - conformer (Conformer): the conformer object of interest
    ### DEFINITION OF A TERMINAL TORSION ###
    # A terminal torsion is a symmetric top at the end of a molecule #

    Returns:
    - terminal_torsions (list): a list of the terminal torsions that are not important in a conformer search
    - non_terminal_torsions (list): a list of the internal torsions that are needed in a conformer search
    """
    terminal_torsions = []
    non_terminal_torsions = []
    for torsion in conformer.torsions:

        _, j, k, _ = torsion.atom_indices
        rmg_mol = conformer.rmg_molecule

        assert rmg_mol

        atom_j = rmg_mol.atoms[j]
        atom_k = rmg_mol.atoms[k]

        terminal = False

        if (atom_j.is_carbon()) and (len(atom_j.bonds) == 4):
            num_hydrogens = 0
            for atom_other in list(atom_j.bonds.keys()):
                if atom_other.is_hydrogen():
                    num_hydrogens += 1

            if num_hydrogens == 3:
                terminal = True

        if (atom_k.is_carbon()) and (len(atom_k.bonds) == 4):
            num_hydrogens = 0
            for atom_other in list(atom_k.bonds.keys()):
                if atom_other.is_hydrogen():
                    num_hydrogens += 1

            if num_hydrogens == 3:
                terminal = True

        if terminal:
            terminal_torsions.append(torsion)
        else:
            non_terminal_torsions.append(torsion)

    return terminal_torsions, non_terminal_torsions
