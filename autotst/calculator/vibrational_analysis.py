#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2020 Richard H. West (r.west@northeastern.edu)
#   and the AutoTST Team
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

import os
import logging
import pandas as pd
import numpy as np
import ase
import cclib.io 
from ..reaction import Reaction, TS
from ..species import Species, Conformer


def percent_change(original, new):
    "A function to calculate the percent change between two values"
    percent_change = (abs(new - original) / original) * 100
    return percent_change


class VibrationalAnalysis():
    """
    A class that allows one to perform vibrational analysis. It takes a
    TS and uses it to parse out the finalized geometry from a
    corresponding log file, and then compares the geometry before and after
    displacement from the imaginary frequency.
    """

    def __init__(self, transitionstate=None, directory="."):
        """
        Variables:
        - ts (TS): A TS object that you want to run vibratinal analysis on
        - scratch (str): the directory where log files of the TS are located
        """
        self.ts = transitionstate
        self.ts.get_geometries()
        self.directory = directory

        self.log_file = None

    def __repr__(self):
        if self.ts is None:
            label = None
        else:
            label = self.ts.reaction_label
        return f'<Vibrational Analysis "{label}">'

    def get_log_file(self):
        """
        This method obtains the logfile name from the TS

        Variables:
        - scratch (str): The directory where the log files are located in
        - ts (TS): the TS object of interest

        Returns:
        - log_file (str): a path to the log file of interest
        """
        self.log_file = os.path.join(
            self.directory,
            "ts",
            self.ts.reaction_label,
            "conformers",
            f"{self.ts.reaction_label}_{self.ts.direction}_{self.ts.index}.log")
        return self.log_file

    def parse_vibrations(self):
        """
        This method obtains the vibrations from the log file of interest using
        cclib. It then creates a zipped list with the vibrational frequencies
        and their corresponding displacements.

        Variables:
        - log_file (str): the log file you want to obtain vibrations from. Often found from `get_log_file`

        Returns:
        - vibrations (list): a list of the vibrations and their corresponding displacements in XYZ coordinates
        """

        if not self.log_file:
            self.get_log_file()

        assert os.path.exists(self.log_file), "Log file provided does not exist"

        log_file_info = cclib.io.ccread(self.log_file, loglevel=logging.ERROR)
        self.vibrations = list(zip(log_file_info.vibfreqs, log_file_info.vibdisps))
        
        return self.vibrations

    def obtain_geometries(self):
        """
        This method obtains the previbrational geometry (the geometry returned
        by a quantum optimizer), and the postvibrational geometry.

        Variables:
        - ts (TS): A transition state object of interest
        - vibrations (list): a list of the vibrations and their corresponding displacements in XYZ coords. Often from `parse_vibrations`

        Returns:
        - pre_geometry (ase.Atoms): an ase.Atoms object containing the geometry before applying the vibrations
        - post_geometry (ase.Atoms): an ase.Atoms object containing the geometry after applying the vibrations
        """

        assert isinstance(self.ts, TS)

        
        symbol_dict = {
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }
        atoms = []

        parser = cclib.io.ccread(self.log_file, loglevel=logging.ERROR)

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
            atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))
        
        self.ts._ase_molecule = ase.Atoms(atoms)
        self.ts.update_coords_from("ase")

        self.pre_geometry = self.ts.ase_molecule.copy()
        self.post_geometry = self.ts.ase_molecule.copy()

        for vib, displacements in self.vibrations:
            if vib < 0:  # Finding the imaginary frequency
                self.post_geometry.arrays["positions"] -= displacements

        return self.pre_geometry, self.post_geometry

    def obtain_percent_changes(self):
        """
        This method takes the connectivity of a TS and then uses
        that to identify the percent change of the bonds
        between the before_geometry and the post_geometryself.

        This returns a dataframe with the type of geometry, the indicies in it,
        if it is close to the reaction center, and the percent change of that
        geometry.

        Variables:
        - ts (TS): the transition state of interest
        - berfore_geometry (ase.Atoms): the ase.Atoms object describing the geometry before applying the vibrations
        - post_geometry (ase.Atoms): the ase.Atoms object describing the geometry after applying the vibrations

        Returns:
        - results (DataFrame): a pandas dataframe containing the bond index, atom indicies of those bonds,
            if the bond is in the reaction center, and percent change of all bonds

        """

        results = []
        
        for bond in self.ts.bonds:
            i, j = bond.atom_indices
            before = self.pre_geometry.get_distance(i, j)
            after = self.post_geometry.get_distance(i, j)
            results.append([bond.index, bond.atom_indices,
                            bond.reaction_center, percent_change(before, after)])

        self.percent_changes = pd.DataFrame(results, columns=["index", "atom_indices", "center", "percent_change"])

        return self.percent_changes

    def validate_ts(self):
        """
        A method designed to run the above and return a bool that states if we
        have arrived at a TS. We say we have arrived at a TS if the average
        change of geometries in the reaction center is one order of magnitude
        geater elsewhere.

        Variables:
        - scratch (str): the path to where the log file of the optimized transition state is
        - ts (TS): the transition state of interest

        Returns:
        - (bool): True if the TS can be validated via vibrational analysis and False if it cannot
        """
        try:
            self.get_log_file()

            self.parse_vibrations()

            self.obtain_geometries()

            self.percent_changes = self.obtain_percent_changes()


            center_values = np.log(
                self.percent_changes[self.percent_changes.center].percent_change.mean())
            shell_values = np.log(
                self.percent_changes[self.percent_changes.center != True].percent_change.mean())

            if center_values > shell_values + 1:
                logging.info("Vibrational analysis was successful")
                return True
            else:
                logging.info(
                    "Cannot reasonably say that we have arrived at a TS through vibrational analysis.")
                return False
        except AssertionError:
            logging.info("Something went wrong when attempting vibrational analysis...")
            logging.info("Cannot verify via vibrational analysis")
            return False
