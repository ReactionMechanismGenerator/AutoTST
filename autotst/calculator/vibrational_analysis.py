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

import os
import logging
import pandas as pd
import numpy as np
from ase import Atom, Atoms
from cclib.io import ccread
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer


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

    def __init__(self, transitionstate=None, log_file=None, directory="."):
        """
        Variables:
        - ts (TS): A TS object that you want to run vibratinal analysis on
        - scratch (str): the directory where log files of the TS are located
        """
        self.ts = transitionstate
        self.ts.get_geometries()
        self.directory = directory

        if log_file:
            self.log_file = log_file
            if not os.path.exists(log_file):
                logging.warning('log_file path does not exist')
        else:
            self.log_file = None

    def __repr__(self):
        if self.ts is None:
            label = None
        else:
            label = self.ts.reaction_label
        return '<Vibrational Analysis "{0}">'.format(
            label)

    def get_log_file(self):
        """
        This method obtains the logfile name from the TS

        Variables:
        - scratch (str): The directory where the log files are located in
        - ts (TS): the TS object of interest

        Returns:
        - log_file (str): a path to the log file of interest
        """

        log_path = os.path.join(
            self.directory,
            "ts",
            self.ts.reaction_label,
            "conformers",
            "{}_{}_{}.log".format(self.ts.reaction_label,
                                  self.ts.direction, self.ts.index))

        assert os.path.exists(log_path),"Could not retrieve log file from path {}".format(path)

        self.log_file = log_path
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

        log_file_info = ccread(self.log_file, loglevel=logging.ERROR)

        try:
            self.vibrations = list(zip(log_file_info.vibfreqs, log_file_info.vibdisps))
        except:
            self.vibrations = None
            logging.info("Could not parse vibrations from {}".format(self.log_file))
        
        return self.vibrations

    def obtain_geometries(self):
        """
        This method obtains the previbrational geometry (the geometry returned
        by a quantum optimizer), and the postvibrational geometry.

        Variables:
        - ts (TS): A transition state object of interest
        - vibrations (list): a list of the vibrations and their corresponding displacements in XYZ coords. Often from `parse_vibrations`

        Returns:
        - pre_geometry (ASEAtoms): an ASEAtoms object containing the geometry before applying the vibrations
        - post_geometry (ASEAtoms): an ASEAtoms object containing the geometry after applying the vibrations
        """

        assert isinstance(self.ts, TS)
        assert self.vibrations is not None,"Failed to get vibrations from {}".format(self.log_file)

        
        symbol_dict = {
            35: "Br",
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }
        atoms = []

        parser = ccread(self.log_file, loglevel=logging.ERROR)

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
            atoms.append(Atom(symbol=symbol_dict[atom_num], position=coords))
        
        self.ts.ase_molecule = Atoms(atoms)
        self.ts.update_coords_from("ase")
        self.ts.get_geometries()

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
        - berfore_geometry (ASEAtoms): the ASEAtoms object describing the geometry before applying the vibrations
        - post_geometry (ASEAtoms): the ASEAtoms object describing the geometry after applying the vibrations

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
            if self.log_file is None:
                self.get_log_file()

            self.parse_vibrations()

            self.obtain_geometries()

            self.percent_changes = self.obtain_percent_changes()


            center_values = np.log(
                self.percent_changes[self.percent_changes.center].percent_change.mean())
            shell_values = np.log(
                self.percent_changes[self.percent_changes.center != True].percent_change.mean())
            logging.info("center_values = {}".format(center_values))
            logging.info("shell_values = {}".format(shell_values))

            if center_values > shell_values + 1:
                logging.info("center_values ({0}) is greater than shell_values + 1 ({1})".format(center_values,shell_values+1))
                logging.info("Vibrational analysis was successful")
                return True
            else:
                logging.info("center_values ({0}) is not greater than shell_values + 1 ({1})".format(center_values, shell_values+1))
                logging.info(
                    "Cannot reasonably say that we have arrived at a TS through vibrational analysis.")
                return False
        except:
            logging.info("Something went wrong when attempting vibrational analysis...")
            logging.info("Cannot verify via vibrational analysis")
            return False
