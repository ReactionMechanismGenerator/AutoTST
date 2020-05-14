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
import rmgpy
import rmgpy.molecule


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

    def __init__(self, transitionstate=None, directory=".", log_file=None):
        """
        Variables:
        - ts (TS): A TS object that you want to run vibratinal analysis on
        - directory (str): the directory where log files of the TS are located
        - log_file (str): path to the TS log file
        """
        self.ts = transitionstate
        self.ts.get_geometries()
        self.directory = directory

        self.log_file = log_file
        
        if self.log_file is None:
            self.get_log_file()
        
        try:
            self.parser = cclib.io.ccread(self.log_file, loglevel=logging.ERROR)
        except:
            self.parser = None

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
            self.parser = cclib.io.ccread(
                self.log_file, loglevel=logging.ERROR)

        self.vibrations = list(zip(self.parser.vibfreqs, self.parser.vibdisps))
        
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
            35: "Br",
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }
        atoms = []

        for atom_num, coords in zip(self.parser.atomnos, self.parser.atomcoords[-1]):
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
            if self.log_file is None:
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

    def _get_complexes_from_rxn_label(self):
        """
        Returns a tuple of rmgpy.molecule.Molecule objects of
        reactant and product complexes from the TS reaction label.
        Returns:
            (Reactants Complex (rmgpy.molecule.Molecule),
            Products Complex (rmgpy.molecule.Molecule))
        """

        complexes = []
        for i in (0, 1):
            smiles = self.ts.reaction_label.split('_')[i].split('+')
            spcs = [rmgpy.molecule.Molecule(
                smiles=smi).to_single_bonds() for smi in smiles]
            spcs_complex = spcs.pop(0)
            for s in spcs:
                spcs_complex = spcs_complex.merge(s)
            spcs_complex.update_connectivity_values()
            spcs_complex.update_lone_pairs()
            spcs_complex.update_multiplicity()
            complexes.append(spcs_complex.copy(deep=True))
        
        return complexes


    def validate_by_connecting_the_dots(self):
        """
        Validates TS by:
        1) Creating 2 'freq' rmgpy.molecule complexes from negative frequency atom displacements
        2) Checking isomophism of 'freq' complexes to reactant and product complexes from reaction label
        Returns True if Validated and False if Invalid TS
        """

        try:
            self.parser.vibfreqs
        except AttributeError:
            logging.warning(f"{self.ts} has no vibrational frequencies. Vibrational analysis FAILED")
            return False

        freq = self.parser.vibfreqs[0]
        if freq > 0:
            logging.warning(
                f"{self.ts} does not have a negative frequency. Vibrational analysis FAILED")
            return False

        rxn_label_complexes = self._get_complexes_from_rxn_label()

        def make_complexes(coeff):

            complex1 = rmgpy.molecule.Molecule()
            complex1.from_xyz(
                self.parser.atomnos,
                self.parser.atomcoords[-1] + coeff*self.parser.vibdisps[0],
                raise_atomtype_exception=False,
            )
            complex2 = rmgpy.molecule.Molecule()
            complex2.from_xyz(
                self.parser.atomnos,
                self.parser.atomcoords[-1] - coeff*self.parser.vibdisps[0],
                raise_atomtype_exception=False
            )

            for c in (complex1, complex2):
                c.update_connectivity_values()
                c.update_lone_pairs()
                c.update_multiplicity()

            return (complex1, complex2)

        for i in np.linspace(0.25, 2, 8):

            complex1, complex2 = make_complexes(i)

            if complex1.is_isomorphic(rxn_label_complexes[0]) \
            and complex2.is_isomorphic(rxn_label_complexes[1]):
                logging.info(f"Validated {self.ts} through vibrational analysis!")
                return True

            if complex1.is_isomorphic(rxn_label_complexes[1]) \
            and complex2.is_isomorphic(rxn_label_complexes[0]):
                logging.info(
                    f"Validated {self.ts} through vibrational analysis!")
                return True

        logging.warning(f"{self.ts} FAILED vibrational analysis with freq {freq}")
        return False
