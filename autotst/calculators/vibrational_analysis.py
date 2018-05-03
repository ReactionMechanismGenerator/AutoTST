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

from autotst.reaction import *
from cclib.io import ccread

def percent_change(original,new):
    "A function to calculate the percent change between two values"
    percent_change = (abs(new - original) / original) * 100
    return percent_change


class Vibrational_Analysis():
    """
    A class that allows one to perform vibrational analysis. It takes an
    AutoTST_Reaction and uses it to parse out the finalized geometry from a
    corresponding log file, and then compares the geometry before and after
    displacement from the imaginary frequency.
    """

    def __init__(self, reaction):
        """
        reaction: (AutoTST_Reaction) a reaction that proives the connectivity
        and label name for this analysis
        """
        self.reaction = reaction

    def __repr__(self):
        return '<AutoTST Vibrational Analysis "{0}">'.format(self.reaction.label)

    def get_log_file(self, reaction):
        """
        This method obtains the logfile name from the AutoTST_Reaction
        """

        self.log_file = reaction.label + "_overall.log"

    def parse_vibrations(self):
        """
        This method obtains the vibrations from the log file of interest using
        cclib. It then creates a zipped list with the vibrational frequencies
        and their corresponding displacements.
        """

        assert os.path.exists(self.log_file)

        log_file_info = ccread(self.log_file)
        self.vibrations = zip(log_file_info.vibfreqs, log_file_info.vibdisps)

        return self.vibrations

    def obtain_geometries(self, reaction):
        """
        This method obtains the previbrational geometry (the geometry returned
        by a quantum optimizer), and the postvibrational geometry.
        """

        assert isinstance(reaction, AutoTST_Reaction)

        self.before_geometry = reaction.ts.ase_ts.copy()
        self.post_geometry = reaction.ts.ase_ts.copy()

        for vib, displacements in reaction.vibrations:
            if vib < 0: # Finding the imaginary frequency
                got_imaginary_frequency = True
                self.post_geometry.arrays["positions"] -= displacements

        return self.before_geometry, self.post_geometry

    def obtain_percent_changes(self, reaction):
        """
        This method takes the connectivity of an AutoTST_Reaction and then uses
        that to identify the percent change of the bonds, angles, and dihedrals
        between the before_geometry and the post_geometryself.

        This returns a dataframe with the type of geometry, the indicies in it,
        if it is close to the reaction center, and the percent change of that
        geometry.
        """

        results = []
        for torsion in reaction.ts.torsions:
            i,j,k,l = torsion.indices
            before = self.before_geometry.get_dihedral(i,j,k,l)
            after = self.post_geometry.get_dihedral(i,j,k,l)
            results.append(["Tor", torsion.indices, torsion.reaction_center, percent_change(before, after)])

        for angle in reaction.ts.angles:
            i,j,k = angle.indices
            before = self.before_geometry.get_angle(i,j,k)
            after = self.post_geometry.get_angle(i,j,k)
            results.append(["Ang", angle.indices, angle.reaction_center, percent_change(before, after)])

        for bond in reaction.ts.bonds:
            i,j = bond.indices
            before = self.before_geometry.get_distance(i,j)
            after = self.post_geometry.get_distance(i,j)
            results.append(["Bond", bond.indices, bond.reaction_center, percent_change(before, after)])

        results = pd.DataFrame(results)
        results.columns = ["type", "index", "center", "percent_change"]

        self.percent_changes = results

    def validate_ts(self):
        """
        A method designed to run the above and return a bool that states if we
        have arrived at a TS. We say we have arrived at a TS if the average
        change of geometries in the reaction center is one order of magnitude
        geater elsewhere.
        """

        self.get_log_file(self.reaction)

        self.parse_vibrations()

        self.obtain_geometries(self.reaction)

        self.obtain_percent_changes(self.reaction)

        if (np.log10(((self.percent_changes[self.percent_changes.center == "Yes"].mean()))) > np.log10(((self.percent_changes[self.percent_changes.center != "Yes"].mean()))) + 1).all():
            logging.info("Vibrational analysis was successful")
            return True

        else:
            logging.info("Cannot reasonably say that we have arrived at a TS through vibrational analysis.\nRunning an IRC calc.")
            return False
