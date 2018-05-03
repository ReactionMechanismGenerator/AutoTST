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

class Torsion:
    """
    A class that acts as a container for torsion information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
        i.e. (atoms i, j, and k) or (atoms j, k, and l) are labeled (*1, *2, *3)
    """

    def __init__(self, indices, dihedral, left_mask, right_mask, reaction_center="No"):
        self.indices = indices
        self.dihedral = dihedral
        self.left_mask = left_mask
        self.right_mask = right_mask
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<AutoTST Torsion "{0}">'.format(self.indicies)


class Angle:
    """
    A class that acts as a container for angle information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
    * reaction_center (bool): a bool that states if this angle is in the reaction center
        i.e. atoms i, j and, k are labeled (*1, *2, *3)
    """

    def __init__(self, indices, degree, left_mask, right_mask, reaction_center="No"):
        self.indices = indices
        self.degree = degree
        self.left_mask = left_mask
        self.right_mask = right_mask
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<AutoTST Angle "{0}">'.format(self.indicies)

class Bond:
    """
    A class that acts as a container for bond information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * length (float): the angle (in degrees) of the dihedral
    * reaction_center (bool): a bood that states if this bond is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    """

    def __init__(self, indices, length, reaction_center="No"):
        self.indices = indices
        self.length = length
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<AutoTST Bond "{0}">'.format(self.indicies)
