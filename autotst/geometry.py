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


class Bond:
    """
    A class that acts as a container for bond information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * length (float): the angle (in degrees) of the dihedral
    * reaction_center (bool): a bood that states if this bond is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    """

    def __init__(self, index, atom_indices, length, reaction_center=False):
        self.index = index
        self.atom_indices = atom_indices
        self.length = length
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<Bond "{}">'.format(self.atom_indices)


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

    def __init__(
            self,
            index,
            atom_indices,
            degree,
            mask,
            reaction_center=False):
        self.index = index
        self.atom_indices = atom_indices
        self.degree = degree
        self.mask = mask
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<Angle "{0}">'.format(self.atom_indices)


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

    def __init__(
            self,
            index,
            atom_indices,
            dihedral,
            mask,
            reaction_center=False):
        self.index = index
        self.atom_indices = atom_indices
        self.dihedral = dihedral
        self.mask = mask
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<Torsion "{0}">'.format(self.atom_indices)


class CisTrans():
    """
    A class that acts as a container for CisTrans information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
        i.e. (atoms i, j, and k) or (atoms j, k, and l) are labeled (*1, *2, *3)
    """

    def __init__(
            self,
            index,
            atom_indices,
            dihedral,
            mask,
            stero=None,
            reaction_center=False):
        self.index = index
        self.atom_indices = atom_indices
        self.dihedral = dihedral
        self.mask = mask
        self.stero = stero
        self.reaction_center = reaction_center

    def __repr__(self):
        return '<CisTrans "{0} - {1}">'.format(self.atom_indices, self.stero)


class ChiralCenter():

    def __init__(self, index, atom_index, chirality):
        self.index = index
        self.atom_indices = atom_index
        self.chirality = chirality

    def __repr__(self):
        return '<ChiralCenter "{0} - {1}">'.format(
            self.atom_index, self.chirality)
