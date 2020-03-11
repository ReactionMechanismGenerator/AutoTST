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


class Bond:
    """
    A class that acts as a container for bond information.

    Variables:
    - index (int): the unique bond index
    - atom_indicies (list): indicies of atoms that describe the bond
    - length (float): the length of the bond in angstroms
    - reaction_center (bool): a bool that states if this bond is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    - mask (list): a list of bools that describe if atoms are on the right side of the bond
    """

    def __init__(self, index, atom_indices, length, mask=None, reaction_center=False):
        self.index = index
        self.atom_indices = atom_indices
        self.length = length
        self.reaction_center = reaction_center
        self.mask = mask

    def __repr__(self):
        return f'<Bond "{self.atom_indices}">'


class Angle():
    """
    A class that acts as a container for angle information.

    Variables:
    - index (int): the unique angle index
    - atom_indicies (list): indicies of atoms that describe the angle
    - degree (float): the angle of the bond in degrees
    - reaction_center (bool): a bool that states if this angle is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    - mask (list): a list of bools that describe if atoms are on the right side of the angle
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
        return f'<Angle "{self.atom_indices}">'


class Torsion():
    """
    A class that acts as a container for torsion information.

    Variables:
    - index (int): the unique torsion index
    - atom_indicies (list): indicies of atoms that describe the torsion
    - dihedral (float): the angle of the dihedral in degrees
    - reaction_center (bool): a bood that states if this torsion is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    - mask (list): a list of bools that describe if atoms are on the right side of the torsion
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
        return f'<Torsion "{self.atom_indices}">'


class CisTrans():
    """
    A class that acts as a container for CisTrans information.

    Variables:
    - index (int): the unique cistrans index
    - atom_indicies (list): indicies of atoms that describe the cistrans
    - dihedral (float): the angle of the dihedral in degrees
    - reaction_center (bool): a bood that states if this cistrans is in the reaction center
        i.e. atoms i and j are labeled (*1, *2, *3)
    - mask (list): a list of bools that describe if atoms are on the right side of the cistrans
    - stero (str): the sterochemistry of this cistrans bond
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
        return f'<CisTrans "{self.atom_indices} - {self.stero}">'


class ChiralCenter():
    """
    A class that acts as a container for ChiralCenter information.

    Variables:
    - index (int): the unique chiralcenter index
    - atom_index (int): index of atom that describe the chiralcenter
    - chirality (str): the sterochemistry of this chiralcenter bond
    """

    def __init__(self, index, atom_index, chirality):
        self.index = index
        self.atom_index = atom_index
        self.chirality = chirality

    def __repr__(self):
        return '<ChiralCenter "{0} - {1}">'.format(
            self.atom_index, self.chirality)
