

class Torsion:
    """
    A class that acts as a container for torsion information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
    """

    def __init__(self, indices, dihedral, left_mask, right_mask):
        self.indices = indices
        self.dihedral = dihedral
        self.left_mask = left_mask
        self.right_mask = right_mask


class Angle:
    """
    A class that acts as a container for torsion information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
    """

    def __init__(self, indices, degree, left_mask, right_mask):
        self.indices = indices
        self.degree = degree
        self.left_mask = left_mask
        self.right_mask = right_mask