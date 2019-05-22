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
import numpy as np
from copy import deepcopy

import rdkit
from rdkit import DistanceGeometry
from rdkit.Chem import rdDistGeom


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm3D import EmbedLib

import ase

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.species import Species as RMGSpecies
from rmgpy.reaction import Reaction as RMGReaction
from rmgpy.reaction import ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

import autotst
from autotst.data.base import DistanceData, TransitionStateDepository, TSGroups, TransitionStates
from autotst.species import Species, Conformer
from autotst.geometry import Torsion, Angle, Bond, CisTrans, ChiralCenter

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

try:
    import py3Dmol
except ImportError:
    logging.info("Error importing py3Dmol")


# Currently this is set up to only work with H_Abstraction
# TODO: Edit this so it works with other reaction families

class Reaction():

    rmg_database = None   # will be an RMGDatabase instance, once loaded.
    # a dictionary will have reaction family names as keys and
    # TransitionStates instances as values, once loaded.
    ts_databases = dict()

    def __init__(
            self,
            label=None,
            rmg_reaction=None,
            reaction_family="H_Abstraction"):

        self.possible_families = [  # These families (and only these) will be loaded from both RMG and AutoTST databases
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration"
        ]

        self.label = label
        self.rmg_reaction = rmg_reaction
        self.reaction_family = reaction_family
        self.rmg_database = None
        self.ts_databases = None
        self._ts = None
        self._distance_data = None

    def __repr__(self):
        return '<Reaction "{}">'.format(self.label)

    @property
    def ts(self):
        """ #TODO: DOESN'T WORK RN
        The TS transition state for this reaction.

        Calls create_ts_geometries() if it has not previously been found.
        To update, call create_ts_geometries() manually.

        Variables:
        - None

        Return:
        - _ts (dict): the dictionary containing the forward and reverse TSs
        """
        if self._ts is None:
            ts_dict = {}
            for direction, complex in self.get_rmg_complexes().items():
                ts = TS(
                    reaction_label=self.label,
                    #smiles=complex.toSMILES()
                    direction=direction,
                    rmg_molecule=complex,
                    reaction_family=self.reaction_family,
                    distance_data=self.distance_data
                )
                ts_dict[direction] = [ts]

            self._ts = ts_dict

        return self._ts

    @property
    def distance_data(self):
        """
        The distance data.

        Calls generate_distance_data() if it has not previously been found.
        To update, call create_distance_data()

        Variables:
        - None

        Return:
        - _distance_data (DistanceData): a container for the key distance data
        """
        if self._distance_data is None:
            self.generate_distance_data()
        return self._distance_data

    @classmethod
    def load_databases(self, force_reload=False):
        """
        Load the RMG and AutoTST databases, if they have not already been loaded,
        into the class level variables where they are stored.

        Variables:
        - force_reload (bool):if set to True then forces a reload, even if already loaded.

        Returns:
        - None
        """
        if self.rmg_database and self.ts_databases and not force_reload:
            return self.rmg_database, self.ts_databases

        rmg_database = RMGDatabase()
        database_path = rmgpy.settings['database.directory']

        logging.info("Loading RMG database from '{}'".format(database_path))


        self.possible_families = [  # These families (and only these) will be loaded from both RMG and AutoTST databases
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration"
        ]
        try:
            rmg_database.load(
                database_path,
                kineticsFamilies=self.possible_families,
                transportLibraries=[],
                reactionLibraries=[],
                seedMechanisms=[],
                thermoLibraries=[
                    'primaryThermoLibrary',
                    'thermo_DFT_CCSDTF12_BAC',
                    'CBS_QB3_1dHR'],
                solvation=False,
            )
        except IOError:
            logging.info(
                "RMG database not found at '{}'. This can occur if a git repository of the database is being"
                "used rather than the binary version".format(database_path))
            database_path = os.path.join(database_path, 'input')
            logging.info(
                "Loading RMG database instead from '{}'".format(database_path))
            rmg_database.load(
                database_path,
                kineticsFamilies=self.possible_families,
                transportLibraries=[],
                reactionLibraries=[],
                seedMechanisms=[],
                thermoLibraries=[
                    'primaryThermoLibrary',
                    'thermo_DFT_CCSDTF12_BAC',
                    'CBS_QB3_1dHR'],
                solvation=False,
            )

        self.rmg_database = rmg_database

        self.ts_databases = dict()
        for reaction_family in self.possible_families:
            ts_database = TransitionStates()
            path = os.path.join(
                autotst.settings['tst_database_path'],
                reaction_family)
            global_context = {'__builtins__': None}
            local_context = {'DistanceData': DistanceData}
            family = rmg_database.kinetics.families[reaction_family]
            ts_database.family = family
            ts_database.load(path, local_context, global_context)

            self.ts_databases[reaction_family] = ts_database

        return self.rmg_database, self.ts_databases

    def generate_distance_data(self):
        """
        Generates the distance estimates using group additivity.
        Requires self.rmg_reaction
        Stores it in self.distance_data.

        Variables:
        - rmg_reaction (RMGReaction): The RMGReaction of interest

        Returns:
        - None
        """
        if not (self.rmg_database or self.ts_databases):
            self.rmg_database, self.ts_databases = self.load_databases()
        self.get_labeled_reaction()
        assert self.rmg_reaction, "try calling get_rmg_reaction() first"
        self._distance_data = self.ts_databases[self.reaction_family].groups.estimateDistancesUsingGroupAdditivity(
            self.rmg_reaction)

        if np.isclose(self._distance_data.distances["d12"] + self._distance_data.distances["d23"],
                      self._distance_data.distances["d13"],
                      atol=0.01):
            logging.info(
                "Distance between *1 and *3 is too small, setting it to lower bound of uncertainty")

            self._distance_data.distances["d13"] -= self._distance_data.uncertainties["d13"] / 2

        logging.info("The distance data is as follows: \n{}".format(
            self._distance_data))

        return self._distance_data

    def generate_reactants_and_products(self):
        """
        A module that will generate AutoTST Species for a given reaction's 
        reactants and products

        Variabels:
        - rmg_reaction (RMGReaction): the RMGReaction of interest

        Returns:
        - reactants (list): a list of AutoTST Species corresponding to the reactnats
        - products (list): a list of AutoTST Species corresponding to the products
        """
        if self.rmg_reaction is None:
            self.get_rmg_reaction()

        self.reactants = []
        self.products = []
        for react in self.rmg_reaction.reactants:
            mol = Species(rmg_species=react)
            mol.generate_structures()
            self.reactants.append(mol)

        for prod in self.rmg_reaction.products:
            mol = Species(rmg_species=prod)
            mol.generate_structures()
            self.products.append(mol)

        return self.reactants, self.products

    def get_labeled_reaction(self):
        """
        A method that will return a labeled reaction given a reaction label or rmg_reaction
        A label or an rmg_reaction needs to be provided in order for this method to work.
        If both are provided, we assert that the label matches the reaction.

        Variables:
        - label (str): the reaction label of interest
        - rmg_reaction (RMGReaction): the reaction of interest

        Returns:
        - reaction (RMGReaction): An RMGReaction with labeled reactants and products
        - name (str): The string corresponding to the reaction family matched to the reaction of interest
        """
        if not (self.rmg_database or self.ts_databases):
            self.rmg_database, self.ts_databases = self.load_databases()

        assert (
            self.label or self.rmg_reaction), "You must provide a reaction or a reaction label"

        match = False
        if self.label:
            rmg_reactants = []
            rmg_products = []
            r, p = self.label.split("_")
            for react in r.split("+"):
                s = RMGMolecule(SMILES=react)
                rmg_reactants.append(s)
            for prod in p.split("+"):
                s = RMGMolecule(SMILES=prod)
                rmg_products.append(s)

            test_reaction = RMGReaction(
                reactants=rmg_reactants, 
                products=rmg_products)

            if self.rmg_reaction:
                assert self.rmg_reaction.isIsomorphic(
                    test_reaction), "The reaction label provided does not match the RMGReaction provided..."

            for name, family in list(self.rmg_database.kinetics.families.items()):
                if match:
                    break

                labeled_r, labeled_p = family.getLabeledReactantsAndProducts(
                    test_reaction.reactants, test_reaction.products)
                if not (labeled_r and labeled_p):
                    continue

                if ((len(labeled_r) > 0) and (len(labeled_p) > 0)):
                    logging.info("Matched reaction to {} family".format(name))

                    labeled_reactants = deepcopy(labeled_r)
                    labeled_products = deepcopy(labeled_p)
                    test_reaction.reactants = labeled_r[:]
                    test_reaction.products = labeled_p[:]
                    match = True
                    final_name = name
                    break

        elif self.rmg_reaction: #RMGReaction but no label
            rmg_reactants = []
            rmg_products = []
            for react in self.rmg_reaction.reactants:
                if isinstance(react, RMGSpecies):
                    rmg_reactants.append(react.molecule)
                elif isinstance(react, RMGMolecule):
                    rmg_reactants.append([react])

            for prod in self.rmg_reaction.products:
                if isinstance(prod, RMGSpecies):
                    rmg_products.append(prod.molecule)
                elif isinstance(prod, RMGMolecule):
                    rmg_products.append([prod])

            test_reactants = []
            test_products = []
            if len(rmg_reactants) == 1:
                test_reactants = rmg_reactants
            elif len(rmg_reactants) == 2:
                l1, l2 = rmg_reactants
                for i in l1:
                    for j in l2:
                        test_reactants.append([i, j])

            if len(rmg_products) == 1:
                test_reactants = rmg_products
            elif len(rmg_products) == 2:
                l1, l2 = rmg_products
                for i in l1:
                    for j in l2:
                        test_products.append([i, j])
            for name, family in list(self.rmg_database.kinetics.families.items()):
                if match:
                    break
                for test_reactant in test_reactants:
                    for test_products in test_products:

                        if match:
                            continue

                        test_reaction = RMGReaction(
                            reactants=test_reactant, products=test_products)

                        labeled_r, labeled_p = family.getLabeledReactantsAndProducts(
                            test_reaction.reactants, test_reaction.products)
                        if not (labeled_r and labeled_p):
                            continue

                        if ((len(labeled_r) > 0) and (len(labeled_p) > 0)):
                            logging.info(
                                "Matched reaction to {} family".format(name))

                            labeled_reactants = deepcopy(labeled_r)
                            labeled_products = deepcopy(labeled_p)
                            test_reaction.reactants = labeled_r[:]
                            test_reaction.products = labeled_p[:]
                            logging.info("\n{}".format(labeled_r))
                            logging.info("\n{}".format(labeled_p))
                            match = True
                            final_name = name
                            break
        assert match, "Could not idetify labeled reactants and products"

        reaction_list = self.rmg_database.kinetics.generate_reactions_from_families(
            test_reaction.reactants, test_reaction.products, only_families=[final_name])

        assert reaction_list, "Could not match a reaction to a reaction family..."

        for reaction in reaction_list:
            if test_reaction.isIsomorphic(reaction):
                reaction.reactants = labeled_reactants
                reaction.products = labeled_products
                break
        self.rmg_reaction = reaction
        self.reaction_family = name
        return self.rmg_reaction, self.reaction_family

    def get_label(self):
        """
        A method to get the reaction label corresponding to an rmg_reaction

        Variables:
        - rmg_reaction (RMGReaction): The RMGReaction of interest

        Returns:
        - string (str): the reaction label in the format r1+r2_p1+p2
        """
        if self.label:
            return self.label
        
        string = ""
        for react in self.rmg_reaction.reactants:
            if isinstance(react, RMGSpecies):
                string += "{}+".format(react.molecule[0].toSMILES())
            elif isinstance(react, RMGMolecule):
                string += "{}+".format(react.toSMILES())
        string = string[:-1]
        string += "_"
        for prod in self.rmg_reaction.products:
            if isinstance(prod, RMGSpecies):
                string += "{}+".format(prod.molecule[0].toSMILES())
            elif isinstance(prod, RMGMolecule):
                string += "{}+".format(prod.toSMILES())
        self.label = string[:-1]
        return self.label

    def get_rmg_reaction(self):

        if self.rmg_reaction:
            return self.rmg_reaction
        
        r, p = self.label.split("_")

        reactants = []
        products = []
        for react in r.split("+"):
            reactants.append(RMGMolecule(SMILES=react))
        for prod in p.split("+"):
            products.append(RMGMolecule(SMILES=prod))

        self.rmg_reaction = RMGReaction(reactants=reactants, products=products)
        return self.rmg_reaction

    def get_rmg_complexes(self):
        """
        A method to create a forward and reverse TS complexes used to initialize transition state geometries

        Variables:
        - rmg_reaction (RMGReaction): The RMGReaction of interest

        Returns:
        - complexes (dict): a dictionary containing RMGMolecules of the forward and reverse reaction complexes
        """

        if self.rmg_reaction is None:
            self.get_rmg_reaction()

        reactant_complex = RMGMolecule()
        for react in self.rmg_reaction.reactants:
            if isinstance(react, RMGMolecule):
                reactant_complex = reactant_complex.merge(react)
            elif isinstance(react, RMGSpecies):
                for mol in react.molecule:
                    if len(mol.getLabeledAtoms()) > 0:
                        reactant_complex = reactant_complex.merge(mol)

        product_complex = RMGMolecule()
        for prod in self.rmg_reaction.products:
            if isinstance(prod, RMGMolecule):
                product_complex = product_complex.merge(prod)
            elif isinstance(prod, RMGSpecies):
                for mol in prod.molecule:
                    if len(mol.getLabeledAtoms()) > 0:
                        product_complex = product_complex.merge(mol)

        reactant_complex.updateMultiplicity()
        product_complex.updateMultiplicity()

        if len(reactant_complex.getLabeledAtoms()) == 0 or len(product_complex.getLabeledAtoms()) == 0:
            logging.warning("REACTING ATOMS LABELES NOT PROVIDED. Please call `Reaction.get_labeled_reaction` to generate labeled complexes")

        self.complexes = {
            "forward": reactant_complex,
            "reverse": product_complex}

        return self.complexes

    def generate_conformers(self, ase_calculator=None):
        """
        A method to generate an ensemble of low energy conformers.
        Currently only supports a systematic search with the goal of adding evolutionary searches

        Variables: 
        - method (str): the method of the conformer search. Currently only systematic is supported
        - calculator (ASECalculator): the calculator you want to evaluate your conformers with.

        Returns:
        - ts (dict): a dictionary containing an ensemble of low energy transition state geometries in 
            the forward and reverse direction
        """

        assert ase_calculator, "Please provide an ASE calculator object"
        
        from autotst.conformer.systematic import systematic_search

        for direction, conformers in self.ts.items():
            conformer = conformers[0]
            conformer.ase_molecule.set_calculator(ase_calculator)
            conformers = systematic_search(conformer, delta=120)
            for conformer in conformers:
                conformer.direction = direction
            self.ts[direction] = conformers

        return self.ts


class TS(Conformer):
    """
    A class that defines the 3D geometry of a transition state (TS)
    """

    def __init__(
            self,
            smiles=None,
            reaction_label=None,
            direction='forward',
            rmg_molecule=None,
            reaction_family="H_Abstraction",
            distance_data=None,
            index=0):

        self.energy = None
        self.reaction_label = reaction_label
        self.direction = direction.lower()
        self.reaction_family = reaction_family
        self.distance_data = distance_data
        self.index = index
        self.bm = None

        assert direction in ["forward",
                             "reverse"], "Please provide a valid direction"

        self._rdkit_molecule = None
        self._ase_molecule = None

        if (smiles or rmg_molecule):
            if smiles and rmg_molecule:
                assert rmg_molecule.isIsomorphic(RMGMolecule(
                    SMILES=smiles)), "SMILES string did not match RMG Molecule object"
                self.smiles = smiles
                self.rmg_molecule = rmg_molecule

            elif rmg_molecule:
                self.rmg_molecule = rmg_molecule
                self.smiles = rmg_molecule.toSMILES()

            else:
                self.smiles = smiles
                self.rmg_molecule = RMGMolecule(SMILES=smiles)

            self.rmg_molecule.updateMultiplicity()
            self._symmetry_number = None

        else:
            self.smiles = None
            self.rmg_molecule = None
            self.rdkit_molecule = None
            self._pseudo_geometry = None
            self.ase_molecule = None
            self.bonds = []
            self.angles = []
            self.torsions = []
            self.cistrans = []
            self.chiral_centers = []
            self._symmetry_number = None

    def __repr__(self):
        return '<TS "{}">'.format(self.smiles)

    def copy(self):
        copy_conf = TS(
            reaction_label=self.reaction_label,
            reaction_family=self.reaction_family)
        copy_conf.smiles = self.smiles
        copy_conf.rmg_molecule = self.rmg_molecule.copy()
        copy_conf.rdkit_molecule = self.rdkit_molecule.__copy__()
        copy_conf._pseudo_geometry = self._pseudo_geometry.__copy__()
        copy_conf.ase_molecule = self.ase_molecule.copy()
        copy_conf.get_geometries()
        copy_conf.energy = self.energy
        copy_conf._symmetry_number = self._symmetry_number
        return copy_conf

    @property
    def symmetry_number(self):

        if not self._symmetry_number:
            self._symmetry_number = self.calculate_symmetry_number()
        return self._symmetry_number

    @property
    def rdkit_molecule(self):
        if (self._rdkit_molecule is None) and self.distance_data:
            self._rdkit_molecule = self.get_rdkit_mol()
        return self._rdkit_molecule

    @property
    def ase_molecule(self):
        if (self._ase_molecule is None):
            self._ase_molecule = self.get_ase_mol()
        return self._ase_molecule

    def get_rdkit_mol(self):
        """
        A method to create an rdkit geometry... slightly different than that of the conformer method
        returns both the rdkit_molecule and the bm
        """
        self.rdkit_molecule = Conformer(rmg_molecule=self.rmg_molecule).get_rdkit_mol()

        self.get_labels()
        for i, atom in enumerate(self.rmg_molecule.atoms):
            assert atom.number == self.rdkit_molecule.GetAtoms()[i].GetAtomicNum()

        if len(self.labels) == 3:

            rd_copy = Chem.RWMol(self.rdkit_molecule.__copy__())

            lbl1, lbl2, lbl3 = self.labels

            if not rd_copy.GetBondBetweenAtoms(lbl1, lbl2):
                rd_copy.AddBond(lbl1, lbl2,
                                order=rdkit.Chem.rdchem.BondType.SINGLE)
            else:
                rd_copy.AddBond(lbl2, lbl3,
                                order=rdkit.Chem.rdchem.BondType.SINGLE)

            self._pseudo_geometry = rd_copy

        logging.info("Initially embedded molecule")

        self.bm = None

        if self.distance_data:
            logging.info("Getting bounds matrix")
            self.bm = self.get_bounds_matrix()

            if len(self.labels) > 0:
                logging.info("Editing bounds matrix")
                self.bm = self.edit_matrix()

            logging.info("Performing triangle smoothing on bounds matrix.")
            DistanceGeometry.DoTriangleSmoothing(self.bm)

            logging.info("Now attempting to embed using edited bounds matrix.")

            self.rd_embed()
        return self.rdkit_molecule

    def get_bounds_matrix(self):
        """
        A method to obtain the bounds matrix
        """
        self.bm = rdDistGeom.GetMoleculeBoundsMatrix(self.rdkit_molecule)
        return self.bm

    def set_limits(self, lbl1, lbl2, value, uncertainty):
        """
        A method to set the limits of a particular distance between two atoms

        :param bm: an array of arrays corresponding to the bounds matrix
        :param lbl1: the label of one atom
        :param lbl2: the label of another atom
        :param value: the distance from a distance data object (float)
        :param uncertainty: the uncertainty of the `value` distance (float)
        :return bm: an array of arrays corresponding to the edited bounds matrix
        """
        logging.info(
            "For atoms {0} and {1} we have a distance of: \t {2}".format(
                lbl1, lbl2, value))
        if lbl1 > lbl2:
            self.bm[lbl2][lbl1] = value + uncertainty / 2
            self.bm[lbl1][lbl2] = max(0, value - uncertainty / 2)
        else:
            self.bm[lbl2][lbl1] = max(0, value - uncertainty / 2)
            self.bm[lbl1][lbl2] = value + uncertainty / 2

        return self.bm

    def bm_pre_edit(self, sect):
        """
        Clean up some of the atom distance limits before attempting triangle smoothing.
        This ensures any edits made do not lead to unsolvable scenarios for the molecular
        embedding algorithm.

        sect is the list of atom indices belonging to one species.
        """
        others = list(range(len(self.bm)))
        for idx in sect:
            others.remove(idx)

        for i in range(len(self.bm)):  # sect:
            for j in range(i):  # others:
                if i < j:
                    continue
                for k in range(len(self.bm)):
                    if k == i or k == j or i == j:
                        continue
                    Uik = self.bm[i, k] if k > i else self.bm[k, i]
                    Ukj = self.bm[j, k] if k > j else self.bm[k, j]

                    maxLij = Uik + Ukj - 0.1
                    if self.bm[i, j] > maxLij:
                        logging.info(
                            "Changing lower limit {0} to {1}".format(self.bm[i, j], maxLij))
                        self.bm[i, j] = maxLij

        return self.bm

    def get_labels(self):
        """
        A method to get the labeled atoms from a reaction

        :param reactants: a combined rmg_molecule object
        :return labels: the atom labels corresponding to the reaction center
        :return atomMatch: a tuple of tuples the atoms labels corresponding to the reaction center
        """

        if len(self.rmg_molecule.getLabeledAtoms()) == 0:
            labels = []
            atomMatch = ()

        if self.reaction_family.lower() in [
            'h_abstraction',
            'r_addition_multiplebond',
            'intra_h_migration']:
            # for i, atom in enumerate(reactants.atoms):
            lbl1 = self.rmg_molecule.getLabeledAtoms()["*1"].sortingLabel
            lbl2 = self.rmg_molecule.getLabeledAtoms()["*2"].sortingLabel
            lbl3 = self.rmg_molecule.getLabeledAtoms()["*3"].sortingLabel
            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,), (lbl2,), (lbl3,))
        elif self.reaction_family.lower() in ['disproportionation']:
            lbl1 = self.rmg_molecule.getLabeledAtoms()["*2"].sortingLabel
            lbl2 = self.rmg_molecule.getLabeledAtoms()["*4"].sortingLabel
            lbl3 = self.rmg_molecule.getLabeledAtoms()["*1"].sortingLabel

            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,), (lbl2,), (lbl3,))

        #logging.info("The labled atoms are {}.".format(labels))
        self.labels = labels
        self.atom_match = atomMatch
        return self.labels, self.atom_match

    def edit_matrix(self):
        """
        A method to edit the bounds matrix using labels and distance data
        """

        lbl1, lbl2, lbl3 = self.labels

        sect = []

        for atom in self.rmg_molecule.split()[0].atoms:
            sect.append(atom.sortingLabel)

        uncertainties = {'d12': 0.02, 'd13': 0.02, 'd23': 0.02}
        self.bm = self.set_limits(
            lbl1,
            lbl2,
            self.distance_data.distances['d12'],
            uncertainties['d12'])
        self.bm = self.set_limits(
            lbl2,
            lbl3,
            self.distance_data.distances['d23'],
            uncertainties['d23'])
        self.bm = self.set_limits(
            lbl1,
            lbl3,
            self.distance_data.distances['d13'],
            uncertainties['d13'])

        self.bm = self.bm_pre_edit(sect)

        return self.bm

    def optimize_rdkit_molecule(self):
        """
        Optimizes the rdmol object using UFF.
        Determines the energy level for each of the conformers identified in rdmol.GetConformer.


        :param rdmol:
        :param boundsMatrix:
        :param atomMatch:
        :return rdmol, minEid (index of the lowest energy conformer)
        """

        energy = 0.0
        minEid = 0
        lowestE = 9.999999e99  # start with a very high number, which would never be reached

        for conf in self.rdkit_molecule.GetConformers():
            if (self.bm is None) or (self.atom_match is None):
                AllChem.UFFOptimizeMolecule(self.rdkit_molecule, confId=conf.GetId())
                energy = AllChem.UFFGetMoleculeForceField(
                    self.rdkit_molecule, confId=conf.GetId()).CalcEnergy()
            else:
                _, energy = EmbedLib.OptimizeMol(
                    self.rdkit_molecule, self.bm, atomMatches=self.atom_match, forceConstant=100000.0)

            if energy < lowestE:
                minEid = conf.GetId()
                lowestE = energy

        return self.rdkit_molecule, minEid

    def rd_embed(self):
        """
        This portion of the script is literally taken from rmgpy but hacked to work without defining a geometry object

        Embed the RDKit molecule and create the crude molecule file.
        """
        numConfAttempts = 10000
        if (self.bm is None) or (self.atom_match is None):
            AllChem.EmbedMultipleConfs(self.rdkit_molecule, numConfAttempts, randomSeed=1)

            self.rdkit_molecule, minEid = self.optimize_rdkit_molecule()
        else:
            """
            Embed the molecule according to the bounds matrix. Built to handle possible failures
            of some of the embedding attempts.
            """
            self.rdkit_molecule.RemoveAllConformers()
            for i in range(0, numConfAttempts):
                try:
                    EmbedLib.EmbedMol(self.rdkit_molecule, self.bm, atomMatch=self.atom_match)
                    break
                except ValueError:
                    logging.info(
                        "RDKit failed to embed on attempt {0} of {1}".format(
                            i + 1, numConfAttempts))
                except RuntimeError:
                    logging.info("RDKit failed to embed.")
            else:
                logging.error("RDKit failed all attempts to embed")
                return None, None

            """
            RDKit currently embeds the conformers and sets the id as 0, so even though multiple
            conformers have been generated, only 1 can be called. Below the id's are resolved.
            """
            for i in range(len(self.rdkit_molecule.GetConformers())):
                self.rdkit_molecule.GetConformers()[i].SetId(i)

            self.rdkit_molecule, minEid = self.optimize_rdkit_molecule()

        return self.rdkit_molecule, minEid

    def get_bonds(self):
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        test_conf.rdkit_molecule = self._pseudo_geometry
        test_conf.ase_molecule = self.ase_molecule
        return test_conf.get_bonds()

    def get_torsions(self):
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        test_conf.rdkit_molecule = self._pseudo_geometry
        test_conf.ase_molecule = self.ase_molecule
        return test_conf.get_torsions()

    def get_angles(self):
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        test_conf.rdkit_molecule = self._pseudo_geometry
        test_conf.ase_molecule = self.ase_molecule
        return test_conf.get_angles()
