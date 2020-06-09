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

import os, itertools, logging
import numpy as np
from copy import deepcopy
from typing import List, Set, Dict, Tuple, Optional, Callable, Iterator, Union

import rdkit
import rdkit.Chem 
import rdkit.Chem.AllChem
import rdkit.Chem.Pharm3D.EmbedLib
import rdkit.DistanceGeometry

import ase 
import ase.calculators.gaussian

import rmgpy
import rmgpy.molecule
import rmgpy.species 
import rmgpy.reaction 
import rmgpy.data.rmg 
import rmgpy.exceptions 

import autotst
from .data.base import DistanceData, TransitionStateDepository, TSGroups, TransitionStates
from .species import Species, Conformer
from .geometry import Torsion, Angle, Bond, CisTrans, ChiralCenter

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

try:
    import py3Dmol
except ImportError:
    logging.info("Error importing py3Dmol")

class Reaction():

    rmg_database = None   # will be an RMGDatabase instance, once loaded.
    # a dictionary will have reaction family names as keys and
    # TransitionStates instances as values, once loaded.
    ts_databases = dict()

    def __init__(
            self,
            label:Optional[str] = None,
            rmg_reaction:Optional[rmgpy.reaction.Reaction] = None,
            reaction_family:str = "H_Abstraction"):

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
        return f'<Reaction "{self.label}">'

    @property
    def ts(self)->Dict[str, TS]:
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
                    #smiles=complex.to_smiles()
                    direction=direction,
                    rmg_molecule=complex,
                    reaction_family=self.reaction_family,
                    distance_data=self.distance_data,
                    action=self.possible_actions[self.reaction_family][direction]
                )
                ts_dict[direction] = [ts]

            self._ts = ts_dict

        return self._ts

    @property
    def distance_data(self)->DistanceData:
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
    def load_databases(self, 
                      rmg_database_path:str = None, 
                      force_reload: str = False)-> Tuple[rmgpy.data.base.Database]:
        """
        Load the RMG and AutoTST databases, if they have not already been loaded,
        into the class level variables where they are stored.

        Variables:
        - force_reload (bool):if set to True then forces a reload, even if already loaded.
        - rmg_database_path (str): path to rmg database directory. If None, database will be 
          loaded from rmgpy.settings['database.directory']

        Returns:
        - Tuple of RMG Database classes
        """
        if self.rmg_database and self.ts_databases and not force_reload:
            return self.rmg_database, self.ts_databases

        rmg_database = rmgpy.data.rmg.RMGDatabase()

        if rmg_database_path is None:
	        database_path = rmgpy.settings['database.directory']
        else:
            database_path = rmg_database_path

        logging.info(f"Loading RMG database from '{database_path}'")

        self.possible_families = [  # These families (and only these) will be loaded from both RMG and AutoTST databases
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration"
        ]
        try:
            rmg_database.load(
                database_path,
                kinetics_families=self.possible_families,
                transport_libraries=[],
                reaction_libraries=[],
                seed_mechanisms=[],
                thermo_libraries=[
                    'primaryThermoLibrary',
                    'thermo_DFT_CCSDTF12_BAC',
                    'CBS_QB3_1dHR'],
                solvation=False,
            )
        except IOError:
            logging.info(
                f"RMG database not found at '{database_path}'. "
                "This can occur if a git repository of the database is being used rather than the binary version")
            database_path = os.path.join(database_path, 'input')
            logging.info(
                f"Loading RMG database instead from '{database_path}'")
            rmg_database.load(
                database_path,
                kinetics_families=self.possible_families,
                transport_libraries=[],
                reaction_libraries=[],
                seed_mechanisms=[],
                thermo_libraries=[
                    'primaryThermoLibrary',
                    'thermo_DFT_CCSDTF12_BAC',
                    'CBS_QB3_1dHR'],
                solvation=False,
            )

        self.rmg_database = rmg_database

        actions = {}
        for possible_family in self.possible_families:
            actions[possible_family] = {}
            action = self.rmg_database.kinetics.families[possible_family].forward_recipe.actions
            actions[possible_family]['forward'] = action

            action = self.rmg_database.kinetics.families[possible_family].reverse_recipe.actions
            actions[possible_family]['reverse'] = action

        self.possible_actions = actions

        self.ts_databases = dict()
        for reaction_family in self.possible_families:
            ts_database = TransitionStates()
            path = os.path.join(
                autotst.settings['tst_database_path'],
                reaction_family)
            global_context = {'__builtins__': None}
            local_context = {'DistanceData': DistanceData}
            family = self.rmg_database.kinetics.families[reaction_family]
            family_copy = deepcopy(family)
            ts_database.family = family_copy
            ts_database.load(path, local_context, global_context)

            self.ts_databases[reaction_family] = ts_database

        return self.rmg_database, self.ts_databases

    def generate_distance_data(self)->DistanceData:
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
        self._distance_data = self.ts_databases[self.reaction_family].groups.estimate_distances_using_group_additivity(
            self.rmg_reaction)
        if ((np.isclose(self._distance_data.distances["d12"] + self._distance_data.distances["d23"],
                      self._distance_data.distances["d13"],
                      atol=0.05)) or (self._distance_data.distances["d12"] + self._distance_data.distances["d23"] < self._distance_data.distances["d13"])):
            logging.info(
                "Distance between *1 and *3 is too small, setting it to lower bound of uncertainty")

            self._distance_data.distances["d13"] -= self._distance_data.uncertainties["d13"] / 2

        logging.info(f"The distance data is as follows: {self._distance_data}")

        return self._distance_data

    def generate_reactants_and_products(self)->Tuple[List[rmg.species.Species]]:
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

    def get_labeled_reaction(self)->Tuple[rmgpy.reaction.Reaction, str]:
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

        match = False # We have not found a labeled reaction that matches the one provided.
        if self.label: # A reaction label was provided
            
            # Generating lists of lists. Main list is the reactans or products
            # Secondary list is composed of the resonance structures for that species
            r, p = self.label.split("_") 

            smiles_conversions = {
                "[CH]":"[CH...]"
            }
            def get_rmg_mol(smile):
                if smile.upper() in list(smiles_conversions.keys()):
                    smile = smiles_conversions[smile.upper()]
                return rmgpy.molecule.Molecule(smiles=smile).generate_resonance_structures()
            rmg_reactants = [get_rmg_mol(smile) for smile in r.split("+")]
            rmg_products = [get_rmg_mol(smile) for smile in p.split("+")]

            combos_to_try = list(itertools.product(
                list(itertools.product(*rmg_reactants)),
                list(itertools.product(*rmg_products))
            ))

            # looping though each reaction family and each combination of reactants and products
            for name, family in list(self.rmg_database.kinetics.families.items()):
                logging.info(f"Trying to match reacction to {family}")
                for rmg_reactants, rmg_products in combos_to_try:
                    # Making a test reaction
                    test_reaction = rmgpy.reaction.Reaction(
                        reactants=list(rmg_reactants), 
                        products=list(rmg_products))

                    try: # Trying to label the reaction
                        labeled_r, labeled_p = family.get_labeled_reactants_and_products(
                            test_reaction.reactants, test_reaction.products)
                    except: # Failed to match a reaction to the family
                        logging.error(f"Couldn't match {test_reaction} to {name}, trying different combination...")
                        continue
                        
                    if not (labeled_r and labeled_p):
                        # Catching an error where it matched the reaction but the reactants and products we not return...
                        # A weird bug in RMG that I can explain
                        continue

                    if ((len(labeled_r) > 0) and (len(labeled_p) > 0)): # We found a match.
                        if match: 
                            # We found a match already, but we were able to label the reaction again... 
                            # This would happen if different combinations of resonance structures match up 
                            logging.warning("This reaction has been already labeled... \
                                it seems that resonance structures \for reactants and products exist.")
                            logging.warning("Skipping this duplicate for now")
                            continue

                        logging.info(f"Matched reaction to {name} family")
                        
                        # Copying labed reactions and saving them for later
                        labeled_reactants = deepcopy(labeled_r)
                        labeled_products = deepcopy(labeled_p)
                        match_reaction = rmgpy.reaction.Reaction(
                            reactants=deepcopy(labeled_r), 
                            products=deepcopy(labeled_p))
                        # Setting it true that we matche the reaction and remembering the
                        # RMG reaction family and family name 
                        match = True
                        final_family = family
                        final_name = name
                        

        elif self.rmg_reaction: #RMGReaction but no label
            rmg_reactants = []
            rmg_products = []
            for react in self.rmg_reaction.reactants:
                if isinstance(react, rmgpy.species.Species):
                    rmg_reactants.append(react.molecule)
                elif isinstance(react, rmgpy.molecule.Molecule):
                    rmg_reactants.append(react.generate_resonance_structures())

            for prod in self.rmg_reaction.products:
                if isinstance(prod, rmgpy.species.Species):
                    rmg_products.append(prod.molecule)
                elif isinstance(prod, rmgpy.molecule.Molecule):
                    rmg_products.append(prod.generate_resonance_structures())

            combos_to_try = list(itertools.product(
                list(itertools.product(*rmg_reactants)),
                list(itertools.product(*rmg_products))
            ))

            for name, family in list(self.rmg_database.kinetics.families.items()):
                logging.info(f"Trying to match reaction to {name}")
                for rmg_reactants, rmg_products in combos_to_try:
                    # Making a test reaction
                    test_reaction = rmgpy.reaction.Reaction(
                        reactants=list(rmg_reactants), 
                        products=list(rmg_products))
                    
                    try: # Trying to label the reaction
                        labeled_r, labeled_p = family.get_labeled_reactants_and_products(
                            test_reaction.reactants, test_reaction.products)
                        if not (labeled_r and labeled_p):
                            logging.error("Unable to determine a reaction for the forward direction. Trying the reverse direction.")
                            raise rmgpy.exceptions.ActionError
                    except:
                        try: 
                            # Trying the reverse reaction if the forward reaction doesn't work
                            # This is useful for R_Addition reactions
                            labeled_r, labeled_p = family.get_labeled_reactants_and_products(
                                test_reaction.products, test_reaction.reactants)
                        except:
                            logging.error(f"Couldn't match {test_reaction} to {name}, trying different combination...")
                            continue
                        
                    if not (labeled_r and labeled_p):
                        # Catching an error where it matched the reaction but the reactants and products we not return...
                        # A weird bug in RMG that I can explain
                        continue

                    if ((len(labeled_r) > 0) and (len(labeled_p) > 0)): # We found a match.
                        if match: 
                            # We found a match already, but we were able to label the reaction again... 
                            # This would happen if different combinations of resonance structures match up 
                            logging.warning("This reaction has been already labeled... \
                                it seems that resonance structures \for reactants and products exist.")
                            logging.warning("Skipping this duplicate for now")
                            continue

                        logging.info(f"Matched reaction to {name} family")
                        
                        # Copying labed reactions and saving them for later
                        labeled_reactants = deepcopy(labeled_r)
                        labeled_products = deepcopy(labeled_p)
                        match_reaction = rmgpy.reaction.Reaction(
                            reactants=deepcopy(labeled_r), 
                            products=deepcopy(labeled_p))
                        # Setting it true that we matche the reaction and remembering the
                        # RMG reaction family and family name 
                        match = True
                        final_family = family
                        final_name = name
                            
        assert match, "Could not idetify labeled reactants and products"

        reaction_list = final_family.generate_reactions(
            match_reaction.reactants, match_reaction.products)

        assert reaction_list, "Could not match a reaction to a reaction family..."

        for reaction in reaction_list:
            if match_reaction.is_isomorphic(reaction):
                reaction.reactants = labeled_reactants
                reaction.products = labeled_products
                break
        self.rmg_reaction = reaction
        self.reaction_family = final_name
        return self.rmg_reaction, self.reaction_family

    def get_label(self)->str:
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
            if isinstance(react, rmgpy.species.Species):
                string += f"{react.molecule[0].to_smiles()}+"
            elif isinstance(react, rmgpy.molecule.Molecule):
                string += f"{react.to_smiles()}+"
        string = string[:-1]
        string += "_"
        for prod in self.rmg_reaction.products:
            if isinstance(prod, rmgpy.species.Species):
                string += f"{prod.molecule[0].to_smiles()}+"
            elif isinstance(prod, rmgpy.molecule.Molecule):
                string += f"{prod.to_smiles()}+"
        self.label = string[:-1]
        return self.label

    def get_rmg_reaction(self)->rmgpy.reaction.Reaction:

        if self.rmg_reaction:
            return self.rmg_reaction
        
        r, p = self.label.split("_")

        reactants = []
        products = []
        for react in r.split("+"):
            reactants.append(rmgpy.molecule.Molecule(smiles=react))
        for prod in p.split("+"):
            products.append(rmgpy.molecule.Molecule(smiles=prod))

        self.rmg_reaction = rmgpy.reaction.Reaction(reactants=reactants, products=products)
        return self.rmg_reaction

    def get_rmg_complexes(self)->Dict[str, rmgpy.molecule.Molecule]:
        """
        A method to create a forward and reverse TS complexes used to initialize transition state geometries

        Variables:
        - rmg_reaction (RMGReaction): The RMGReaction of interest

        Returns:
        - complexes (dict): a dictionary containing rmgpy.molecule.Molecules of the forward and reverse reaction complexes
        """

        if self.rmg_reaction is None:
            self.get_labeled_reaction()

        reactant_complex = rmgpy.molecule.Molecule()
        for react in self.rmg_reaction.reactants:
            if isinstance(react, rmgpy.molecule.Molecule):
                reactant_complex = reactant_complex.merge(react)
            elif isinstance(react, rmgpy.species.Species):
                for mol in react.molecule:
                    if len(mol.get_all_labeled_atoms()) > 0:
                        reactant_complex = reactant_complex.merge(mol)

        product_complex = rmgpy.molecule.Molecule()
        for prod in self.rmg_reaction.products:
            if isinstance(prod, rmgpy.molecule.Molecule):
                product_complex = product_complex.merge(prod)
            elif isinstance(prod, rmgpy.species.Species):
                for mol in prod.molecule:
                    if len(mol.get_all_labeled_atoms()) > 0:
                        product_complex = product_complex.merge(mol)

        reactant_complex.update_multiplicity()
        product_complex.update_multiplicity()

        if len(reactant_complex.get_all_labeled_atoms()) == 0 or len(product_complex.get_all_labeled_atoms()) == 0:
            logging.warning("REACTING ATOMS LABELES NOT PROVIDED. Please call `Reaction.get_labeled_reaction` to generate labeled complexes")

        self.complexes = {
            "forward": reactant_complex,
            "reverse": product_complex}

        return self.complexes

    def generate_conformers(self, ase_calculator:Optional[ase.calculators.gaussian.Gaussian]=None)->Dict[str, TS]:
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
        
        from .conformer.systematic import systematic_search

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
            smiles:Optional[str] = None,
            reaction_label:Optional[str] = None,
            direction: str = 'forward',
            rmg_molecule:str = None,
            reaction_family:str = "H_Abstraction",
            distance_data:Optional[DistanceData] = None,
            action:List = [],
            index:int = 0):

        self.energy = None
        self.reaction_label = reaction_label
        self.direction = direction.lower()
        self.reaction_family = reaction_family
        self.distance_data = distance_data
        self.index = index
        self.bm = None
        self.action = action

        assert direction in ["forward",
                             "reverse"], "Please provide a valid direction"

        self._rdkit_molecule = None
        self._ase_molecule = None

        if (smiles or rmg_molecule):
            if smiles and rmg_molecule:
                assert rmg_molecule.is_isomorphic(rmgpy.molecule.Molecule(
                    smiles=smiles)), "smiles string did not match RMG Molecule object"
                self.smiles = smiles
                self.rmg_molecule = rmg_molecule

            elif rmg_molecule:
                self.rmg_molecule = rmg_molecule
                self.smiles = rmg_molecule.to_smiles()

            else:
                self.smiles = smiles
                self.rmg_molecule = rmgpy.molecule.Molecule(smiles=smiles)

            self.rmg_molecule.update_multiplicity()
            self._symmetry_number = None

        else:
            self.smiles = None
            self.rmg_molecule = None
            self._rdkit_molecule = None
            self._pseudo_geometry = None
            self._ase_molecule = None
            self.bonds = []
            self.angles = []
            self.torsions = []
            self.cistrans = []
            self.chiral_centers = []
            self._symmetry_number = None

    def __repr__(self):
        return f'<TS "{self.smiles}">'

    def copy(self):
        copy_conf = TS(
            reaction_label=self.reaction_label,
            reaction_family=self.reaction_family)
        copy_conf.smiles = self.smiles
        copy_conf.rmg_molecule = self.rmg_molecule.copy()
        copy_conf._rdkit_molecule = self.rdkit_molecule.__copy__()
        copy_conf._pseudo_geometry = self._pseudo_geometry.__copy__()
        copy_conf._ase_molecule = self.ase_molecule.copy()
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

    def get_rdkit_mol(self)->rdkit.Chem.Mol:
        """
        A method to create an rdkit geometry... slightly different than that of the conformer method
        returns both the rdkit_molecule and the bm
        """
        self._rdkit_molecule = Conformer(rmg_molecule=self.rmg_molecule).get_rdkit_mol()

        self.get_labels()
        for i, atom in enumerate(self.rmg_molecule.atoms):
            assert atom.number == self.rdkit_molecule.GetAtoms()[i].GetAtomicNum()

        if len(self.labels) == 3:

            rd_copy = rdkit.Chem.RWMol(self.rdkit_molecule.__copy__())

            for a in self.action:
                if 'BOND' in a[0]:
                    _, label1, _, label2 = a
                    sorting_label_1 = self.rmg_molecule.get_labeled_atoms(label1)[0].sorting_label
                    sorting_label_2 = self.rmg_molecule.get_labeled_atoms(label2)[0].sorting_label
                    if not rd_copy.GetBondBetweenAtoms(sorting_label_1, sorting_label_2):
                        rd_copy.AddBond(sorting_label_1, sorting_label_2, order=rdkit.Chem.rdchem.BondType.SINGLE)

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
            rdkit.DistanceGeometry.DoTriangleSmoothing(self.bm)

            logging.info("Now attempting to embed using edited bounds matrix.")

            self.rd_embed()
        return self.rdkit_molecule

    def get_bounds_matrix(self)->np.ndarray:
        """
        A method to obtain the bounds matrix
        """
        self.bm = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(self.rdkit_molecule)
        return self.bm

    def set_limits(self, lbl1: str, 
                   lbl2: str, 
                   value: float, 
                   uncertainty:float)->np.ndarray:
        """
        A method to set the limits of a particular distance between two atoms

        :param lbl1: the label of one atom
        :param lbl2: the label of another atom
        :param value: the distance from a distance data object (float)
        :param uncertainty: the uncertainty of the `value` distance (float)
        :return bm: an array of arrays corresponding to the edited bounds matrix
        """
        logging.info(
            f"For atoms {lbl1} and {lbl2} we have a distance of: \t {value}")
        if lbl1 > lbl2:
            self.bm[lbl2][lbl1] = value + uncertainty / 2
            self.bm[lbl1][lbl2] = max(0, value - uncertainty / 2)
        else:
            self.bm[lbl2][lbl1] = max(0, value - uncertainty / 2)
            self.bm[lbl1][lbl2] = value + uncertainty / 2

        return self.bm

    def bm_pre_edit(self, sect: List[int])->np.ndarray:
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
                    u_ik = self.bm[i, k] if k > i else self.bm[k, i]
                    u_kj = self.bm[j, k] if k > j else self.bm[k, j]

                    max_lij = u_ik + u_kj - 0.1
                    if self.bm[i, j] > max_lij:
                        logging.info(
                            f"Changing lower limit {self.bm[i, j]} to {max_lij}")
                        self.bm[i, j] = max_lij

        return self.bm

    def get_labels(self)->Tuple[List, List[Tuple]]:
        """
        A method to get the labeled atoms from a reaction

        :return labels: the atom labels corresponding to the reaction center
        :return atomMatch: a tuple of tuples the atoms labels corresponding to the reaction center
        """

        if len(self.rmg_molecule.get_all_labeled_atoms()) == 0:
            labels = []
            atom_match = ()

        if self.reaction_family.lower() in [
            'h_abstraction',
            'r_addition_multiplebond',
            'intra_h_migration']:
            # for i, atom in enumerate(reactants.atoms):
            lbl1 = self.rmg_molecule.get_all_labeled_atoms()["*1"].sorting_label
            lbl2 = self.rmg_molecule.get_all_labeled_atoms()["*2"].sorting_label
            lbl3 = self.rmg_molecule.get_all_labeled_atoms()["*3"].sorting_label
            labels = [lbl1, lbl2, lbl3]
            atom_match = ((lbl1,), (lbl2,), (lbl3,))
        elif self.reaction_family.lower() in ['disproportionation']:
            lbl1 = self.rmg_molecule.get_all_labeled_atoms()["*2"].sorting_label
            lbl2 = self.rmg_molecule.get_all_labeled_atoms()["*4"].sorting_label
            lbl3 = self.rmg_molecule.get_all_labeled_atoms()["*1"].sorting_label

            labels = [lbl1, lbl2, lbl3]
            atom_match = ((lbl1,), (lbl2,), (lbl3,))

        #logging.info(f"The labled atoms are {labels}.")
        self.labels = labels
        self.atom_match = atom_match
        return self.labels, self.atom_match

    def edit_matrix(self)->np.ndarray:
        """
        A method to edit the bounds matrix using labels and distance data
        """

        lbl1, lbl2, lbl3 = self.labels

        sect = []

        for atom in self.rmg_molecule.split()[0].atoms:
            sect.append(atom.sorting_label)

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

    def optimize_rdkit_molecule(self)->Tuple[rdkit.Chem.Mol,int]:
        """
        Optimizes the rdmol object using UFF.
        Determines the energy level for each of the conformers identified in rdmol.GetConformer.


        :return rdmol, minEid (index of the lowest energy conformer)
        """

        energy = 0.0
        min_eid = 0
        lowest_e = 9.999999e99  # start with a very high number, which would never be reached

        for conf in self.rdkit_molecule.GetConformers():
            if (self.bm is None) or (self.atom_match is None):
                rdkit.Chem.AllChem.UFFOptimizeMolecule(self._rdkit_molecule, confId=conf.GetId())
                energy = rdkit.Chem.AllChem.UFFGetMoleculeForceField(
                    self._rdkit_molecule, confId=conf.GetId()).CalcEnergy()
            else:
                _, energy = rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(
                    self._rdkit_molecule, self.bm, atomMatches=self.atom_match, forceConstant=100000.0)

            if energy < lowest_e:
                min_eid = conf.GetId()
                lowest_e = energy

        return self.rdkit_molecule, min_eid

    def rd_embed(self)->Union[Tuple[None], Tuple[rdkit.Chem.Mol,int]]:
        """
        This portion of the script is literally taken from rmgpy but hacked to work without defining a geometry object

        Embed the RDKit molecule and create the crude molecule file.
        """
        num_conf_attempts = 10000
        if (self.bm is None) or (self.atom_match is None):
            rdkit.Chem.AllChem.EmbedMultipleConfs(self._rdkit_molecule, num_conf_attempts, randomSeed=1)

            self._rdkit_molecule, minEid = self.optimize_rdkit_molecule()
        else:
            """
            Embed the molecule according to the bounds matrix. Built to handle possible failures
            of some of the embedding attempts.
            """
            self._rdkit_molecule.RemoveAllConformers()
            for i in range(0, num_conf_attempts):
                try:
                    rdkit.Chem.Pharm3D.EmbedLib.EmbedMol(self._rdkit_molecule, self.bm, atomMatch=self.atom_match)
                    break
                except ValueError:
                    logging.info(
                        f"RDKit failed to embed on attempt {i + 1} of {num_conf_attempts}")
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

            self._rdkit_molecule, min_eid = self.optimize_rdkit_molecule()

        return self._rdkit_molecule, min_eid

    def get_bonds(self):
        #Why are there two methods? : Krishna
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        try:
            test_conf._rdkit_molecule = self._pseudo_geometry
        except:
            self.get_rdkit_mol()
            test_conf._rdkit_molecule = self._pseudo_geometry
        test_conf._ase_molecule = self.ase_molecule

    def get_bonds(self)-> List[Bond]:
        """
        A method for identifying all of the bonds in a conformer
        """
        bond_list = []
        for bond in self._pseudo_geometry.GetBonds():
            bond_list.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

        bonds = []
        for index, indices in enumerate(bond_list):
            i, j = indices

            length = self.ase_molecule.get_distance(i, j)
            center = False

            if self.rmg_molecule.atoms[i].label != '' and self.rmg_molecule.atoms[j].label != '':
                for a in self.action:
                    if not 'BOND' in a[0]:
                        # Not a bond action
                        continue
                    elif self.rmg_molecule.atoms[i].label in a and self.rmg_molecule.atoms[j].label in a:
                        # Both atom labels are participating in a bond being formed or breaking
                        center = True

            bond = Bond(index=index,
                        atom_indices=indices,
                        length=length,
                        reaction_center=center)
            mask = self.get_mask(bond)
            bond.mask = mask

            bonds.append(bond)

        self.bonds = bonds

        return self.bonds

    def get_torsions(self)->List[Torsion]:
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        try:
	        test_conf._rdkit_molecule = self._pseudo_geometry
        except:
            self.get_rdkit_mol()
            test_conf._rdkit_molecule = self._pseudo_geometry
        test_conf._ase_molecule = self.ase_molecule
        return test_conf.get_torsions()

    def get_angles(self)->List[Angle]:
        test_conf = Conformer()
        test_conf.rmg_molecule = self.rmg_molecule
        try:
	        test_conf._rdkit_molecule = self._pseudo_geometry
        except:
            self.get_rdkit_mol()
            test_conf._rdkit_molecule = self._pseudo_geometry
        test_conf._ase_molecule = self.ase_molecule
        return test_conf.get_angles()
