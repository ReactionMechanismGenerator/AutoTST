import os
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import re
import imp
import itertools
import random
import numpy as np
from numpy import array
import pandas as pd

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem import AllChem

from rdkit.Chem.Pharm3D import EmbedLib

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction, _isomorphicSpeciesList
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

# AutoTST imports
from database import DistanceData, TransitionStateDepository, TSGroups, TransitionStates
from molecule import *
from geometry import *


rmg_database = RMGDatabase()
database_path = os.path.join(os.path.expanduser('~'), 'Code',  'RMG-database', 'input')
rmg_database.load(database_path,
                 kineticsFamilies=['H_Abstraction'],
                 transportLibraries=[],
                 reactionLibraries=[],
                 seedMechanisms=[],
                 thermoLibraries=['primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR' ],
                 solvation=False,
                 )

# TODO: Edit this so it works with multiple databases

ts_database = TransitionStates()
path = "../database/H_Abstraction"
global_context = { '__builtins__': None }
local_context={'DistanceData': DistanceData}
family = rmg_database.kinetics.families["H_Abstraction"]
ts_database.family = family
ts_database.load(path, local_context, global_context)



class AutoTST_Reaction():

    """
    A class to describe reactions in 3D in ase, rdkit and rmg.

    INPUTS:
    * label (str): a label to describe the reaction in the following format: r1+r2_p1+p2. Where r1, r2, p1, p2 are SMILES strings for the reactants and products.
    * reaction_family (str): a string to describe the rmg_reaction family
    * rmg_reaction (RMG Reaction) [optional]: an RMG Reaction object that describes the reaction of interest

    """


    def __init__(self, label=None, reaction_family=None, rmg_reaction=None):

        #assert label, "Please provde a reaction in the following format: r1+r2_p1+p2. Where r1, r2, p1, p2 are SMILES strings for the reactants and products."
        #assert reaction_family, "Please provide a reaction family."

        self.label = label
        self.reaction_family = reaction_family


        if rmg_reaction:

            reactant_mols = []
            product_mols = []

            reactants = rmg_reaction.reactants
            products = rmg_reaction.products

            for reactant in reactants:
                if type(reactant) == rmgpy.species.Species:
                    reactant_mols.append(AutoTST_Molecule(rmg_molecule=reactant.molecule[0]))
                elif type(reactant) == rmgpy.molecule.molecule.Molecule:
                    reactant_mols.append(AutoTST_Molecule(rmg_molecule=reactant))

            for product in products:
                if type(product) == rmgpy.species.Species:
                    product_mols.append(AutoTST_Molecule(rmg_molecule=product.molecule[0]))
                elif type(product) == rmgpy.molecule.molecule.Molecule:
                    product_mols.append(AutoTST_Molecule(rmg_molecule=product))

            self.reactant_mols = reactant_mols
            self.product_mols = product_mols
        elif label and reaction_family:
            self.get_reactants_and_products()


        self.get_rmg_reactions()
        self.create_ts_geometries()

    def get_reactants_and_products(self):

        """
        This uses the reaction label to creat multi_molecule objects of

        :return:
        """
        reactants, products = self.label.split("_")

        if "+" in reactants:
            reactants = reactants.split("+")

        if "+" in products:
            products = products.split("+")

        reactant_mols = []
        product_mols = []

        for reactant in reactants:
            reactant_mols.append(AutoTST_Molecule(reactant))

        for product in products:
            product_mols.append(AutoTST_Molecule(product))

        self.reactant_mols = reactant_mols
        self.product_mols = product_mols


    def get_rmg_reactions(self):
        """
        This method creates a labeled rmg_reaction from the reaction string

        the reaction string should look as follows: r1+r2_p1+p2
        """

        rmg_reactants = []
        rmg_products = []

        for reactant_mol in self.reactant_mols:
            rmg_reactants.append(reactant_mol.rmg_molecule)

        for product_mol in self.product_mols:
            rmg_products.append(product_mol.rmg_molecule)



        labeled_r, labeled_p = family.getLabeledReactantsAndProducts(rmg_reactants, rmg_products)

        test_reaction = Reaction(reactants=labeled_r, products=labeled_p, reversible=True)


        reaction_list = rmg_database.kinetics.generate_reactions_from_families(
            rmg_reactants,
            rmg_products)

        assert reaction_list

        for reaction in reaction_list:
            if reaction.isIsomorphic(test_reaction):
                if (_isomorphicSpeciesList(reaction.reactants, test_reaction.reactants)) and (_isomorphicSpeciesList(reaction.products, test_reaction.products)):
                    reaction.reactants = test_reaction.reactants
                    reaction.products = test_reaction.products

                elif (_isomorphicSpeciesList(reaction.products, test_reaction.reactants)) and (_isomorphicSpeciesList(reaction.reactants, test_reaction.products)):
                    reaction.products = test_reaction.reactants
                    reaction.reactants = test_reaction.products

        """
        #NOTE::: This was how we did it before, but now it should be working for master

        print rmg_reactants

        print reaction_list

        for reaction in reaction_list:
            print reaction
            # Check if any of the RMG proposed reactions matches the reaction in the mechanism
            if test_reaction.isIsomorphic(reaction):
                atom_labels_reactants = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
                atom_labels_products = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

                for reactant in reaction.reactants:
                    reactant.clearLabeledAtoms()
                    for atom in reactant.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_reactants[atom_label[0]] = True

                for product in reaction.products:
                    product.clearLabeledAtoms()
                    for atom in product.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_products[atom_label[0]] = True

                if all(atom_labels_reactants.values()) and all(atom_labels_products.values()):
                    # We successfully labeled all of the atoms
                    break
        """
        self.rmg_reaction = reaction
        self.distance_data = ts_database.groups.estimateDistancesUsingGroupAdditivity(reaction)

    def create_ts_geometries(self):
        """
        A method to use the tools in rmg / autotst to create a reasonable TS geometry
        This will create the geometry in both rdkit and ase

        :return:
        self.multi_ts: a multi_ts object that contains geometries of a ts in
                        rdkit, ase, and rmg molecules
        """
        self.ts = AutoTST_TS(self)


class AutoTST_TS():
    def __init__(self, autotst_reaction):

        self.autotst_reaction = autotst_reaction
        self.label = autotst_reaction.label  # make sure that both the reaction and TS have same label

        self.create_rdkit_ts_geometry()
        self.create_ase_ts_geometry()
        self.create_rmg_ts_geometry()
        self.get_ts_torsion_list()
        self.get_ts_torsions()
        self.get_ts_angle_indices()
        self.angle_indices
        self.get_ts_angle()
        i, j, k =self.angle_indices
        r_mask = self.angle.right_mask

        self.ase_ts.set_angle(a1=i, a2=j, a3=k, angle=float(180), mask=r_mask)
        labels = []
        for label, atom in self.rmg_ts.getLabeledAtoms().iteritems():
            labels.append(atom.sortingLabel)

        for torsion in self.torsions:
            if set(labels).issubset(torsion.indices[:3]):
                i, j, k, l = torsion.indices
                r_mask = torsion.right_mask

                self.ase_ts.set_dihedral(a1=i,
                                         a2=j,
                                         a3=k,
                                         a4=l,
                                         mask=r_mask,
                                         angle=float(180))
        self.update_from_ase_ts()

    def create_rdkit_ts_geometry(self):

        self.rmg_ts, product = self.setup_molecules()

        self.rmg_ts.updateMultiplicity()

        labels, atom_match = self.get_labels(self.rmg_ts)

        """combined = rdkit.Chem.Mol()

        for mol in self.autotst_reaction.reactant_mols:
            combined = rdkit.Chem.CombineMols(combined, mol.rdkit_molecule)"""

        combined = self.rmg_ts.toRDKitMol(removeHs=False)#, returnMapping=True)
        Chem.rdDistGeom.EmbedMolecule(combined)
        self.rdkit_ts = combined
        bm = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(combined)
        #print "Before editing"
        #print bm
        #print
        bm = self.edit_matrix(self.rmg_ts, bm, labels)
        #print "After editing"
        #print bm
        #print

        rdkit.DistanceGeometry.DoTriangleSmoothing(bm)
        #print "After smoothing"
        #print bm
        self.bm = bm

        self.rdkit_ts = self.rd_embed(self.rdkit_ts, 10000, bm=bm, match=atom_match)[0]

    def setup_molecules(self):

        merged_reacts = None
        merged_prods = None

        if len(self.autotst_reaction.rmg_reaction.reactants) == 2:
            merged_reacts = Molecule.merge(self.autotst_reaction.rmg_reaction.reactants[0],
                                           self.autotst_reaction.rmg_reaction.reactants[1])

        if len(self.autotst_reaction.rmg_reaction.products) == 2:
            merged_prods = Molecule.merge(self.autotst_reaction.rmg_reaction.products[0],
                                           self.autotst_reaction.rmg_reaction.products[1])

        return merged_reacts, merged_prods

    def get_labels(self, reactants):
        """
        A method to get the labeled atoms from a reaction

        :param reactants: a combined rmg_molecule object
        :return labels: the atom labels corresponding to the reaction center
        :return atomMatch: a tuple of tuples the atoms labels corresponding to the reaction center
        """

        if self.autotst_reaction.rmg_reaction.family.lower() in ['h_abstraction', 'r_addition_multiplebond', 'intra_h_migration']:
            for i, atom in enumerate(reactants.atoms):
                if atom.label == "*1":
                    lbl1 = i
                if atom.label == "*2":
                    lbl2 = i
                if atom.label == "*3":
                    lbl3 = i
            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,), (lbl2,), (lbl3,))
        elif self.autotst_reaction.rmg_reaction.family.lower() in ['disproportionation']:
            for i, atom in enumerate(reactants.atoms):
                if atom.label == "*2":
                    lbl1 = i
                if atom.label == "*4":
                    lbl2 = i
                if atom.label == "*1":
                    lbl3 = i
            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,), (lbl2,), (lbl3,))

        return labels, atomMatch

    def set_limits(self, bm, lbl1, lbl2, value, uncertainty):
        """
        A method to set the limits of a particular distance between two atoms

        :param bm: an array of arrays corresponding to the bounds matrix
        :param lbl1: the label of one atom
        :param lbl2: the label of another atom
        :param value: the distance from a distance data object (float)
        :param uncertainty: the uncertainty of the `value` distance (float)
        :return bm: an array of arrays corresponding to the edited bounds matrix
        """
        if lbl1 > lbl2:
            bm[lbl2][lbl1] = value + uncertainty / 2
            bm[lbl1][lbl2] = max(0, value - uncertainty / 2)
        else:
            bm[lbl2][lbl1] = max(0, value - uncertainty / 2)
            bm[lbl1][lbl2] = value + uncertainty / 2

        return bm


    def edit_matrix(self, rmg_ts, bm, labels):
        """
        A method to set the limits of the reaction center using the `set_limits` method

        :param rmg_ts: an rmg_molecule corresponding the the combined reactants or products (i.e. r1 and r2 or p1 and p2)
        :param bm: an array of arrays that corresponds to the bounds matrix
        :param labels: the atom indices corresponding to the reaction center
        :return bm: the edited bounds matrix
        """
        lbl1, lbl2, lbl3 = labels

        sect = []
        for atom in rmg_ts.split()[1].atoms: sect.append(atom.sortingLabel)

        uncertainties = {'d12': 0.02, 'd13': 0.02, 'd23': 0.02}  # distanceData.uncertainties or {'d12':0.1, 'd13':0.1, 'd23':0.1 } # default if uncertainty is None
        bm = self.set_limits(bm, lbl1, lbl2, self.autotst_reaction.distance_data.distances['d12'], uncertainties['d12'])
        bm = self.set_limits(bm, lbl2, lbl3, self.autotst_reaction.distance_data.distances['d23'], uncertainties['d23'])
        bm = self.set_limits(bm, lbl1, lbl3, self.autotst_reaction.distance_data.distances['d13'], uncertainties['d13'])

        bm = self.bm_pre_edit(bm, sect)
        #bm = (bm < 5) * bm + (bm >= 5) * 5
        return bm

    def bm_pre_edit(self, bm, sect):
        """
        Clean up some of the atom distance limits before attempting triangle smoothing.
        This ensures any edits made do not lead to unsolvable scenarios for the molecular
        embedding algorithm.

        sect is the list of atom indices belonging to one species.
        """
        others = range(len(bm))
        for idx in sect: others.remove(idx)

        for i in range(len(bm)):#sect:
            for j in range(i):#others:
                if i<j: continue
                for k in range(len(bm)):
                    if k==i or k==j or i==j: continue
                    Uik = bm[i,k] if k>i else bm[k,i]
                    Ukj = bm[j,k] if k>j else bm[k,j]

                    maxLij = Uik + Ukj - 0.1
                    if bm[i,j] >  maxLij:
                        logging.info("Changing lower limit {0} to {1}".format(bm[i, j], maxLij))
                        bm[i,j] = maxLij

        return bm

    def optimize(self, rdmol, boundsMatrix=None, atomMatch=None):
        """

        Optimizes the rdmol object using UFF.
        Determines the energy level for each of the conformers identified in rdmol.GetConformer.


        :param rdmol:
        :param boundsMatrix:
        :param atomMatch:
        :return rdmol, minEid (index of the lowest energy conformer)
        """

        energy = 0.0
        minEid = 0;
        lowestE = 9.999999e99;  # start with a very high number, which would never be reached
        crude = Chem.Mol(rdmol.ToBinary())

        for conf in rdmol.GetConformers():
            if boundsMatrix is None:
                AllChem.UFFOptimizeMolecule(rdmol, confId=conf.GetId())
                energy = AllChem.UFFGetMoleculeForceField(rdmol, confId=conf.GetId()).CalcEnergy()
            else:
                eBefore, energy = EmbedLib.OptimizeMol(rdmol, boundsMatrix, atomMatches=atomMatch,
                                                       forceConstant=100000.0)

            if energy < lowestE:
                minEid = conf.GetId()
                lowestE = energy

        return rdmol, minEid

    def rd_embed(self, rdmol, numConfAttempts, bm=None, match=None):
        """
        This portion of the script is literally taken from rmgpy but hacked to work without defining a geometry object

        Embed the RDKit molecule and create the crude molecule file.
        """
        if bm is None:  # bm = bounds matrix?
            AllChem.EmbedMultipleConfs(rdmol, numConfAttempts, randomSeed=1)

            rdmol, minEid = self.optimize(rdmol)
        else:
            """
            Embed the molecule according to the bounds matrix. Built to handle possible failures
            of some of the embedding attempts.
            """
            rdmol.RemoveAllConformers()
            for i in range(0, numConfAttempts):
                try:
                    EmbedLib.EmbedMol(rdmol, bm, atomMatch=match)
                    break
                except ValueError:
                    x = 3
                    logging.info("RDKit failed to embed on attempt {0} of {1}".format(i + 1, numConfAttempts))
                    # What to do next (what if they all fail?) !!!!!
                except RuntimeError:
                    logging.info("RDKit failed to embed.")
            else:
                logging.error("RDKit failed all attempts to embed")
                return None, None

            """
            RDKit currently embeds the conformers and sets the id as 0, so even though multiple
            conformers have been generated, only 1 can be called. Below the id's are resolved.
            """
            for i in range(len(rdmol.GetConformers())):
                rdmol.GetConformers()[i].SetId(i)


            rdmol, minEid = self.optimize(rdmol, boundsMatrix=bm, atomMatch=match)

        return rdmol, minEid

    def create_ase_ts_geometry(self):

        mol_list = AllChem.MolToMolBlock(self.rdkit_ts).split('\n')
        ase_atoms = []
        for i, line in enumerate(mol_list):

            if i > 3:

                try:
                    atom0, atom1, bond, rest = line
                    atom0 = int(atom0)
                    atom0 = int(atom1)
                    bond = float(bond)

                except ValueError:
                    try:
                        x, y, z, symbol = line.split()[0:4]
                        x = float(x)
                        y = float(y)
                        z = float(z)
                        # print symbol

                        ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

                    except:
                        continue

        self.ase_ts = Atoms(ase_atoms)

    def create_rmg_ts_geometry(self):

        for i, position in enumerate(self.ase_ts.get_positions()):
            self.rmg_ts.atoms[i].coords = position

    def view_ts(self, mol=None):
        """
        A method designed to create a 3D figure of the Multi_Molecule with py3Dmol
        """
        if not mol:
            mol = self.rdkit_ts

        mb = Chem.MolToMolBlock(mol)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mb, "sdf")
        p.setStyle({'stick':{}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

    def create_pseudo_geometry(self):

        rdmol_copy = self.rdkit_ts.__copy__()
        rdmol_copy = Chem.RWMol(rdmol_copy)
        for atom in rdmol_copy.GetAtoms():
            idx = atom.GetIdx()
            rmg_atom = self.rmg_ts.atoms[idx]

            if rmg_atom.label:
                if rmg_atom.label == "*1":
                    atom1_star = atom
                if rmg_atom.label == "*2":
                    atom2_star = atom
                if rmg_atom.label == "*3":
                    atom3_star = atom

        try:
            rdmol_copy.AddBond(atom1_star.GetIdx(), atom2_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)
        except RuntimeError:
            rdmol_copy.AddBond(atom2_star.GetIdx(), atom3_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)

        return rdmol_copy

    def get_ts_angle_indices(self):
        rdmol_copy = self.rdkit_ts.__copy__()
        rdmol_copy = Chem.RWMol(rdmol_copy)
        for atom in rdmol_copy.GetAtoms():
            idx = atom.GetIdx()
            rmg_atom = self.rmg_ts.atoms[idx]

            if rmg_atom.label:
                if rmg_atom.label == "*1":
                    atom1_star = atom
                if rmg_atom.label == "*2":
                    atom2_star = atom
                if rmg_atom.label == "*3":
                    atom3_star = atom

        if rdmol_copy.GetBondBetweenAtoms(atom1_star.GetIdx(), atom2_star.GetIdx()):
            bond_list = list(atom3_star.GetBonds())

            for bond in bond_list:
                atomX = bond.GetOtherAtom(atom3_star)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2_star.GetIdx():
                    other_atom = atomX
            angle_indices = [atom2_star.GetIdx(), atom3_star.GetIdx(), other_atom.GetIdx()]

        else:
            bond_list = list(atom1_star.GetBonds())

            for bond in bond_list:
                atomX = bond.GetOtherAtom(atom1_star)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2_star.GetIdx():
                    other_atom = atomX

            angle_indices = [atom2_star.GetIdx(), atom1_star.GetIdx(), other_atom.GetIdx()]

        self.angle_indices = angle_indices

    def get_ts_angle(self):

        i, j, k = self.angle_indices
        ang = self.ase_ts.get_angle(i, j, k)
        angle = Angle(indices=self.angle_indices, degree=ang, left_mask=[], right_mask=[])
        left_mask = self.get_ts_left_mask(angle)
        right_mask = self.get_ts_right_mask(angle)

        self.angle = Angle(self.angle_indices, ang, left_mask, right_mask)

        return self.angle

    def get_ts_torsion_list(self):

        rdmol_copy = self.create_pseudo_geometry()

        torsion_list = []
        for bond1 in rdmol_copy.GetBonds():
            atom1 = bond1.GetBeginAtom()
            atom2 = bond1.GetEndAtom()
            if atom1.IsInRing() or atom2.IsInRing():
                # Making sure that bond1 we're looking at are in a ring
                continue

            bond_list1 = list(atom1.GetBonds())
            bond_list2 = list(atom2.GetBonds())

            if not len(bond_list1) > 1 and not len(bond_list2) > 1:
                # Making sure that there are more than one bond attached to
                # the atoms we're looking at
                continue

            # Getting the 0th and 3rd atom and insuring that atoms
            # attached to the 1st and 2nd atom are not terminal hydrogens
            # We also make sure that all of the atoms are properly bound together

            # If the above are satisfied, we append a tuple of the torsion our torsion_list
            got_atom0 = False
            got_atom3 = False

            for bond0 in bond_list1:
                atomX = bond0.GetOtherAtom(atom1)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2.GetIdx():
                    got_atom0 = True
                    atom0 = atomX

            for bond2 in bond_list2:
                atomY = bond2.GetOtherAtom(atom2)
                if atomY.GetAtomicNum() == 1 and len(atomY.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    continue
                if atomY.GetIdx() != atom1.GetIdx():
                    got_atom3 = True
                    atom3 = atomY

            if not (got_atom0 and got_atom3):
                # Making sure atom0 and atom3 were not found
                continue

            # Looking to make sure that all of the atoms are properly bonded to eached
            if (
                            rdmol_copy.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
                            rdmol_copy.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
                        rdmol_copy.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())):
                torsion_tup = (atom0.GetIdx(), atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx())
                torsion_list.append(torsion_tup)

        self.torsion_list = torsion_list

    def get_ts_torsions(self):
        torsions = []
        for indices in self.torsion_list:
            i, j, k, l = indices

            dihedral = self.ase_ts.get_dihedral(i, j, k, l)
            tor = Torsion(indices=indices, dihedral=dihedral, left_mask=[], right_mask=[])
            left_mask = self.get_ts_left_mask(tor)
            right_mask = self.get_ts_right_mask(tor)

            torsions.append(Torsion(indices, dihedral, left_mask, right_mask))
        self.torsions = torsions
        return self.torsions

    def get_ts_right_mask(self, torsion_or_angle):

        rdmol_copy = self.create_pseudo_geometry()

        rdkit_atoms = rdmol_copy.GetAtoms()

        if "Torsion" in str(torsion_or_angle.__class__):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif "Angle" in str(torsion_or_angle.__class__):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]

        complete_RHS = False
        i = 0
        atom_index = RHS_atoms_index[0]
        while complete_RHS == False:
            try:
                RHS_atom = rdkit_atoms[atom_index]
                for neighbor in RHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in RHS_atoms_index) or (neighbor.GetIdx() in LHS_atoms_index):
                        continue
                    else:
                        RHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = RHS_atoms_index[i]

            except IndexError:
                complete_RHS = True

        right_mask = [index in RHS_atoms_index for index in range(len(self.ase_ts))]

        return right_mask

    def get_ts_left_mask(self, torsion_or_angle):

        rdmol_copy = self.create_pseudo_geometry()

        rdkit_atoms = rdmol_copy.GetAtoms()

        if "Torsion" in str(torsion_or_angle.__class__):

            L1, L0, R0, R1 = torsion_or_angle.indices

            # trying to get the left hand side of this torsion
            LHS_atoms_index = [L0, L1]
            RHS_atoms_index = [R0, R1]

        elif "Angle" in str(torsion_or_angle.__class__):
            a1, a2, a3 = torsion_or_angle.indices
            LHS_atoms_index = [a2, a1]
            RHS_atoms_index = [a2, a3]


        complete_LHS = False
        i = 0
        atom_index = LHS_atoms_index[0]
        while complete_LHS == False:
            try:
                LHS_atom = rdkit_atoms[atom_index]
                for neighbor in LHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in LHS_atoms_index) or (neighbor.GetIdx() in RHS_atoms_index):
                        continue
                    else:
                        LHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = LHS_atoms_index[i]

            except IndexError:
                complete_LHS = True

        left_mask = [index in LHS_atoms_index for index in range(len(self.ase_ts))]

        return left_mask

    def set_rmg_ts_coords(self, molecule_base):

        if molecule_base == "RDKit":
            mol_list = AllChem.MolToMolBlock(self.rdkit_ts).split('\n')
            for i, atom in enumerate(self.rmg_ts.atoms):

                j = i + 4
                coords = mol_list[j].split()[:3]

                for k, coord in enumerate(coords):

                    coords[k] = float(coord)
                atom.coords = np.array(coords)

        elif molecule_base == "ASE":
            for i, position in enumerate(self.ase_ts.get_positions()):
                self.rmg_ts.atoms[i].coords = position

    def update_from_rdkit_ts(self):
        # In order to update the ase molecule you simply need to rerun the get_ase_molecule method
        self.create_ase_ts_geometry()
        self.set_rmg_ts_coords("RDKit")

        # Getting the new torsion angles
        self.get_ts_torsions()

    def update_from_ase_ts(self):

        self.set_rmg_ts_coords("ASE")

        # setting the geometries of the rdkit molecule

        positions = self.ase_ts.get_positions()

        conf = self.rdkit_ts.GetConformers()[0]

        for i, atom in enumerate(self.rdkit_ts.GetAtoms()):
            conf.SetAtomPosition(i, positions[i])

        # Getting the new torsion angles
        self.get_ts_torsions()


    def update_from_rmg_ts(self):

        conf = self.rdkit_ts.GetConformers()[0]
        ase_atoms = []
        for i, atom in enumerate(self.rmg_ts.atoms):
            x, y, z = atom.coords
            symbol = atom.symbol

            conf.SetAtomPosition(i, [x, y, z])

            ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

        self.ase_ts = Atoms(ase_atoms)

        # Getting the new torsion angles
        self.get_ts_torsion_list()
