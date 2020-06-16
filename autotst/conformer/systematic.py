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

import itertools, logging, os, time
import pandas as pd
import numpy as np
import multiprocessing
from copy import deepcopy

import ase
import ase.units
import ase.calculators.calculator 
import ase.optimize
import ase.constraints

import rdkit.Chem

import rmgpy.exceptions
import rmgpy.molecule

import autotst
from ..species import Conformer
from ..reaction import TS
from .utilities import get_energy, find_terminal_torsions

def find_all_combos(
        conformer,
        delta=float(120),
        cistrans=True,
        chiral_centers=True):
    """
    A function to find all possible conformer combinations for a given conformer

    Params:
    - conformer (`Conformer`) an AutoTST `Conformer` object of interest
    - delta (int or float): a number between 0 and 180 or how many conformers to generate per dihedral
    - cistrans (bool): indication of if one wants to consider cistrans bonds
    - chiral_centers (bool): indication of if one wants to consider chiral centers bonds

    Returns:
    - all_combos (list): a list corresponding to the number of unique conformers created.
    """

    conformer.get_geometries()

    _, torsions = find_terminal_torsions(conformer)
    chiral_centers = conformer.chiral_centers

    torsion_angles = np.arange(0, 360, delta)
    torsion_combos = list(itertools.product(
        torsion_angles, repeat=len(torsions)))

    if cistrans:
        cistranss = []
        cistrans_options = ["E", "Z"]
        try:
            ring_info = conformer._pseudo_geometry.GetRingInfo()
        except AttributeError:
            ring_info = conformer.rdkit_molecule.GetRingInfo()

        for cistrans in conformer.cistrans:
            i,j,k,l = cistrans.atom_indices
            if (ring_info.NumAtomRings(i) != 0) or (ring_info.NumAtomRings(k) != 0):
                continue
            cistranss.append(cistrans)

        cistrans_combos = list(itertools.product(
            cistrans_options, repeat=len(cistranss)))

    else:
        cistrans_combos = [()]

    if chiral_centers:
        chiral_options = ["R", "S"]
        chiral_combos = list(itertools.product(
            chiral_options, repeat=len(chiral_centers)))

    else:
        chiral_combos = [()]

    all_combos = list(
        itertools.product(
            torsion_combos,
            cistrans_combos,
            chiral_combos))
    return all_combos


def opt_conf(i):
        """
        A helper function to optimize the geometry of a conformer.
        param i: index of the conformer or a conformer object
        """
        try:   
            conformer = conformers[i] #use the global object
        except:
            # When running tests for this single function, it's hard to create a global 
            # conformers dict. This step allows users to pass in conformer objects rather
            # than specify a global dict
            conformer = i
            
        if not isinstance(conformer, TS):
            reference_mol = conformer.rmg_molecule.copy(deep=True)
            reference_mol = reference_mol.to_single_bonds()
        calculator = conformer.ase_molecule.get_calculator()
        calculator.__init__()
        calculator = deepcopy(calculator)
        labels = []
        for bond in conformer.get_bonds():
            labels.append(bond.atom_indices)
    
        if isinstance(conformer, TS):
            label = conformer.reaction_label
            ind1 = conformer.rmg_molecule.get_labeled_atoms("*1")[0].sorting_label
            ind2 = conformer.rmg_molecule.get_labeled_atoms("*3")[0].sorting_label
            labels.append([ind1, ind2])
            type = 'ts'
        else:
            label = conformer.smiles
            type = 'species'

        if isinstance(calculator, ase.calculators.calculator.FileIOCalculator):
            if calculator.directory:
                directory = calculator.directory 
            else: 
                directory = 'conformer_logs'
            calculator.label = f"{conformer.smiles}_{conformer.index}"
            calculator.directory = os.path.join(directory, label,f'{conformer.smiles}_{conformer.index}')
            if not os.path.exists(calculator.directory):
                try:
                    os.makedirs(calculator.directory)
                except OSError:
                    logging.info(f"An error occured when creating {calculator.directory}")

            calculator.atoms = conformer.ase_molecule
        conformer.ase_molecule.set_calculator(calculator)
        opt = ase.optimize.BFGS(conformer.ase_molecule, logfile=None)
        if type == 'species':
            if isinstance(conformer.index,int):
                c = ase.constraints.FixBondLengths(labels)
                conformer.ase_molecule.set_constraint(c)
            try:
                opt.run(steps=1e6)
            except RuntimeError:
                logging.info("Optimization failed...we will use the unconverged geometry")
                pass
            if str(conformer.index) == 'ref':
                conformer.update_coords_from("ase")
                try:
                    rmg_mol = rmgpy.molecule.Molecule()
                    rmg_mol.from_xyz(
                        conformer.ase_molecule.arrays["numbers"],
                        conformer.ase_molecule.arrays["positions"]
                    )
                    if not rmg_mol.is_isomorphic(reference_mol):
                        logging.info(f"{conformer}_{str(conformer.index)} is not isomorphic with reference mol")
                        return False
                except rmgpy.exceptions.AtomTypeError:
                    logging.info("Could not create a RMG Molecule from optimized conformer coordinates...assuming not isomorphic")
                    return False
        converged = False
        if type == 'ts':
            c = ase.constraints.FixBondLengths(labels)
            conformer.ase_molecule.set_constraint(c)
            try:
                opt.run(fmax=0.20, steps=1e6)
            except RuntimeError:
                logging.info("Optimization failed...we will use the unconverged geometry")
                converged = True
                pass
       
        conformer.update_coords_from("ase")  
        try:
            energy = conformer.ase_molecule.get_potential_energy()
        except:
            if not converged:
                logging.error("Unable to parse energy from unconverged geometry")
            else:
                logging.error("Unable to parse energy from geometry")
            energy = 1e5
        try:
            conformers[i] = conformer #update the conformer from old object
        except:
            logging.error('Could not add updated conformer to conformers dict')
        return energy #return energy
def systematic_search(conformer,
                      delta=float(120),
                      energy_cutoff = 10.0, #kcal/mol
                      rmsd_cutoff = 0.5, #angstroms
                      cistrans = True,
                      chiral_centers = True,
                      multiplicity = False,
                      ):
    """
    Perfoms a systematic conformer analysis of a `Conformer` or a `TS` object

    Variables:
    - conformer (`Conformer` or `TS`): a `Conformer` or `TS` object of interest
    - delta (int or float): a number between 0 and 180 or how many conformers to generate per dihedral
    - energy_cutoff (str or float): energy in kcal/mol 
    - rmsd_cutoff (str or float): root mean square deviation of inter atomic positions 
    - cistrans (bool): indication of if one wants to consider cistrans bonds
    - chiral_centers (bool): indication of if one wants to consider chiral centers bonds

    Returns:
    - confs (list): a list of unique `Conformer` objects within 10 kcal/mol of the lowest energy conformer determined
    """
    
    rmsd_cutoff_options = {
        'loose' : 1.0,
        'default': 0.5,
        'tight': 0.1
    }

    energy_cutoff_options = {
        'high' : 50.0,
        'default' : 10.0,
        'low' : 5.0
    }

    if isinstance(rmsd_cutoff,str):
        rmsd_cutoff = rmsd_cutoff.lower()
        assert rmsd_cutoff in rmsd_cutoff_options.keys(), 'rmsd_cutoff options are loose, default, and tight'
        rmsd_cutoff = rmsd_cutoff_options[rmsd_cutoff]

    if isinstance(energy_cutoff,str):
        energy_cutoff = energy_cutoff.lower()
        assert energy_cutoff in energy_cutoff_options.keys(), 'energy_cutoff options are low, default, and high'
        energy_cutoff = energy_cutoff_options[energy_cutoff]
    
    if not isinstance(conformer, TS):
        reference_mol = conformer.rmg_molecule.copy(deep=True)
        reference_mol = reference_mol.to_single_bonds()
   
    #if not isinstance(conformer,TS):
    #    calc = conformer.ase_molecule.get_calculator()
    #    reference_conformer = conformer.copy()
    #    if opt_conf(reference_conformer, calc, 'ref', rmsd_cutoff):
    #        conformer = reference_conformer

    combos = find_all_combos(
        conformer,
        delta=delta,
        cistrans=cistrans,
        chiral_centers=chiral_centers)

    if len(combos) == 0:
        logging.info(
            "This species has no torsions, cistrans bonds, or chiral centers")
        logging.info("Returning origional conformer")
        return [conformer]

    _, torsions = find_terminal_torsions(conformer)

    calc = conformer.ase_molecule.get_calculator()
    if isinstance(calc, ase.calculators.calculator.FileIOCalculator):
        logging.info("The calculator generates input and output files.")

    results = []
    global conformers
    conformers = {}
    combinations = {}
    logging.info(f"There are {len(combos)} possible conformers to investigate...")
    for index, combo in enumerate(combos):

        combinations[index] = combo

        torsions, cistrans, chiral_centers = combo
        copy_conf = conformer.copy()

        for i, torsion in enumerate(torsions):

            tor = copy_conf.torsions[i]
            i, j, k, l = tor.atom_indices
            mask = tor.mask

            copy_conf.ase_molecule.set_dihedral(
                a1=i,
                a2=j,
                a3=k,
                a4=l,
                angle=torsion,
                mask=mask
            )
            copy_conf.update_coords()

        for i, e_z in enumerate(cistrans):
            ct = copy_conf.cistrans[i]
            copy_conf.set_cistrans(ct.index, e_z)

        for i, s_r in enumerate(chiral_centers):
            center = copy_conf.chiral_centers[i]
            copy_conf.set_chirality(center.index, s_r)

        copy_conf.update_coords_from("ase")
        copy_conf.ase_molecule.set_calculator(calc)
  
        conformers[index] = copy_conf

    num_threads = multiprocessing.cpu_count() - 1 or 1
    pool = multiprocessing.Pool(processes=num_threads)
    """
    to_calculate_list = []
    for i, conformer in list(conformers.items()):
        to_calculate_list.append(conformer)
    """
    results = pool.map(opt_conf,range(len(conformers)))
    pool.close()
    pool.join()
    energies = []
    for i,energy in enumerate(results):
        energies.append((conformers[i],energy))

    df = pd.DataFrame(energies,columns=["conformer","energy"])
    df = df[df.energy < df.energy.min() + (energy_cutoff * ase.units.kcal / ase.units.mol /
            ase.units.eV)].sort_values("energy").reset_index(drop=True)

    redundant = []
    conformer_copies = [conf.copy() for conf in df.conformer]
    for i,j in itertools.combinations(range(len(df.conformer)),2):
        copy_1 = conformer_copies[i].rdkit_molecule
        copy_2 = conformer_copies[j].rdkit_molecule
        rmsd = rdkit.Chem.rdMolAlign.GetBestRMS(copy_1,copy_2)
        if rmsd <= rmsd_cutoff:
            redundant.append(j)

    redundant = list(set(redundant))
    df.drop(df.index[redundant], inplace=True)

    if multiplicity and conformer.rmg_molecule.multiplicity > 2:
        rads = conformer.rmg_molecule.get_radical_count()
        if rads % 2 == 0:
            multiplicities = range(1,rads+2,2)
        else:
            multiplicities = range(2,rads+2,2)
    else:
        multiplicities = [conformer.rmg_molecule.multiplicity]

    confs = []
    i = 0
    for conf in df.conformer:
        if multiplicity:
            for mult in multiplicities:
                conf_copy = conf.copy()
                conf_copy.index = i
                conf_copy.rmg_molecule.multiplicity = mult
                confs.append(conf_copy)
                i += 1
        else:
            conf.index = i
            confs.append(conf)
            i += 1

    logging.info(f"We have identified {len(confs)} unique, low-energy conformers for {conformer}")
    
    return confs
