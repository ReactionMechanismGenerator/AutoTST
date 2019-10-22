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
import itertools
import logging
import pandas as pd
import numpy as np
import os
import multiprocessing
from multiprocessing import Process, Manager
import time

import ase
from ase import Atoms
from ase import calculators
from ase.calculators.calculator import FileIOCalculator
from ase.optimize import BFGS
from ase import units
from ase.constraints import FixBondLengths

from rdkit.Chem import rdMolAlign

import autotst
from autotst.species import Conformer
from autotst.reaction import TS
from autotst.conformer.utilities import get_energy, find_terminal_torsions

from rmgpy.exceptions import AtomTypeError
from rmgpy.molecule import Molecule

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
    cistranss = conformer.cistrans
    chiral_centers = conformer.chiral_centers

    torsion_angles = np.arange(0, 360, delta)
    torsion_combos = list(itertools.product(
        torsion_angles, repeat=len(torsions)))

    if cistrans:
        cistrans_options = ["E", "Z"]
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
    - cistrans (bool): indication of if one wants to consider cistrans bonds
    - chiral_centers (bool): indication of if one wants to consider chiral centers bonds

    Returns:
    - confs (list): a list of unique `Conformer` objects within 1 kcal/mol of the lowest energy conformer determined
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
    manager = Manager()
    return_dict = manager.dict()
    pool = multiprocessing.Pool()

    def opt_conf(i, rmsd_cutoff):
        """
        A helper function to optimize the geometry of a conformer.
        Only for use within this parent function
        """
        conformer = conformers[i]

        calculator = conformer.ase_molecule.get_calculator()

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

        if isinstance(calc, FileIOCalculator):
            if calculator.directory:
                directory = calculator.directory 
            else: 
                directory = 'conformer_logs'
            calculator.label = "{}_{}".format(conformer.smiles, i)
            calculator.directory = os.path.join(directory, label,'{}_{}'.format(conformer.smiles, i))
            if not os.path.exists(calculator.directory):
                try:
                    os.makedirs(calculator.directory)
                except OSError:
                    logging.info("An error occured when creating {}".format(calculator.directory))

            calculator.atoms = conformer.ase_molecule

        conformer.ase_molecule.set_calculator(calculator)
        opt = BFGS(conformer.ase_molecule, logfile=None)

        if type == 'species':
            if isinstance(i,int):
                c = FixBondLengths(labels)
                conformer.ase_molecule.set_constraint(c)
            try:
                opt.run(steps=1e6)
            except RuntimeError:
                logging.info("Optimization failed...we will use the unconverged geometry")
                pass
            if str(i) == 'ref':
                conformer.update_coords_from("ase")
                try:
                    rmg_mol = Molecule()
                    rmg_mol.from_xyz(
                        conformer.ase_molecule.arrays["numbers"],
                        conformer.ase_molecule.arrays["positions"]
                    )
                    if not rmg_mol.is_isomorphic(reference_mol):
                        logging.info("{}_{} is not isomorphic with reference mol".format(conformer,str(i)))
                        return False
                except AtomTypeError:
                    logging.info("Could not create a RMG Molecule from optimized conformer coordinates...assuming not isomorphic")
                    return False
        
        if type == 'ts':
            c = FixBondLengths(labels)
            conformer.ase_molecule.set_constraint(c)
            try:
                opt.run(fmax=0.20, steps=1e6)
            except RuntimeError:
                logging.info("Optimization failed...we will use the unconverged geometry")
                pass

        conformer.update_coords_from("ase")  
        energy = get_energy(conformer)
        conformer.energy = energy
        if len(return_dict)>0:
            conformer_copy = conformer.copy()
            for index,post in return_dict.items():
                conf_copy = conformer.copy()
                conf_copy.ase_molecule.positions = post
                conf_copy.update_coords_from("ase")
                rmsd = rdMolAlign.GetBestRMS(conformer_copy.rdkit_molecule,conf_copy.rdkit_molecule)
                if rmsd <= rmsd_cutoff:
                    return True
        if str(i) != 'ref':
            return_dict[i] = conformer.ase_molecule.get_positions()
        return True

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
    if isinstance(calc, FileIOCalculator):
        logging.info("The calculator generates input and output files.")

    results = []
    global conformers
    conformers = {}
    combinations = {}
    logging.info("There are {} possible conformers to investigate...".format(len(combos)))
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


    processes = []
    for i, conf in list(conformers.items()):
        p = Process(target=opt_conf, args=(i, rmsd_cutoff))
        processes.append(p)

    active_processes = []
    for process in processes:
        if len(active_processes) < multiprocessing.cpu_count():
            process.start()
            active_processes.append(process)
            continue

        else:
            one_done = False
            while not one_done:
                for i, p in enumerate(active_processes):
                    if not p.is_alive():
                        one_done = True
                        break

            process.start()
            active_processes[i] = process
    complete = np.zeros_like(active_processes, dtype=bool)
    while not np.all(complete):
        for i, p in enumerate(active_processes):
            if not p.is_alive():
                complete[i] = True

    energies = []
    for positions in list(return_dict.values()):
        conf = conformer.copy()
        conf.ase_molecule.positions = positions
        conf.ase_molecule.set_calculator(calc)
        energy = conf.ase_molecule.get_potential_energy()
        conf.update_coords_from("ase")
        energies.append((conf,energy))

    df = pd.DataFrame(energies,columns=["conformer","energy"])
    df = df[df.energy < df.energy.min() + (energy_cutoff * units.kcal / units.mol /
            units.eV)].sort_values("energy").reset_index(drop=True)

    redundant = []
    conformer_copies = [conf.copy() for conf in df.conformer]
    for i,j in itertools.combinations(range(len(df.conformer)),2):
        copy_1 = conformer_copies[i].rdkit_molecule
        copy_2 = conformer_copies[j].rdkit_molecule
        rmsd = rdMolAlign.GetBestRMS(copy_1,copy_2)
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

    logging.info("We have identified {} unique, low-energy conformers for {}".format(
        len(confs), conformer))
    
    return confs