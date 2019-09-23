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

from rmgpy.exceptions import AtomTypeError
from rmgpy.molecule import Molecule

import autotst
from autotst.species import Conformer
from autotst.reaction import TS
from autotst.conformer.utilities import get_energy, find_terminal_torsions

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
                      rmsd_cutoff = 1.0, #angstroms
                      cistrans = True,
                      chiral_centers = True,
                      multiplicity = True
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
    # Takes each of the molecule objects
    
    if not isinstance(conformer, TS):
        reference_mol = conformer.rmg_molecule.copy(deep=True)
        reference_mol = reference_mol.toSingleBonds()
    manager = Manager()
    return_dict = manager.dict()

    def opt_conf(conformer, calculator, i, rmsd_cutoff):
        """
        A helper function to optimize the geometry of a conformer.
        Only for use within this parent function
        """

        labels = []
        for bond in conformer.get_bonds():
            labels.append(bond.atom_indices)
    
        if isinstance(conformer, TS):
            label = conformer.reaction_label
            ind1 = conformer.rmg_molecule.getLabeledAtom("*1").sortingLabel
            ind2 = conformer.rmg_molecule.getLabeledAtom("*3").sortingLabel
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
                    rmg_mol.fromXYZ(
                        conformer.ase_molecule.arrays["numbers"],
                        conformer.ase_molecule.arrays["positions"]
                    )
                    if not rmg_mol.isIsomorphic(reference_mol):
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
            for index,conf in return_dict.items():
                conf_copy = conf.copy()
                rmsd = rdMolAlign.GetBestRMS(conformer_copy.rdkit_molecule,conf_copy.rdkit_molecule)
                if rmsd <= rmsd_cutoff:
                    return True
        if str(i) != 'ref':
            return_dict[i] = conformer
        return True

    if not isinstance(conformer,TS):
        calc = conformer.ase_molecule.get_calculator()
        reference_conformer = conformer.copy()
        if opt_conf(reference_conformer, calc, 'ref', rmsd_cutoff):
            conformer = reference_conformer

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
    conformers = {}
    combinations = {}
    for index, combo in enumerate(combos):

        combinations[index] = combo

        torsions, cistrans, chiral_centers = combo

        for i, torsion in enumerate(torsions):

            tor = conformer.torsions[i]
            i, j, k, l = tor.atom_indices
            mask = tor.mask

            conformer.ase_molecule.set_dihedral(
                a1=i,
                a2=j,
                a3=k,
                a4=l,
                angle=torsion,
                mask=mask
            )
            conformer.update_coords()

        for i, e_z in enumerate(cistrans):
            ct = conformer.cistrans[i]
            conformer.set_cistrans(ct.index, e_z)

        for i, s_r in enumerate(chiral_centers):
            center = conformer.chiral_centers[i]
            conformer.set_chirality(center.index, s_r)

        conformer.update_coords_from("ase")
  
        conformers[index] = conformer.copy()


    processes = []
    for i, conf in list(conformers.items()):
        p = Process(target=opt_conf, args=(conf, calc, i, rmsd_cutoff))
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
    for conf in list(return_dict.values()):
        energies.append((conf,conf.energy))

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

    if conformer.rmg_molecule.multiplicity > 2:
        rads = conformer.rmg_molecule.getRadicalCount()
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