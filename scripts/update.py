#!/usr/bin/env python
# coding: utf-8

import os, datetime, logging
from pickle import dump, load
import autotst
from autotst.reaction import Reaction
import autotst.data.update as Update
import ase
import cclib.io
import rmgpy.data.base 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "-d", 
    "--directory", 
    help="The directory containing completed AutoTST jobs. This contains a 'ts' and 'species' directory.", 
    default=".")
parser.add_argument(
    "-m",
    "--method",
    help="The quantum chemistry method and basis used. In the format of 'method/basis'.",
    default="m062x/cc-pVTZ"
)
parser.add_argument(
    "-s",
    "--short-description",
    help="A short description of the source of the data that you're wanting to add.",
    default="m062x/cc-pVTZ calculation via group additive TS generator."
)
args = parser.parse_args()

# Families that can be updated
to_update = {
    "H_Abstraction":[],
    "R_Addition_MultipleBond":[],
    "intra_H_migration":[]
}

def read_log(file_path):
    """
    A helper method that allows one to easily parse log files

    file_path (str): a path to a log file that cclib can read

    returns: 
    Arrays (ase.Atoms().arrays, np.array): arrays that contain all of the relevant information
        needed to construct an ase Atoms object
    """
    symbol_dict = {
        35: "Br",
        17: "Cl",
        9:  "F",
        8:  "O",
        7:  "N",
        6:  "C",
        1:  "H",
    }
    atoms = []
    p = cclib.io.ccread(file_path, loglevel=logging.ERROR)
    for atom_num, coords in zip(p.atomnos, p.atomcoords[-1]):
        atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))

    return ase.Atoms(atoms).arrays

if __name__ == "__main__":
    scratch = args.directory
    for reaction_string in os.listdir(scratch):
        print("")
        print(reaction_string)
        if not os.path.exists(os.path.join(scratch, reaction_string, reaction_string + ".log")):
            print("\tLog file doesn't exist for {}, skipping...".format(reaction_string))
            continue
        if reaction_string.count("_") != 1:
            print("\t{} actually isn't a reaction string... hmmm...".format(reaction_string))
            continue

        if reaction_string.count("+") == 2:
            print("\tMatched to H_Abstraction")
            to_update["H_Abstraction"].append((reaction_string, read_log(os.path.join(scratch, reaction_string, reaction_string +".log" ))))
        elif reaction_string.count("+") == 1:
            print("\tMatched to R_Addition_MultipleBond")
            to_update["R_Addition_MultipleBond"].append((reaction_string, read_log(os.path.join(scratch, reaction_string, reaction_string +".log" ))))
        elif reaction_string.count("+") == 0:
            print("\tMatched to intra_H_migration")
            to_update["intra_H_migration"].append((reaction_string, read_log(os.path.join(scratch, reaction_string, reaction_string +".log" ))))

    print("We need to update the following number of reactions:")
    for key, values in to_update.items():
       print("\t- {:23} : {}".format(key, len(values)))

    for family, family_dict in to_update.items():
        #if family.lower() != "intra_h_migration": continue
        reactions_to_update = []
        print("Updating the {} reaction family...".format(family))
        print("")
        for reaction_string, arrays in family_dict:
            print("Updating {}".format(reaction_string))
            
            reaction = Reaction(reaction_string)
            reaction.get_labeled_reaction()
            ts = reaction.ts["forward"][0]

            if all(arrays["numbers"] == ts.ase_molecule.arrays["numbers"]): # Matches the forward direction
                ts._ase_molecule.arrays = arrays
            else: # Matched the reverse direction
                ts = reaction.ts["reverse"][0]
                ts._ase_molecule.arrays = arrays
                
            labels, _ = ts.get_labels()
            lbl1, lbl2, lbl3 = labels

            d12 = ts.ase_molecule.get_distance(lbl1, lbl2)
            d23 = ts.ase_molecule.get_distance(lbl2, lbl3)
            d13 = ts.ase_molecule.get_distance(lbl1, lbl3)

            reaction.distance_data.distances["d12"] = d12
            reaction.distance_data.distances["d23"] = d23
            reaction.distance_data.distances["d13"] = d13

            reaction.distance_data.comment = "Calculated from AutoTST ({}), m062x/cc-pVTZ".format(datetime.date.today())
            reactions_to_update.append(reaction)
            
        Update.update_all(reactions_to_update, family, method=args.method, short_desc='m062x/cc-pVTZ calculation via group additive TS generator.')
