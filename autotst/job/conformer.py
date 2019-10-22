import sys, os
from autotst.species import Species, Conformer
from autotst.reaction import Reaction, TS
import pickle

try:
    from hotbit import Hotbit as ase_calculator
except:
    from ase.calculators.emt import EMT as ase_calculator

smiles = sys.argv[1]
directory = sys.argv[2]

if "_" in smiles:
    autotst_object = Reaction(smiles)
else:
    autotst_object = Species(smiles.split("/"))

autotst_object.generate_conformers(ase_calculator=ase_calculator())
if isinstance(autotst_object, Species):
    conformers = autotst_object.conformers
elif isinstance(autotst_object, Reaction):
    conformers = autotst_object.ts
    
save_dict = {}
for key, confs in conformers.items():
    key_dict = {}
    for conf in confs:
        key_dict[conf.index] = conf.ase_molecule.arrays

    save_dict[key] = key_dict

f = open(os.path.join(directory, "{}_conformers.pkl".format(smiles.replace(" / ", "_"))), "wb")
pickle.dump(save_dict, f)