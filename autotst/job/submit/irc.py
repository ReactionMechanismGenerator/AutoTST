#!/usr/bin/env python3

import sys, os

from ase.io import read, write
from ase import Atoms, Atom
import itertools
from sella import MinModeAtoms, optimize
from ase.calculators.emt import EMT


### python3 autotst-sella.py
file_path = os.environ["FILE_PATH"]
calc_label = os.environ["CALC_LABEL"]
directory = os.environ["DIRECTORY"]

    
### Reading in the .xyz file containing the atoms of interest
atoms = read(file_path)
      

### Creating the calculator of interest
possible_calculators = [
    "gaussian",
    "nwchem"
]
assert calc_label in possible_calculators, "Calculator requested is not available in AutoTST"

if calc_label.lower() == "gaussian":
    from ase.calculators.gaussian import Gaussian as ASEGaussian
    calc = ASEGaussian(
        label=file_path.strip(".xyz"),
        method="m062x",
        basis="cc-pVTZ",
        directory=directory
    )
    calc.command = calc.command.replace("g09", "g16")
elif calc_label.lower() == "nwchem":
    from ase.calculators.nwchem import NWChem as ASENWChem
    calc = ASENWChem(
        label=file_path.strip(".xyz"),
        xc="m06-2x",
        basis="cc-pVTZ",
    )
    
#calc = EMT()
### Setting the optimization type and the tolerance based on the 

traj = file_path.replace("_input.xyz", ".traj")  
        
### Actually performing the optimization
myminmode = MinModeAtoms(atoms,  # Your Atoms object
                         calc,  # Your calculator
                         trajectory=traj,  # Your trajectory file
                         )

calculation = irc(myminmode,    # Your MinMode object
                  maxiter=1000, # Maximum number of force evaluations
                  ftol=1e-2,    #
                  dx=0.01,      # Norm of the force vector, convergence threshold
                  direction="both",   # Maximum residual for eigensolver convergence (should be <= 1)
                  )

reactants = myminmode.atoms[0]
products = myminmode.atoms[-1]

write(file_path.replace("_input.xyz","_reactants.xyz"),reactants)
write(file_path.replace("_input.xyz","_products.xyz"),reactants)


        
    