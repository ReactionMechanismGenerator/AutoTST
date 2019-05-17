#!/usr/bin/env python3

import sys, os

from ase.io import read, write
from ase import Atoms, Atom
import itertools
from sella import MinModeAtoms, optimize
from ase.calculators.emt import EMT


### python3 autotst-sella.py Path/To/File.xyz CalcLabel OptType
file_path = os.environ["FILE_PATH"]
calc_label = os.environ["CALC_LABEL"]
opt_type = os.environ["OPT_TYPE"]
directory = os.environ["DIRECTORY"]
multiplicity = os.environ["MULT"]

atom_indicies = None
if opt_type in ["shell", "center"]:
    atom_indicies = os.environ["ATOM_INDICIES"][:-1].split(",")

#assert len(sys.argv) > 4, "To perform a partial optimization, you must provide atoms to freeze"

    
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
        directory=directory,
        multiplicity=int(multiplicity)
    )
    calc.command = calc.command.replace("g09", "g16")
elif calc_label.lower() == "nwchem":
    from ase.calculators.nwchem import NWChem as ASENWChem
    calc = ASENWChem(
        label=file_path.strip(".xyz"),
        xc="m062x",
        basis="cc-pVTZ",
        #directory=directory,
        multiplicity=int(multiplicity)
    )
    
#calc = EMT()
### Setting the optimization type and the tolerance based on the 
assert opt_type.lower() in ["overall", "shell", "center", "species"], "Please provide a valid TS optimization type"

if opt_type.lower() in ["overall", "center"]:
    order = 1
    ftol = 1e-2
else:
    order = 0
    ftol = 1e-1

    
### Obtaining the constraints
constraints = []
if atom_indicies:
    constrained_atoms = [int(index) for index in atom_indicies]
    for pair in itertools.permutations(constrained_atoms,2):
        if (pair in constraints) or (pair[::-1] in constraints):
            continue
        constraints.append(pair)
        
        
### Actually performing the optimization
myminmode = MinModeAtoms(atoms,  # Your Atoms object
                         calc,  # Your calculator
                         constraints=dict(bonds=constraints),  # Your constraints
                         )

x1 = optimize(myminmode,    # Your MinMode object
              maxiter=1000,  # Maximum number of force evaluations
              ftol=ftol,    # Norm of the force vector, convergence threshold
              r_trust=0.1,  # Initial trust radius (Angstrom)
              order=order,      # Order of saddle point to find (set to 0 for minimization)
              dxL=1e-4,     # Finite difference displacement magnitude (Angstrom)
              maxres=0.1,   # Maximum residual for eigensolver convergence (should be <= 1)
              )

write(file_path.replace("_input.xyz", ".xyz"), myminmode.atoms)


        
    