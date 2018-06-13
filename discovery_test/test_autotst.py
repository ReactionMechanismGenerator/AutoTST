# To create TS geometry guesses

import numpy as np
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem.Pharm3D import EmbedLib

import ase

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction, _isomorphicSpeciesList, ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

import os
from autotst.reaction import AutoTST_Reaction, AutoTST_TS
from autotst.molecule import AutoTST_Molecule
from autotst.geometry import Bond, Angle, Torsion
from ase.io.gaussian import read_gaussian_out

# To perform TS search
from ase.calculators.gaussian import Gaussian
from autotst.calculators.gaussian import AutoTST_Gaussian
from autotst.calculators.vibrational_analysis import Vibrational_Analysis, percent_change
from autotst.calculators.cantherm import AutoTST_CanTherm

# For conformer analysis
from ase.calculators.lj import LennardJones
from autotst.conformer.utilities  import create_initial_population, select_top_population
from autotst.conformer.ga import perform_ga
from autotst.conformer.simple_es import perform_simple_es

reactants = [Molecule(SMILES="CCCC"), Molecule(SMILES="[O]O")]
products = [Molecule(SMILES="[CH2]CCC"), Molecule(SMILES="OO")]
rmg_reaction = Reaction(reactants=reactants, products=products)

test_reaction = AutoTST_Reaction("CC=C(C)C+[O]O_[CH2]C=C(C)C+OO", "H_Abstraction")
"""
### Performing the conformer analyses ###
## On the reaction ##
test_reaction.ts.ase_ts.set_calculator(LennardJones())
mol = test_reaction
initial_pop = create_initial_population(test_reaction)
top_pop = perform_ga(test_reaction, initial_pop)

for i, dihedral in enumerate(top_pop.iloc[0,1:]):
    tor = mol.ts.torsions[i]
    mol.ts.ase_ts.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)
test_reaction.ts.update_from_ase_ts()

## On the reactants ##
for i, reactant in enumerate(test_reaction.reactant_mols):
    reactant.ase_molecule.set_calculator(LennardJones())
    initial_pop = create_initial_population(reactant)
    top_pop = perform_ga(reactant, initial_pop)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)

    reactant.update_from_ase_mol()

## On the products ##
for i, reactant in enumerate(test_reaction.product_mols):
    reactant.ase_molecule.set_calculator(LennardJones())
    initial_pop = create_initial_population(reactant)
    top_pop = perform_ga(reactant, initial_pop)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)
    reactant.update_from_ase_mol()
"""

### Performing the partial optimizations
tst_calculators = AutoTST_Gaussian(test_reaction, scratch="test")
gaussian_results = tst_calculators.run_all(vibrational_analysis=False)

if gaussian_results:
    ### Running CanTherm ###
    cantherm = AutoTST_CanTherm(tst_calculators.reaction, scratch="test", output_directory="test")
    cantherm.write_files()
    cantherm.run()
    cantherm.set_reactants_and_products()

    ### Printing Results ###
    logging.info("The kinetics of intrest are as follows:")
    logging.info("{0!r}".format(cantherm.kinetics_job.reaction))

else:
    logging.info("Failed gaussian... :(")
