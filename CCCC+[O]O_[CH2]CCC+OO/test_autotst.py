# To create TS geometry guesses
import os
from autotst.reaction import *
from autotst.molecule import *
from autotst.geometry import *
from ase.io.gaussian import read_gaussian, read_gaussian_out

# To perform TS search
from ase.calculators.gaussian import *
from autotst.calculators.gaussian import *
from autotst.calculators.vibrational_analysis import *
from autotst.calculators.cantherm import *

# For conformer analysis
from ase.calculators.lj import *
from autotst.conformer.utilities  import *
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *

reactants = [Molecule(SMILES="CCCC"), Molecule(SMILES="[O]O")]
products = [Molecule(SMILES="[CH2]CCC"), Molecule(SMILES="OO")]
rmg_reaction = Reaction(reactants=reactants, products=products)

test_reaction = AutoTST_Reaction(reaction_family="H_Abstraction", rmg_reaction=rmg_reaction)

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

### Performing the partial optimizations
tst_calculators = AutoTST_Gaussian(test_reaction)
tst_calculators.run_all()

### Running CanTherm ###
cantherm = AutoTST_CanTherm(test_reaction)
cantherm.write_files()
cantherm.run()
cantherm.set_reactants_and_products()

### Printing Results ###
logging.info("{0!r}".format(cantherm.kinetics_job.reaction))
