import os
from autotst.reaction import *
from autotst.molecule import *
from autotst.geometry import *
from ase.io.gaussian import read_gaussian, read_gaussian_out

from autotst.calculators.vibrational_analysis import *
from autotst.calculators.gaussian import *

from ase.calculators.gaussian import *

from ase.optimize import BFGS
from ase.calculators.emt import *
from ase.calculators.gaussian import *

from autotst.conformer.utilities  import *
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *

reactants = [Molecule(SMILES="CC=C(C)C"), Molecule(SMILES="[O]O")]
products = [Molecule(SMILES="[CH2]C=C(C)C"), Molecule(SMILES="OO")]
rmg_reaction = Reaction(reactants=reactants, products=products)

test_reaction = AutoTST_Reaction(reaction_family="H_Abstraction", rmg_reaction=rmg_reaction)

test_reaction.ts.ase_ts.set_calculator(EMT())

mol = test_reaction

initial_pop = create_initial_population(test_reaction)

top_pop = perform_simple_es(test_reaction, initial_pop)

for i, dihedral in enumerate(top_pop.iloc[0,1:]):
    tor = mol.ts.torsions[i]
    mol.ts.ase_ts.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)
test_reaction.ts.update_from_ase_ts()

for i, reactant in enumerate(test_reaction.reactant_mols):

    reactant.ase_molecule.set_calculator(EMT())

    initial_pop = create_initial_population(reactant)
    top_pop = perform_simple_es(reactant, initial_pop)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)

    reactant.update_from_ase_mol()

for i, reactant in enumerate(test_reaction.product_mols):

    reactant.ase_molecule.set_calculator(EMT())

    initial_pop = create_initial_population(reactant)
    top_pop = perform_simple_es(reactant, initial_pop)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)

    reactant.update_from_ase_mol()

tst_calculators = AutoTST_Gaussian(test_reaction)

tst_calculators.run_all()

print "We did it???"
