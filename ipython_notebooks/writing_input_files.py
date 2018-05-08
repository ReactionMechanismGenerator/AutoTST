import os

from autotst.reaction import *
from autotst.molecule import *
from autotst.geometry import *
from ase.io.gaussian import read_gaussian, read_gaussian_out

from ase.optimize import BFGS
from ase.calculators.emt import *
from ase.calculators.gaussian import *

from autotst.conformer.utilities  import *
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *


### Creating rmg_reaction

reactants = [Molecule(SMILES="CC=C(C)C"), Molecule(SMILES="[O]O")]
products = [Molecule(SMILES="[CH2]C=C(C)C"), Molecule(SMILES="OO")]
rmg_reaction = Reaction(reactants=reactants, products=products)

test_reaction = AutoTST_Reaction(reaction_family="H_Abstraction", rmg_reaction=rmg_reaction)

### Conformer Analysis

test_reaction.ts.ase_ts.set_calculator(EMT())

mol = test_reaction

initial_pop = create_initial_population(test_reaction)

top_pop = perform_ga(test_reaction, initial_pop, tolerance=1e-5)

for i, dihedral in enumerate(top_pop.iloc[0,1:]):
    tor = mol.ts.torsions[i]
    mol.ts.ase_ts.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)
test_reaction.ts.update_from_ase_ts()


### Getting reactant and product geometries

def reactants_or_products_calc(autotst_mol, scratch, label):

    autotst_mol.rmg_molecule.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label=label,
                    scratch=scratch,
                    extra="opt=(verytight,gdiis) freq IOP(2/16=3)",
                    multiplicity = autotst_mol.rmg_molecule.multiplicity
                    )
    del calc.parameters['force']

    return calc


for i, reactant in enumerate(test_reaction.reactant_mols):

    label = "r" + str(i)
    reactant.ase_molecule.set_calculator(EMT())

    initial_pop = create_initial_population(reactant)
    top_pop = perform_ga(reactant, initial_pop, tolerance=1e-5)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)

    calc = reactants_or_products_calc(reactant, "/gss_gpfs_scratch/harms.n/QMkin", label)
    calc.calculate(reactant.ase_molecule)

    reactant.ase_molecule = read_gaussian_out(label + ".log")
    reactant.update_from_ase_mol()

    os.rename("./"+label + ".log", "./" + reactant.rmg_molecule.toSMILES() + ".log")
    os.rename("./"+label + ".ase", "./" + reactant.rmg_molecule.toSMILES() + ".ase")
    os.rename("./"+label + ".com", "./" + reactant.rmg_molecule.toSMILES() + ".com")



for i, reactant in enumerate(test_reaction.product_mols):

    label = "p" + str(i)
    reactant.ase_molecule.set_calculator(EMT())

    initial_pop = create_initial_population(reactant)
    top_pop = perform_ga(reactant, initial_pop, tolerance=1e-5)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = reactant.torsions[i]
        reactant.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)

    calc = reactants_or_products_calc(reactant, "/gss_gpfs_scratch/harms.n/QMkin", label)
    calc.calculate(reactant.ase_molecule)
    reactant.ase_molecule = read_gaussian_out(label + ".log")
    reactant.update_from_ase_mol()

    os.rename("./"+label + ".log", "./" + reactant.rmg_molecule.toSMILES() + ".log")
    os.rename("./"+label + ".ase", "./" + reactant.rmg_molecule.toSMILES() + ".ase")
    os.rename("./"+label + ".com", "./" + reactant.rmg_molecule.toSMILES() + ".com")


### Making partial optimization calculators

def rxn_shell_calc(autotst_rxn, scratch):

    indicies = []
    for i, atom in enumerate(autotst_rxn.ts.rmg_ts.atoms):
        if not (atom.label == ""):
            indicies.append(i)

    combos = ""
    for combo in list(itertools.combinations(indicies, 2)):
        a,b = combo
        combos += "{0} {1} F\n".format(a+1,b+1)

    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label="rxn_shell",
                    scratch=scratch,
                    extra="Opt=(ModRedun,Loose) Int(Grid=SG1)",
                    multiplicity = autotst_rxn.ts.rmg_ts.multiplicity,
                    addsec = [combos[:-1]])

    del calc.parameters['force']
    return calc


def rxn_center_calc(autotst_rxn, scratch):

    indicies = []
    for i, atom in enumerate(autotst_rxn.ts.rmg_ts.atoms):
        if (atom.label == ""):
            indicies.append(i)

    combos = ""
    for combo in list(itertools.combinations(indicies, 2)):
        a,b = combo
        combos += "{0} {1} F\n".format(a+1,b+1)

    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label="rxn_center",
                    scratch=scratch,
                    extra="opt=(ts,calcfc,noeigentest,ModRedun)",
                    multiplicity = autotst_rxn.ts.rmg_ts.multiplicity,
                    addsec = [combos[:-1]])

    del calc.parameters['force']
    return calc


def ts_calc(autotst_rxn, scratch):

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label="rxn_overall",
                    scratch=scratch,
                    extra="opt=(ts,calcfc,noeigentest) freq",
                    multiplicity = autotst_rxn.ts.rmg_ts.multiplicity)

    del calc.parameters['force']
    return calc


def irc_calc(autotst_rxn, scratch):

    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label="rxn_irc",
                    scratch=scratch,
                    extra="irc=(calcall)",
                    multiplicity = test_reaction.ts.rmg_ts.multiplicity)

    del calc.parameters['force']
    return calc


calc1 = rxn_shell_calc(test_reaction, "/gss_gpfs_scratch/harms.n/QMkin")
calc2 = rxn_center_calc(test_reaction, "/gss_gpfs_scratch/harms.n/QMkin")
calc3 = ts_calc(test_reaction, "/gss_gpfs_scratch/harms.n/QMkin")
calc4 = irc_calc(test_reaction, "/gss_gpfs_scratch/harms.n/QMkin")


### Shell Optimization
print "Running {}".format(calc1.label)
calc1.calculate(test_reaction.ts.ase_ts)
test_reaction.ts.ase_ts = read_gaussian_out(calc1.label+ ".log")

os.rename("./rxn_shell.log", "./" + test_reaction.label + "_shell.log")
os.rename("./rxn_shell.ase", "./" + test_reaction.label + "_shell.ase")
os.rename("./rxn_shell.com", "./" + test_reaction.label + "_shell.com")

### Center Optimization
print "Running {}".format(calc2.label + ".log")
calc2.calculate(test_reaction.ts.ase_ts)
test_reaction.ts.ase_ts = read_gaussian_out(calc2.label+ ".log")

os.rename("./rxn_center.log", "./" + test_reaction.label + "_center.log")
os.rename("./rxn_center.ase", "./" + test_reaction.label + "_center.ase")
os.rename("./rxn_center.com", "./" + test_reaction.label + "_center.com")

### Overall Optimization
print "Running {}".format(calc3.label)
calc3.calculate(test_reaction.ts.ase_ts)
test_reaction.ts.ase_ts = read_gaussian_out(calc3.label+ ".log")

os.rename("./rxn_overall.log", "./" + test_reaction.label + "_overall.log")
os.rename("./rxn_overall.ase", "./" + test_reaction.label + "_overall.ase")
os.rename("./rxn_overall.com", "./" + test_reaction.label + "_overall.com")

### Running the IRC calc
print "Running {}".format(calc4.label)
try:
    calc4.calculate(test_reaction.ts.ase_ts)
except:
    # This normally fails because of an issue with ase's `read_results` method.
    pass

### Validating the IRC Calc
f = open(calc4.label + ".log", "r")
file_lines = f.readlines()[-5:]

validated_irc = False
for file_line in file_lines:
    if " Normal termination" in file_line: validated_irc = True

os.rename("./rxn_irc.log", "./" + test_reaction.label + "_irc.log")
os.rename("./rxn_irc.ase", "./" + test_reaction.label + "_irc.ase")
os.rename("./rxn_irc.com", "./" + test_reaction.label + "_irc.com")


### Running CanTherm
if validated_irc:
    test_reaction.ts.update_from_ase_ts()
    print "Complete!!!"
