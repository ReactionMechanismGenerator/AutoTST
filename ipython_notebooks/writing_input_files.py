
# coding: utf-8

# In[1]:


from autotst.reaction import *
from autotst.molecule import *
from autotst.geometry import *


# In[2]:


reactants = [Molecule(SMILES="CCCC"), Molecule(SMILES="[O]O")]
products = [Molecule(SMILES="[CH2]CCC"), Molecule(SMILES="OO")]
rmg_reaction = Reaction(reactants=reactants, products=products)
rmg_reaction


# In[3]:


test_reaction = AutoTST_Reaction(reaction_family="H_Abstraction", rmg_reaction=rmg_reaction)
test_reaction


# In[4]:


from ase.calculators.gaussian import *


# In[5]:


def rxn_shell_calc(autotst_rxn, scratch):

    indicies = []
    for i, atom in enumerate(autotst_rxn.ts.rmg_ts.atoms):
        if not (atom.label == ""):
            indicies.append(i)

    combos = []
    for combo in list(itertools.combinations(indicies, 2)):
        a,b = combo
        combos.append("{0} {1} F".format(a,b))


    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label= autotst_rxn.label + "_shell",
                    scratch=scratch,
                    method="m062x",
                    basis="6-311+g(2df,2p)",
                    extra="opt=(ts,calcfc,noeigentest) freq",
                    multiplicity = test_reaction.ts.rmg_ts.multiplicity,
                    addsec = combos)

    del calc.parameters['force']

    autotst_rxn.ts.ase_ts.set_calculator(calc)
    return autotst_rxn


def rxn_center_calc(autotst_rxn, scratch):

    indicies = []
    for i, atom in enumerate(autotst_rxn.ts.rmg_ts.atoms):
        if (atom.label == ""):
            indicies.append(i)

    combos = []
    for combo in list(itertools.combinations(indicies, 2)):
        a,b = combo
        combos.append("{0} {1} F".format(a,b))

    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label=autotst_rxn.label + "_center",
                    scratch=scratch,
                    method="m062x",
                    basis="6-311+g(2df,2p)",
                    extra="opt=(ts,calcfc,noeigentest)",
                    multiplicity = test_reaction.ts.rmg_ts.multiplicity,
                    addsec = combos)
    del calc.parameters['force']

    autotst_rxn.ts.ase_ts.set_calculator(calc)

    return autotst_rxn



def irc_calc(autotst_rxn, scratch):

    autotst_rxn.ts.rmg_ts.updateMultiplicity()

    calc = Gaussian(mem="5GB",
                    nprocshared="20",
                    label= autotst_rxn.label + "_irc",
                    scratch=scratch,
                    method="m062x",
                    basis="6-311+g(2df,2p)",
                    extra="irc=(calcall)",
                    multiplicity = test_reaction.ts.rmg_ts.multiplicity)
    del calc.parameters['force']

    autotst_rxn.ts.ase_ts.set_calculator(calc)

    return autotst_rxn



# In[6]:


rxn = rxn_shell_calc(test_reaction, ".")
calc1 = rxn.ts.ase_ts.get_calculator()
print calc1.label


# In[7]:


rxn = rxn_center_calc(test_reaction, ".")
calc2 = rxn.ts.ase_ts.get_calculator()
print calc2.label


# In[8]:


rxn = irc_calc(test_reaction, ".")
calc3 = rxn.ts.ase_ts.get_calculator()
print calc3.label


# In[9]:


calc1.calculate(test_reaction.ts.ase_ts)

test_reaction.ts.update_from_ase_ts()
