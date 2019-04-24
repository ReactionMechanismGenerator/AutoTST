#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "M06-2X/cc-pVTZ"
frequencyScaleFactor = 0.982
useHinderedRotors = False
useBondCorrections = False

species('react_0', '../../species/CC=C(C)C/CC=C(C)C.py', structure=SMILES('CC=C(C)C'))
species('react_1', '../../species/[O]O/[O]O.py', structure=SMILES('[O]O'))
species('prod_0', '../../species/OO/OO.py', structure=SMILES('OO'))
species('prod_1', '../../species/[CH2]C=C(C)C/[CH2]C=C(C)C.py', structure=SMILES('[CH2]C=C(C)C'))
transitionState('TS', 'CC=C(C)C+[O]O_[CH2]C=C(C)C+OO.py')

reaction(
    label = 'CC=C(C)C+[O]O_[CH2]C=C(C)C+OO',
    reactants = ['react_0', 'react_1'],
    products = ['prod_0', 'prod_1'],
    transitionState = 'TS',
    tunneling = 'Eckart',
)

statmech('TS')
kinetics('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO')
statmech('react_0')
thermo('react_0', 'NASA')
statmech('react_1')
thermo('react_1', 'NASA')
statmech('prod_0')
thermo('prod_0', 'NASA')
statmech('prod_1')
thermo('prod_1', 'NASA')
