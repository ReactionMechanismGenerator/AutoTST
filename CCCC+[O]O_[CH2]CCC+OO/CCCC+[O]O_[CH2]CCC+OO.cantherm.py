#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "M06-2X/cc-pVTZ"
frequencyScaleFactor = 0.982
useHinderedRotors = False
useBondCorrections = False

species('CCCC', 'CCCC.py')
species('[O]O', '[O]O.py')
species('[CH2]CCC', '[CH2]CCC.py')
species('OO', 'OO.py')
transitionState('TS', 'CCCC+[O]O_[CH2]CCC+OO.py')

reaction(
    label = 'CCCC+[O]O_[CH2]CCC+OO',
    reactants = ['CCCC', '[O]O'],
    products = ['[CH2]CCC', 'OO'],
    transitionState = 'TS',
    tunneling = 'Eckart',
)

statmech('TS')
kinetics('CCCC+[O]O_[CH2]CCC+OO')
