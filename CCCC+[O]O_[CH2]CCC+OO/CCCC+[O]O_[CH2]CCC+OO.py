#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 11,
    'C': 4,
    'O': 2,
}

bonds = {
    'C-C': 3,
    'O-H': 1,
    'O-O': 1,
    'C-H': 10,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('CCCC+[O]O_[CH2]CCC+OO_overall.log'),
}

geometry = GaussianLog('CCCC+[O]O_[CH2]CCC+OO_overall.log')

frequencies = GaussianLog('CCCC+[O]O_[CH2]CCC+OO_overall.log')

rotors = []

