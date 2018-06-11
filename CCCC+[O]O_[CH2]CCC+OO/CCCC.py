#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 10,
    'C': 4,
}

bonds = {
    'C-C': 3,
    'C-H': 10,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('CCCC.log'),
}

geometry = GaussianLog('CCCC.log')

frequencies = GaussianLog('CCCC.log')

rotors = []
