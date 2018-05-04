#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 10,
    'C': 5,
}

bonds = {
    'C=C': 1,
    'C-C': 3,
    'C-H': 10,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('CC=C(C)C.log'),
}

geometry = GaussianLog('CC=C(C)C.log')

frequencies = GaussianLog('CC=C(C)C.log')

rotors = []
