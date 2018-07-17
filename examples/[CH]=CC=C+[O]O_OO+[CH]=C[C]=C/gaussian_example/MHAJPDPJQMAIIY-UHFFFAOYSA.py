#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 2,
    'O': 2,
}

bonds = {
    'O-H': 2,
    'O-O': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('MHAJPDPJQMAIIY-UHFFFAOYSA.log'),
}

geometry = GaussianLog('MHAJPDPJQMAIIY-UHFFFAOYSA.log')

frequencies = GaussianLog('MHAJPDPJQMAIIY-UHFFFAOYSA.log')

rotors = []
