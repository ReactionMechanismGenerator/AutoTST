#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 9,
    'C': 5,
}

bonds = {
    'C=C': 1,
    'C-C': 3,
    'C-H': 9,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('ICPOKTZORRQGLK-UHFFFAOYSA.log'),
}

geometry = GaussianLog('ICPOKTZORRQGLK-UHFFFAOYSA.log')

frequencies = GaussianLog('ICPOKTZORRQGLK-UHFFFAOYSA.log')

rotors = []
