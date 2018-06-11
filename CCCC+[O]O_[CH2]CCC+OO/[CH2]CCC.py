#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 9,
    'C': 4,
}

bonds = {
    'C-C': 3,
    'C-H': 9,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('[CH2]CCC.log'),
}

geometry = GaussianLog('[CH2]CCC.log')

frequencies = GaussianLog('[CH2]CCC.log')

rotors = []
