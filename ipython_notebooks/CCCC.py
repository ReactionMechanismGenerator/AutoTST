#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 10,
    'C': 4,
}

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'M062X/MG3S': GaussianLog('CCCC.log'),
}

geometry = GaussianLog('CCCC.log')

frequencies = GaussianLog('CCCC.log')

rotors = []
