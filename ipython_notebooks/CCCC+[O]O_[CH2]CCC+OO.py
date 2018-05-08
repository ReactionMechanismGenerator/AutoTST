#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 11,
    'C': 4,
    'O': 2,
}

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M062X/MG3S': GaussianLog('CCCC+[O]O_[CH2]CCC+OO.log'),
}

geometry = GaussianLog('CCCC+[O]O_[CH2]CCC+OO.log')

frequencies = GaussianLog('CCCC+[O]O_[CH2]CCC+OO.log')

rotors = []

