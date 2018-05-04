#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 11,
    'C': 5,
    'O': 2,
}

bonds = {
    'O-O': 1,
    'C=C': 1,
    'O-H': 1,
    'C-C': 3,
    'C-H': 10,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO_overall.log'),
}

geometry = GaussianLog('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO_overall.log')

frequencies = GaussianLog('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO_overall.log')

rotors = []

