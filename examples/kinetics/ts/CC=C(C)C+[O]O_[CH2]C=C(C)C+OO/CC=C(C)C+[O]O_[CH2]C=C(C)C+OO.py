#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 11,
    'C': 5,
    'O': 2,
}

bonds = {
    'O-O': 2,
    'C-C': 6,
    'O-H': 2,
    'C=C': 2,
    'C-H': 20,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': Log('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO.log'),
}

geometry = Log('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO.log')

frequencies = Log('CC=C(C)C+[O]O_[CH2]C=C(C)C+OO.log')

rotors = []

