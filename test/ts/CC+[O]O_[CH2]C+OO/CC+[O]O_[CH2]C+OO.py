#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'O': 2,
    'C': 2,
    'H': 7,
}

bonds = {
    'O-O': 1,
    'C-H': 6,
    'C-C': 1,
    'O-H': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': Log('CC+[O]O_[CH2]C+OO.log'),
}

geometry = Log('CC+[O]O_[CH2]C+OO.log')

frequencies = Log('CC+[O]O_[CH2]C+OO.log')

rotors = []

