#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 1,
    'O': 2,
}

bonds = {
    'O-H': 1,
    'O-O': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': Log('[O]O.log'),
}

geometry = Log('[O]O.log')

frequencies = Log('[O]O.log')

