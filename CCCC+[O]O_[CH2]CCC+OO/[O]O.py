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
    'M06-2X/cc-pVTZ': GaussianLog('[O]O.log'),
}

geometry = GaussianLog('[O]O.log')

frequencies = GaussianLog('[O]O.log')

rotors = []
