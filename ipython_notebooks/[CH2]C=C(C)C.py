#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 9,
    'C': 5,
}

bonds = {
    'C-C': 3,
    'C=C': 1,
    'C-H': 9,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('[CH2]C=C(C)C.log'),
}

geometry = GaussianLog('[CH2]C=C(C)C.log')

frequencies = GaussianLog('[CH2]C=C(C)C.log')

rotors = []
