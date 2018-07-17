#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 10,
    'C': 5,
}

bonds = {
    'C=C': 1,
    'C-C': 3,
    'C-H': 10,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'M06-2X/cc-pVTZ': GaussianLog('BKOOMYPCSUNDGP-UHFFFAOYSA.log'),
}

geometry = GaussianLog('BKOOMYPCSUNDGP-UHFFFAOYSA.log')

frequencies = GaussianLog('BKOOMYPCSUNDGP-UHFFFAOYSA.log')

rotors = []
