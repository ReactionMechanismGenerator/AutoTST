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
    'M06-2X/cc-pVTZ': GaussianLog('OUUQCZGPVNCOIJ-UHFFFAOYSA.log'),
}

geometry = GaussianLog('OUUQCZGPVNCOIJ-UHFFFAOYSA.log')

frequencies = GaussianLog('OUUQCZGPVNCOIJ-UHFFFAOYSA.log')

rotors = []
