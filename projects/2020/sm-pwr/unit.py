#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

from scipy.constants import *
second = 1.0
kg = kilo*gram
meter = 1.0
joule = 1.0
pascal = 1.0
watt = 1.0
kelvin = 1.0

cm = centi*meter
ft = foot
cc = (centi*meter)**3
kj = kilo*joule
btu = Btu
barn = 1.0e-28 * meter**2
F = convert_temperature(2,'F','K') - convert_temperature(1,'F','K')
C = convert_temperature(2,'C','K') - convert_temperature(1,'C','K')
K = 1.0
