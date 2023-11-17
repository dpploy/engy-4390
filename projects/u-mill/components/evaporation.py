#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""
Cortix module.
This module is a model of the Evaporator/Calcining process
"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

import unit

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

