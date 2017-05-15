#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
interact.py
-----------

Simple planet-planet occultation flux for a constant
brightness dayside and a constant brightness nightside.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.eyeball import Interact

Interact(n = 31)