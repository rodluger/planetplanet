#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py
----------

Interactive eyeball planet visualizer.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball

eyeball.Interact()