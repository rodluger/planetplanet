#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_scripts.py
---------------

Test all of the scripts.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import sys, os
SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts')
sys.path.insert(1, SCRIPTS)
import scripts

def test_all():
  '''
  Test all scripts in the `scripts/` directory.
  
  '''
  
  for script in dir(scripts):
    if not script.startswith('_'):    
      getattr(scripts, script)._test()

if __name__ == '__main__':
  test_all()