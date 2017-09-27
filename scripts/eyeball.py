#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py |github|
-------------------

Interactive eyeball planet visualizer. See :py:mod:`planetplanet.photo.eyeball`.

  .. plot::
     :align: center
     
     from scripts import eyeball
     eyeball._test()

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/eyeball.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet.photo import eyeball

def _test():
    '''

    '''

    eyeball.Interact()
  
if __name__ == '__main__':

    eyeball.Interact()
