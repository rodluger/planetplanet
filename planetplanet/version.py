#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
version.py |github|
-------------------

Checks for updates to the code on `github`.

    .. role:: raw-html(raw)
         :format: html

    .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/version.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from . import __version__, __url__

def VersionCheck():
    '''
    Checks for updates to the code on `github` and prints a message
    to the screen if there's a new version available.
    
    '''
    
    try:
        
        # Imports
        from six.moves import urllib
        from distutils.version import StrictVersion
        import re
        
        # Attempt to read the __init__.py file from github
        r = urllib.request.Request(__url__)
        handler = urllib.request.urlopen(r, timeout = 1)
        code = handler.getcode()
        txt = handler.read().decode("ascii")

        # Success?
        if int(code) != 200 or "ERROR" in txt:
            raise Exception("Error retrieving current repo version.")

        # Get the version
        tmp = re.findall('__version__ = "([0-9.]*?)"', txt)
        if len(tmp) == 1:
            version = tmp[0]

        # Is an update available?
        if StrictVersion(version) > StrictVersion(__version__):
            print("[planetplanet] A new version of the code " + 
                  "is available (%s)." % version)
            print("[planetplanet] Run `pip install planetplanet " + 
                  "--upgrade` to upgrade.")

    except:

        # Fail silently on all exceptions
        pass