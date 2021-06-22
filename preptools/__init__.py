#!/usr/bin/env python

# global variables
__version__= "1.0"
HEADER="""prep-tools, Version %s"""%(__version__)
print(HEADER)
name = "preptools"

import os
try:
    AMBERHOME = os.environ['AMBERHOME']
except:
    AMBERHOME = None


# from preptools import wrappers
# from preptools import fileclasses

__all__ = ['wrappers', 'fileclasses']

from preptools import test_preptools 
# These tests of installed packages and pathnames are run every time the module is imported,
# so you thepathname preptools.AMBERHOME, e.g. is guaranteed to be valid.




