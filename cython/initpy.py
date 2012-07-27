# ~/.pythonrc
# enable syntax completion
try:
    import readline
except ImportError:
    print "Module readline not available."
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")

import numpy
import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
import espresso as es
