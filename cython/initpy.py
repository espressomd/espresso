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

# Creating a playground
es.glob.box_l = [10.,10.,10.]


# Setting a few particles
es.part[0].pos = [1.,0.,0.]
es.part[1].pos = [0.,1.,0.]
es.part[2].pos = [0.,0.,1.]
es.part[0].type= 0
es.part[1].type= 0
es.part[2].type= 0

