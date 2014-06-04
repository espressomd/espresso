# Define the espresso package
import sys, ctypes

# OpenMPI magic
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
