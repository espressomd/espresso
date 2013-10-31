#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
import espresso as es
print dir(es)

cs=es.cellsystem.Cellsystem()
gh=es.global_variables.GlobalsHandle()

# domain decomposition with verlet list: three equivalent commands
cs.setDomainDecomposition()
cs.setDomainDecomposition(True)
cs.setDomainDecomposition(useVerletList=True)


