#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
import espresso as es
print dir(es)

cs=es.cellsystem.Cellsystem()
gh=es.global_variables.GlobalsHandle()
cs.setDomainDecomposition()
gh.skin=2.
cs.setDomainDecomposition()
