#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
import espresso as es
print dir(es)

th=es.thermostat.Thermostat()
print th.getStatus() # not initialized
th.turnOff()
print th.getStatus() # off
th.setLangevin(11., 12.)
print th.getStatus() # langevin

try:
    th.setLangevin(1,0.1)
except ValueError, err:
    print err # wrong args

try:
    th.setLangevin("lolailo",-1.)
except ValueError, err:
    print err # wrong args type/value
