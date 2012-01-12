#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

print sys.path

import espresso as es
print dir(es)
