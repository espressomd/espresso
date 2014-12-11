#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013,2014 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
from __future__ import print_function
import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))
import espresso as es
print(dir(es))

cs=es.cellsystem.Cellsystem()
gh=es.global_variables.GlobalsHandle()

# domain decomposition with verlet list: three equivalent commands
cs.setDomainDecomposition()
cs.setDomainDecomposition(True)
cs.setDomainDecomposition(useVerletList=True)


