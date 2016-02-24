#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
import espresso as es
print(dir(es))

th=es.thermostat.Thermostat()
print(th.getStatus()) # not initialized
th.turnOff()
print(th.getStatus()) # off
th.setLangevin(11., 12.)
print(th.getStatus()) # langevin

try:
    th.setLangevin(1,0.1)
except ValueError as err:
    print(err) # wrong args

try:
    th.setLangevin("lolailo",-1.)
except ValueError as err:
    print(err) # wrong args type/value
