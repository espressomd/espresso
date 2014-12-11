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
include "myconfig.pxi"


def features():
  """Returns list of compiled-in features"""
  f=[]
  IF ELECTROSTATICS == 1:
     f.append("ELECTROSTATICS")
  IF DIPOLES == 1:
     f.append("DIPOLES") 
  IF LB_GPU == 1:
      f.append("LB_GPU")
  IF ROTATION == 1:
    f.append("ROTATION")
  IF MASS == 1 : 
    f.append("MASS")
  IF VIRTUAL_SITES == 1: 
    f.append("VIRTUAL_SITES")
  IF VIRTUAL_SITES_RELATIVE == 1: 
    f.append("VIRTUAL_SITES_RELATIVE")
   

  return sorted(f)

