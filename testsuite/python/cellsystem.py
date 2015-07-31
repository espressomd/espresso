
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
# Tests particle property setters/getters
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond



class CellSystem(ut.TestCase):
  S=espressomd.System()
  
  def test_cellSystem(self):
     self.S.cellSystem.setNsquare(useVerletLists=False)
     s=self.S.cellSystem.getState()
     self.assertEqual(s,{"useVerletLists":0,"type":"nsquare"})
     self.S.cellSystem.setLayered(nLayers=5)
     s=self.S.cellSystem.getState()
     self.assertEqual(s,{"type":"layered","nLayers":5})
     
     self.S.cellSystem.setDomainDecomposition(useVerletLists=True)
     s=self.S.cellSystem.getState()
     self.assertEqual(s,{"useVerletLists":1,"type":"domainDecomposition"})
     


if __name__ == "__main__":
 print("Features: ",espressomd.features())
 ut.main()
