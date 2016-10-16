from __future__ import print_function
import espressomd
from espressomd import code_info
import os
import numpy as np
import MDAnalysis as mda

from MDAnalysis.topology.base import TopologyReader
from MDAnalysis.coordinates.base import SingleFrameReader
from MDAnalysis.lib import util
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors
from MDAnalysis.coordinates import base
from MDAnalysis.topology.core import get_atom_mass, guess_atom_charge, guess_atom_element
import cStringIO


class Timestep(base.Timestep):
    _ts_order_x = [0, 3, 4]
    _ts_order_y = [5, 1, 6]
    _ts_order_z = [7, 8, 2]

    def _init_unitcell(self):
        return np.zeros(9, dtype=np.float32)

    @property
    def dimensions(self):
        # This information now stored as _ts_order_x/y/z to keep DRY
        x = self._unitcell[self._ts_order_x]
        y = self._unitcell[self._ts_order_y]
        z = self._unitcell[self._ts_order_z]  # this ordering is correct! (checked it, OB)
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        x, y, z = triclinic_vectors(box)
        np.put(self._unitcell, self._ts_order_x, x)
        np.put(self._unitcell, self._ts_order_y, y)
 

class ESPReader(SingleFrameReader):
	format = 'ESP'
	units = {'time': None, 'length': 'nm', 'velocity': 'nm/ps'}
	_Timestep = Timestep
	
	def _read_first_frame(self):
		with util.openany(self.filename, 'rt') as espfile:
			self.n_atoms = len(espfile)
			for pos, line in enumerate(espfile,start=-3):
				if  (pos==-3):
					time = float(line[1:-1])
				elif(pos==-2):
					self.n_atoms = int(line)
					self.ts = ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
					self.ts.time = time
				elif(pos==-1):
					self.ts._unitcell[:3] = np.array(map(float,line[1:-2].split()))
				else:
 					ts._pos[pos] = np.array(list(map(float, line[1:-2].split())))
			ts.positions=ts._pos


