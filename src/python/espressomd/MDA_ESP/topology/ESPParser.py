from __future__ import print_function
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

class ESPParser(TopologyReader):
	"""Parse a list of atoms from a running ESPResSo simulation
	
	
	"""
	format = 'ESP'
	def parse(self):
	        """Access ESPResSo data and return the structure.
	
	        :Returns: MDAnalysis internal *structure* dict.
	        """
		espresso = self._u.kwargs['espresso']
	
		segid = "SYSTEM"
		resid = 1
		resname = "SYSTEM"
	
		atoms = []
		__id=0
		for p in espresso.part:
			name = "O"
	
			elem = guess_atom_element(name)
			mass = p.mass
			charge = p.q
	
			at = Atom(__id, name, elem, resname, resid,
				segid, mass, charge, universe=self._u)
			atoms.append(at)
			__id+=1
	
	        struc = {"atoms": atoms}
	
	        return struc
