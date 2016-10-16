import espressomd
from   MDAnalysis.lib.util import NamedStream
import cStringIO

from topology import ESPParser
from coordinates import ESP

class ESPStream(espressomd.System):
	def __init__(self):
		# add tests here 
		a=0

	def stream(self):
		# time
		__xyz=str(self.time)+'\n'
		# number of particles
		__xyz+=str(len(self.part))+'\n'
		# box edges
		__xyz+=str(self.box_l)+'\n'
		# configuration
		for __p in self.part:
			__xyz+=str(__p.pos)+'\n'
		return  NamedStream(cStringIO.StringIO(__xyz), "__.ESP")

	def system(self):
		return self
	



