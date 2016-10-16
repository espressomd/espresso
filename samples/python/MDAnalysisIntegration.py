## This example shows how to integrate the functionalities of MDAnalysis in ESPResSo
import espressomd
from espressomd import MDA_ESP 
import numpy as np
import MDAnalysis as mda

### set up a minimal sample system 

system = espressomd.System()
system.time_step   		= 0.001
system.cell_system.skin         = 0.1

system.box_l 			= [10,10,10]
for i in range(10):
	system.part.add(id=i, pos=np.random.random(3) * system.box_l)
for i in range(5):
	system.part[i].q = 1.0

### prepare the stream

eos = MDA_ESP.ESPStream()

### MDAnalysis Universe:  we must pass the additional 'espresso' keyword
u =  mda.Universe( eos.stream(), espresso=eos.system() )

### Write the config as a pdb file
print "\n *** Example #1: write configuration on system.pdb"
u.atoms.write("system.pdb")

### Example of a selection + analysis done in MDAnalysis
### RDF of all charged particles
from MDAnalysis.analysis.rdf import  InterRDF

charged = u.select_atoms("prop charge  > 0")

rdf = InterRDF(charged, charged,nbins=7,range=(0,10))

rdf.run()
print "\n *** Example #2: rdf of charged particles:", rdf.rdf 


### Let's add more frames...
### Want to use your favourite gromacs analysis tool? No prob, let's save a .trr
print "\n *** Example #3: write gromacs binary trajectory traj.trr"
from MDAnalysis.coordinates.TRR import TRRWriter
W = TRRWriter("traj.trr",n_atoms=len(system.part))
for i in range(10):
	# integrate (nothing, in this example)
	system.integrator.run(1)
	# replace last frame
	u.load_new( eos.stream(), espresso=eos.system())
	# append it to the .trr trajectory
	W.write_next_timestep(u.trajectory.ts)
 

print ""

