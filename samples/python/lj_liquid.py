from __future__ import print_function
import espressomd.System as es
import numpy
from espressomd import code_info

print("""
=======================================================
=                    lj_liquid.py                     =
=======================================================

Program Information:""")
print(code_info.features())

dev="cpu"

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps   =  1.0
lj_sig   =  1.0
lj_cut   =  1.12246
lj_cap   = 20 

# Integration parameters
#############################################################
es.glob.time_step = 0.01
es.glob.skin      = 0.4
#es._espressoHandle.Tcl_Eval('thermostat langevin 1.0 1.0')
es.thermostat.Thermostat().setLangevin(1.0,1.0)

# warmup integration (with capped LJ potential)
warm_steps   = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist     = 0.9

# integration
int_steps   = 1000
int_n_times = 5


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

es.glob.box_l = [box_l,box_l,box_l]

es.nonBondedInter[0,0].lennardJones.setParams(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
es.nonBondedInter.setForceCap(lj_cap)

print("LJ-parameters:")
print(es.nonBondedInter[0,0].lennardJones.getParams())

# Particle setup
#############################################################

volume = box_l*box_l*box_l
n_part = int(volume*density)

for i in range(n_part):
  es.part[i].pos=numpy.random.random(3)*es.glob.box_l

es.analyze.distto(0)

print("Simulate {} particles in a cubic simulation box {} at density {}."
  .format(n_part, box_l, density).strip())
#print "Interactions:\n"	# Nicht angepasst
#act_min_dist = float(es._espressoHandle.Tcl_Eval('analyze mindist'))
act_min_dist = es.analyze.mindist()
print("Start with minimal distance {}".format(act_min_dist))

es.glob.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

#open Observable file
obs_file = open("pylj_liquid.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\n")
#set obs_file [open "$name$ident.obs" "w"]
#puts $obs_file "\# System: $name$ident"
#puts $obs_file "\# Time\tE_tot\tE_kin\t..."

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
es.nonBondedInter.setForceCap(lj_cap)
print(es.nonBondedInter[0,0].lennardJones)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
  es.integrate(warm_steps)
  # Warmup criterion
#  act_min_dist = float(es._espressoHandle.Tcl_Eval('analyze mindist'))
  act_min_dist = es.analyze.mindist() 
#  print("\rrun %d at time=%f (LJ cap=%f) min dist = %f\r" % (i,es.glob.time,lj_cap,act_min_dist), end=' ')
  i += 1

#   write observables
#    puts $obs_file "{ time [setmd time] } [analyze energy]"

#   Increase LJ cap
  lj_cap = lj_cap + 10
  es.nonBondedInter.setForceCap(lj_cap)

# Just to see what else we may get from the c code
print("""
ro variables:
cell_grid     {0.cell_grid}
cell_size     {0.cell_size} 
local_box_l   {0.local_box_l} 
max_cut       {0.max_cut}
max_part      {0.max_part}
max_range     {0.max_range} 
max_skin      {0.max_skin}
n_nodes       {0.n_nodes}
n_part        {0.n_part}
n_part_types  {0.n_part_types}
periodicity   {0.periodicity}
transfer_rate {0.transfer_rate}
verlet_reuse  {0.verlet_reuse}
""".format(es.glob))

# write parameter file

#polyBlockWrite "$name$ident.set" {box_l time_step skin} "" 
set_file = open("pylj_liquid.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" % (box_l, es.glob.time_step, es.glob.skin))

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0 
es.nonBondedInter.setForceCap(lj_cap)
print(es.nonBondedInter[0,0].lennardJones)

# print initial energies
#energies = es._espressoHandle.Tcl_Eval('analyze energy')
energies = es.analyze.energy()
print(energies)

j = 0
for i in range(0,int_n_times):
  print("run %d at time=%f " % (i,es.glob.time))

#  es._espressoHandle.Tcl_Eval('integrate %d' % int_steps)
  es.integrate(int_steps)
  
#  energies = es._espressoHandle.Tcl_Eval('analyze energy')
  energies = es.analyze.energy()
  print(energies)
  obs_file.write('{ time %s } %s\n' % (es.glob.time,energies))

#   write observables
#    set energies [analyze energy]
#    puts $obs_file "{ time [setmd time] } $energies"
#    puts -nonewline "temp = [expr [lindex $energies 1 1]/(([degrees_of_freedom]/2.0)*[setmd n_part])]\r"
#    flush stdout

#   write intermediate configuration
#    if { $i%10==0 } {
#	polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
#	incr j
#    }

# write end configuration
end_file = open("pylj_liquid.end", "w")
end_file.write("{ time %f } \n { box_l %f }\n" % (es.glob.time, box_l) )
end_file.write("{ particles {id pos type} }")
for i in range(n_part):
	end_file.write("%s\n" % es.part[i].pos)
	# id & type not working yet

obs_file.close()
set_file.close()
end_file.close()
#es._espressoHandle.die()

# terminate program
print("\nFinished.")
