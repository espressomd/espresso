import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd import analyze
from espressomd import electrostatics
from espressomd import lb
from espressomd import galilei

import sys
import numpy as np


# System parameters
#############################################################
box_l = 48. # size of the simulation box

skin      = 0.3
time_step = 0.01
eq_tstep  = 0.001

n_cycle     = 1000
integ_steps = 150

# Interaction parameters (Lennard Jones for raspberry)
#############################################################
radius_col = 3.
harmonic_radius = 3.0

# the subscript c is for colloid and s is for salt (also used for the surface beads)
c_shift = 0.25 # LJ shift.
eps_ss  = 1.  # LJ epsilon between the colloid's surface particles.
sig_ss  = 1.  # LJ sigma between the colloid's surface particles.
eps_cs  = 48. # LJ epsilon between the colloid's central particle and surface particles. 
sig_cs  = radius_col # LJ sigma between the colloid's central particle and surface particles (colloid's radius). 
a_eff   = 0.32 # effective hydrodynamic radius of a bead due to the discreteness of LB.

# System setup
#############################################################
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)
system.box_l = [box_l, box_l, box_l]
system.cell_system.skin = skin
system.periodicity = [1, 1, 1]

#the LJ potential with the central bead keeps all the beads from simply collapsing into the center
system.non_bonded_inter[1, 0].lennard_jones.set_params(
    epsilon=eps_cs, sigma=sig_cs,
    cutoff=sig_cs*np.power(2.,1./6.), shift=c_shift)
#the LJ potential (WCA potential) between surface beads causes them to be roughly equidistant on the colloid surface
system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*np.power(2.,1./6.), shift=c_shift)

#the harmonic potential pulls surface beads towards the central colloid bead
system.bonded_inter[0] = interactions.HarmonicBond(k=3000., r_0=harmonic_radius)

#for the warmup we use a Langevin thermostat with an extremely low temperature and high friction coefficient such that the trajectories roughly follow 
#the gradient of the potential while not accelerating too much
system.thermostat.set_langevin(kT=0.00001, gamma=40.)

print("# Creating raspberry")
center = box_l/2
colPos = np.ones(3)*center

q_col = 40 # note the charge on the colloid is actually negative
n_col_part = int(4*np.pi*np.power(radius_col,2) + 1) # Number of particles making up the raspberry (surface particles + the central particle).

# Place the central particle 
system.part.add(id=0, pos=colPos, type=0, q=-q_col, fix=(1,1,1),rotation=(1,1,1)) # Create central particle

# this loop create n_col_part surface beads randomly placed on a spherical shell Rvec away from the central bead
# create surface beads uniformly distributed over surface of a sphere with radius=radius_col
for i in range(1,n_col_part):
    colSurfPos=np.random.randn(3)
    colSurfPos=colSurfPos/np.linalg.norm(colSurfPos)*radius_col+colPos
    system.part.add(id=i, pos=colSurfPos, type=1)
    system.part[i].add_bond((0, 0))
print("# Number of colloid beads = {}".format(n_col_part))

#here we let the bead positions relax. The LJ potential with the central bead combined with the
#harmonic bond keep the monomers roughly radius_col away from the central bead. The LJ
#between the surface beads cause them to distribute more or less evenly on the surface.
system.force_cap = 1000
system.time_step=eq_tstep

for i in range(n_cycle):
    print("\r # cycle {} of  {}".format(i+1, n_cycle)),
    sys.stdout.flush()
    system.integrator.run(integ_steps)

system.bonded_inter[0] = interactions.HarmonicBond(k=0, r_0=0)

print("\n# min dist to the center pre moving of surface beads = {}".format(analyze.Analysis(system).min_dist([0],[1])))

#this loop moves the surface beads such that they are once again exactly radius_col away from the center
for i in range(1,n_col_part):
    pos = system.part[i].pos
    system.part[i].pos=(pos-colPos)/np.linalg.norm(pos-colPos)*radius_col+colPos
print("# min dist to the center past moving of surface beads = {}".format(analyze.Analysis(system).min_dist([0],[1])))   

# Select the virtual sites scheme to VirtualSitesRelative
from espressomd.virtual_sites import VirtualSitesRelative
system.virtual_sites= VirtualSitesRelative(have_velocity=True)



#setting min_global_cut is necessary when there is no interaction defined with a range larger than the colloid
#such that the virtual particles are able to communicate their forces to the real particle at the center of the colloid
system.min_global_cut = radius_col+1

#here we calculate the center of mass position (com) and the moment of inertia (momI) of the colloid
com=np.zeros(3)
momI=0
for i in range(n_col_part):
    com+=system.part[i].pos
    momI+=np.power(np.linalg.norm(com-system.part[i].pos),2)
com/=n_col_part

for i in range(n_col_part):
    momI+=np.power(np.linalg.norm(com-system.part[i].pos),2)
#Note the 2/3 factor is for a spherical shell
momI*=(2./3.)

#note that the real particle must be at the center of mass of the colloid because of the integrator
print("\n# moving central particle from {} to {}".format(system.part[0].pos, com))
system.part[0].fix=0,0,0
system.part[0].pos=com
system.part[0].mass=n_col_part
system.part[0].rinertia=np.ones(3)*momI

for i in range(1, n_col_part):
    system.part[i].virtual=1
    system.part[i].vs_auto_relate_to(0)
print("# Raspberry made\n")


print("# Adding the positive ions")
salt_rho = 0.001
volume = np.power(box_l,3)
N_counter_ions = int((volume*salt_rho) + q_col + 0.5)

i=0
while i<N_counter_ions:
    pos=np.random.random(3) * system.box_l
    #make sure the ion is placed outside of the colloid
    if (np.power(np.linalg.norm(pos-center),2) > np.power(radius_col,2)+1):
        system.part.add(id=n_col_part+i, pos=pos, type=2, q=1)
        i+=1

nEndPosIon=n_col_part+N_counter_ions
print("# Added {} positive ions".format(N_counter_ions))

print("\n# Adding the negative ions")

N_co_ions = N_counter_ions - q_col
i=0
while i<N_co_ions:
    pos=np.random.random(3) * system.box_l
    #make sure the ion is placed outside of the colloid
    if (np.power(np.linalg.norm(pos-center),2) > np.power(radius_col,2)+1):
        system.part.add(id=nEndPosIon+i, pos=pos, type=3, q=-1)
        i+=1

print("# Added {} negative ions".format(N_co_ions))
nTot=nEndPosIon+N_co_ions

# Exerting the applied electric force on the particles and checking the charge neutrality of the system
Qtot = 0
Efield=0.1 #an electric field of 0.1 is the upper limit of the linear response regime for this model
for i in range(nTot):
    charge=system.part[i].q
    Qtot+=charge
    system.part[i].ext_force=np.array([charge*Efield, 0., 0.])
#make sure that the system has a net charge of zero, otherwise the system as a whole will slowly accelerate
if (Qtot > 0.0001 or Qtot < -0.0001):
   print("# net charge of {} !!! Exiting".format(Qtot))
   print("# Colloid charge {} Positive ions {} Negative ions {}".format(q_col,N_counter_ions,N_co_ions))
   exit(0)


# WCA interactions for the ions, essentially giving them a finite volume
system.non_bonded_inter[0, 2].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*pow(2.,1./6.), shift=c_shift, offset=sig_cs-1+a_eff)
system.non_bonded_inter[0, 3].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*pow(2.,1./6.), shift=c_shift, offset=sig_cs-1+a_eff)
system.non_bonded_inter[2, 2].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*pow(2.,1./6.), shift=c_shift)
system.non_bonded_inter[2, 3].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*pow(2.,1./6.), shift=c_shift)
system.non_bonded_inter[3, 3].lennard_jones.set_params(
    epsilon=eps_ss, sigma=sig_ss,
    cutoff=sig_ss*pow(2.,1./6.), shift=c_shift)


print("\n# Equilibrating the ions (without electrostatics):")
# Langevin thermostat for warmup before turning on the LB.
system.thermostat.set_langevin(kT=1., gamma=1.)

#fix the colloid for equilibrating the ions
for i in range(n_col_part):
    system.part[i].fix=np.ones(3,dtype=np.int)

ljcap = 100
CapSteps = 1000
for i in range(CapSteps):
    system.force_cap = ljcap
    print("\r # step {} of  {}".format(i+1, CapSteps)),
    sys.stdout.flush()
    system.integrator.run(integ_steps)
    ljcap+=5
 
system.force_cap = 0

#let the colloid move now that bad overlaps have been eliminated
for i in range(n_col_part):
    system.part[i].fix=np.zeros(3,dtype=np.int)

# Turning on the electrostatics
errorCoulomb = 0.01
print("\n# p3m starting...")
sys.stdout.flush()
bjerrum = 2.
p3m = electrostatics.P3M(bjerrum_length=bjerrum, accuracy=0.001)
system.actors.add(p3m)
print("# p3m started!")

system.time_step=time_step

#note that it is important to set all of the particle velocities to zero just before initializing the LB so that the total system momentum is zero
system.galilei.kill_particle_motion()

system.thermostat.turn_off()
lb=espressomd.lb.LBFluidGPU(dens=1., visc=3., agrid=1., tau=time_step, fric=20)
system.actors.add(lb)


#replace the Langevin thermostat with the lb thermostat
system.thermostat.set_lb(kT=1.)

posVsTime = open('posVsTime.dat', 'w')# file where the raspberry position will be written
for i in range(1000):
    system.integrator.run(1000)
    posVsTime.write("{} {}\n".format(system.time,system.part[0].pos))
    print("# time: {}, col_pos: {}".format(system.time,system.part[0].pos))
posVsTime.close()
print("\n# Finished")
