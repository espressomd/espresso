from espressomd import lb
from espressomd.observables import ParticlePositions
from espressomd.correlators import Correlator

import numpy as np
import sys

# Constants

loops = 400000
steps = 10
time_step = 0.01
runtime = loops*steps*time_step
box_l = 16
try:
    lb_friction = float(sys.argv[1])
except:
    raise ValueError("Input argument cannot be transformed to float.")

# System setup
system = espressomd.System(box_l=[box_l, box_l, box_l])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
system.time_step = time_step
system.cell_system.skin = 0.4


lbf = lb.LBFluidGPU(agrid=1, dens=1, visc=5, tau=0.01, fric=lb_friction)
system.actors.add(lbf)
system.thermostat.set_lb(kT=1)

system.part.add(pos=[0, 0, 0])


## perform a couple of steps to come to equilbrium
print("Equilibrating the system.")
system.integrator.run(1000)
print("Equlibration finished.")

# Setup observable correlator
pos = ParticlePositions(ids=(0,))
c = Correlator(obs1=pos, tau_lin = 16, tau_max = 1000, dt = time_step,
        corr_operation="square_distance_componentwise", compress1="discard1")
system.auto_update_correlators.add(c)

print("Sampling started.")
for i in range(loops):
    system.integrator.run(steps)

    if i % 1e2 == 0:
        sys.stdout.write("\rSampling: %05i"%i)
        sys.stdout.flush()

print("Sampling finished.")

c.finalize()
corrdata = c.result()
corr = np.zeros((corrdata.shape[0],2))
corr[:,0] = corrdata[:,0]
corr[:,1] = (corrdata[:,2] + corrdata[:,3] + corrdata[:,4]) / 3

np.savetxt("./msd_"+str(lb_friction)+".dat", corr)

