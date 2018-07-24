from espressomd import interactions, lb, polymer
from espressomd.observables import ComPosition
from espressomd.correlators import Correlator

from numpy import savetxt, zeros
import sys

# Setup constant
time_step = 0.01
loops = 50000
step_per_loop = 100

# System setup
system = espressomd.System(box_l=[32, 32, 32])
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
system.cell_system.skin = 0.4

try:
    mpc = int(sys.argv[1])
except:
    raise ValueError("First argument cannot be transformed into integer!")
# Lennard-Jones interaction
system.non_bonded_inter[0,0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0, 
    shift=0.25, cutoff=1.226)

# Fene interaction
fene = interactions.FeneBond(k=7, d_r_max=2)
system.bonded_inter.add(fene)


# Setup polymer of part_id 0 with fene bond

polymer.create_polymer(N_P=1, MPC=mpc, bond=fene, bond_length=1)


print("Warming up the polymer chain.")
## For longer chains (>100) an extensive 
## warmup is neccessary ...
system.time_step = 0.002
system.thermostat.set_langevin(kT=1.0, gamma=10)

for i in range(100):
    system.force_cap = i
    system.integrator.run(1000)

print("Warmup finished.")
system.force_cap = 0
system.integrator.run(10000)
system.time_step = time_step
system.integrator.run(50000)

system.thermostat.turn_off()

system.part[:].v = [0,0,0]

lbf = lb.LBFluidGPU(agrid=1, dens=1, visc=5, tau=time_step, fric=5)
system.actors.add(lbf)
system.thermostat.set_lb(kT=1)

print("Warming up the system with LB fluid.")
system.integrator.run(1000)
print("LB fluid warming finished.")


# configure correlators
com_pos = ComPosition(ids=(0,))
c = Correlator(obs1 = com_pos, tau_lin=16, tau_max=1500, dt=time_step,
        corr_operation="square_distance_componentwise", compress1="discard1")
system.auto_update_correlators.add(c)

print("Sampling started.")
for i in range(loops):
    system.integrator.run(step_per_loop)
    system.analysis.append()
    sys.stdout.write("\rSampling: %05i"%i)
    sys.stdout.flush()

c.finalize()
corrdata = c.result()
corr = zeros((corrdata.shape[0],2))
corr[:,0] = corrdata[:,0]
corr[:,1] = (corrdata[:,2] + corrdata[:,3] + corrdata[:,4]) / 3

savetxt("./msd_nom"+str(mpc)+".dat", corr)

with open("./rh_out.dat","a") as datafile:
    rh = system.analysis.calc_rh(chain_start=0, number_of_chains=1, chain_length=mpc-1)
    datafile.write(str(mpc)+ "    " + str(rh[0])+"\n")

