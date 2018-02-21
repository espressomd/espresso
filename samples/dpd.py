# A minimal DPD fluid
from __future__ import print_function
import espressomd
import numpy as np

# Set up the box and time step
system = espressomd.System(box_l = 3 * [10])
system.time_step = 0.01
system.cell_system.skin = 0.4

# DPD parameters
n_part = 200
kT=1.
gamma=1.5
r_cut=1.
# Repulsive parameter
F_max=1.

# Activate the thermostat
system.thermostat.set_dpd(kT=kT)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

# Set up the DPD friction interaction
system.non_bonded_inter[0,0].dpd.set_params(
    weight_function=0, gamma=gamma, r_cut=r_cut,
    trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

# Set up the repulsive interaction
system.non_bonded_inter[0,0].hat.set_params(F_max=F_max,
                                       cutoff=r_cut)

# Add the particles randomly distributed over the box
system.part.add(pos=system.box_l * np.random.random((n_part,3)))

# As a usage example, we calculate the pressure at serveral
# particle densities.
for V in range(100, 1000, 100):
    # Rescale the system to the new volume
    system.change_volume_and_rescale_particles(V**0.3333)

    # List of samples
    p_samples = []
    for i in range(100):
        system.integrator.run(1000)
        p_samples.append(system.analysis.pressure()['total'])

    # Average pressure
    p_avg = np.mean(p_samples)
    # And std
    p_std = np.std(p_samples)

    print('rho {:.2f} p {:.2f} ({:.2f})'.format(float(n_part) / V, p_avg, p_std))
