"""
This script compares the diffusion coefficient of a single particle obtained 
from the particle's mean square displacement and the auto correllation function
of its velocity to the expected value.
It uses the Observables/Correlators framework.

"""

from __future__ import division, print_function

gamma=2.4
kT=1.37
dt=0.05

import espressomd
from espressomd.correlators import Correlator 
from espressomd.observables import ParticlePositions, ParticleVelocities
import numpy as np

s=espressomd.System(box_l=[1.0, 1.0, 1.0])
p=s.part.add(pos=(0,0,0),id=0)
s.time_step=dt
s.thermostat.set_langevin(kT=kT,gamma=gamma)
s.cell_system.skin =0.4
s.integrator.run(1000)

pos_obs=ParticlePositions(ids=(0,))
vel_obs=ParticleVelocities(ids=(0,))

c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., dt=10*dt,
    corr_operation="square_distance_componentwise", compress1="discard1")
c_vel = Correlator(obs1=vel_obs, tau_lin=16, tau_max=20., dt=dt,
    corr_operation="scalar_product", compress1="discard1")
s.auto_update_correlators.add(c_pos)
s.auto_update_correlators.add(c_vel)

s.integrator.run(1000000)

c_pos.finalize()
c_vel.finalize()

np.savetxt("msd.dat",c_pos.result())
np.savetxt("vacf.dat",c_vel.result())

# Integral of vacf via Green-Kubo
#D= 1/3 int_0^infty <v(t_0)v(t_0+t)> dt

vacf=c_vel.result()
#Integrate w. trapez rule
I=np.trapz(vacf[:,2],vacf[:,0])
print("Ratio of measured and expected diffusion coeff from Green-Kubo:", 1./3.*I/(kT/gamma))

# Check msd
msd=c_pos.result()

expected_msd = lambda x: 2.*kT/gamma *x

print("Ratio of expected and measured msd")
print("#time ratio_x ratio_y ratio_z")
for i in range(4,msd.shape[0],4):
    print(msd[i,0],msd[i,2:5]/expected_msd(msd[i,0]))
