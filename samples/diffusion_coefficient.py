# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Compare the diffusion coefficient of a single thermalized particle obtained
from the particle's mean square displacement and the auto correlation
function of its velocity to the expected value. Uses the
Observables/Correlators framework.
"""

import espressomd
from espressomd.accumulators import Correlator
from espressomd.observables import ParticlePositions, ParticleVelocities
import numpy as np

gamma = 2.4
kT = 1.37
dt = 0.05

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
np.random.seed(seed=42)

p = system.part.add(pos=(0, 0, 0), id=0)
system.time_step = dt
system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
system.cell_system.skin = 0.4
system.integrator.run(1000)

pos_obs = ParticlePositions(ids=(0,))
vel_obs = ParticleVelocities(ids=(0,))

c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., delta_N=10,
                   corr_operation="square_distance_componentwise",
                   compress1="discard1")
c_vel = Correlator(obs1=vel_obs, tau_lin=16, tau_max=20., delta_N=1,
                   corr_operation="scalar_product", compress1="discard1")
system.auto_update_accumulators.add(c_pos)
system.auto_update_accumulators.add(c_vel)

system.integrator.run(1000000)

c_pos.finalize()
c_vel.finalize()

np.savetxt("msd.dat", c_pos.result())
np.savetxt("vacf.dat", c_vel.result())

# Integral of vacf via Green-Kubo
# D= 1/3 int_0^infty <v(t_0)v(t_0+t)> dt

vacf = c_vel.result()
# Integrate with trapezoidal rule
I = np.trapz(vacf[:, 2], vacf[:, 0])
ratio = 1. / 3. * I / (kT / gamma)
print("Ratio of measured and expected diffusion coefficients from Green-Kubo:",
      ratio)

# Check MSD
msd = c_pos.result()


def expected_msd(x):
    return 2. * kT / gamma * x


print("Ratio of expected and measured msd")
print("#time ratio_x ratio_y ratio_z")
for i in range(4, msd.shape[0], 4):
    print(msd[i, 0], msd[i, 2:5] / expected_msd(msd[i, 0]))
