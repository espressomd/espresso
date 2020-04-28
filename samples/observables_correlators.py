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
Measure mean square displacements using the Observables/Correlators framework.
"""

import numpy as np

import espressomd
import espressomd.observables
import espressomd.accumulators

# System setup
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
np.random.seed(seed=42)

system.part.add(pos=(0, 0, 0), v=(1, 2, 3))
system.time_step = 0.01
system.cell_system.skin = 0
system.cell_system.set_n_square(use_verlet_lists=False)
system.thermostat.set_langevin(kT=1, gamma=10, seed=42)
system.integrator.run(1000)

# Initialize observable for a particle with id 0
p = espressomd.observables.ParticlePositions(ids=(0,))
# Ask the observable for its parameters
print(p.get_params())
# Calculate and return current value
print(p.calculate())
# Return stored current value
print(p.calculate())


# Instantiate a correlator correlating the p observable with itself,
# calculating the mean squared displacement (msd).
c = espressomd.accumulators.Correlator(
    tau_lin=16, tau_max=1000, delta_N=1, obs1=p,
    corr_operation="square_distance_componentwise", compress1="discard1")
# Instantiate a correlator calculating the FCS autocorrelation function from
# particle positions, using the symmetric focal spot with wx=wy=wz=10
# (sigma)
fcs = espressomd.accumulators.Correlator(
    tau_lin=16, tau_max=10000, delta_N=10, obs1=p,
    corr_operation="fcs_acf", args=[10, 10, 10], compress1="discard2")
# Ask the correlator for its parameters
print(c.get_params())

# Register the correlator for auto updating at the interval given by its
# dt (currently every timestep)
system.auto_update_accumulators.add(c)
system.auto_update_accumulators.add(fcs)

# Integrate
system.integrator.run(300000)

# Finalize the correlation calculation and write the results to a file
c.finalize()
np.savetxt("res.dat", c.result())
fcs.finalize()
np.savetxt("fcs.dat", fcs.result())
