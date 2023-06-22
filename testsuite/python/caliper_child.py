#
# Copyright (C) 2023 The ESPResSo project
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
#
import numpy as np
import espressomd
import espressomd.profiler
import espressomd.electrostatics

np.random.seed(seed=42)

required_features = ["P3M", "WCA", "CALIPER"]
espressomd.assert_features(required_features)

system = espressomd.System(box_l=[10., 10., 10.])
system.time_step = 0.01
system.cell_system.skin = 0.4
system.non_bonded_inter[0, 0].wca.set_params(epsilon=10., sigma=1.)
system.part.add(pos=np.random.random((300, 3)) * system.box_l,
                q=np.repeat((-1., 1.), 150))
p3m = espressomd.electrostatics.P3M(prefactor=1., accuracy=1e-3,
                                    mesh=[16, 16, 16], cao=6)

# warmup
system.integrator.set_steepest_descent(f_max=0., gamma=1e-3,
                                       max_displacement=1. / 100.)
system.integrator.run(100)
system.actors.add(p3m)
system.integrator.run(20)

# production
system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.integrator.run(100)

cali = espressomd.profiler.Caliper()
for i in range(100):
    cali.begin_section(label="calc_energies")
    energies = system.analysis.energy()
    cali.end_section(label="calc_energies")
