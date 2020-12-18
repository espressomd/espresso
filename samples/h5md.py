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
Write ESPResSo trajectories in the H5MD format. See :ref:`Writing H5MD-files`.
"""

import espressomd
from espressomd.io.writer import h5md  # pylint: disable=import-error
from espressomd import polymer
from espressomd import interactions

system = espressomd.System(box_l=[100.0, 100.0, 100.0])

system.time_step = 0.01
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.cell_system.skin = 0.4

fene = interactions.FeneBond(k=10, d_r_max=2)
system.bonded_inter.add(fene)

positions = polymer.linear_polymer_positions(n_polymers=5,
                                             beads_per_chain=50,
                                             bond_length=1.0,
                                             seed=1234)
for polymer in positions:
    monomers = system.part.add(pos=polymer)
    previous_part = None
    for part in monomers:
        if previous_part:
            part.add_bond((fene, previous_part))
        previous_part = part

h5_units = h5md.UnitSystem(time='ps', mass='u', length='nm', charge='e')
h5_file = h5md.H5md(file_path="sample.h5", unit_system=h5_units)

for i in range(2):
    h5_file.write()
    system.integrator.run(steps=10)

h5_file.flush()
h5_file.close()
