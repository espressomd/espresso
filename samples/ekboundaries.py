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
Set up an electrokinetics (LB) fluid confined between charged walls.
"""

import espressomd

required_features = ["ELECTROKINETICS", "EK_BOUNDARIES", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

from espressomd import System, shapes, electrokinetics, ekboundaries
import os

system = System(box_l=[10, 10, 10])

system.cell_system.skin = 0.4
system.time_step = 0.1


ek = electrokinetics.Electrokinetics(
    lb_density=1, friction=1, agrid=1, viscosity=1, T=1, prefactor=1)

pos = electrokinetics.Species(
    density=0.05, D=0.1, valency=1, ext_force_density=[0, 0, 1.])
neg = electrokinetics.Species(
    density=0.05, D=0.1, valency=-1, ext_force_density=[0, 0, -1.])
ek.add_species(pos)
ek.add_species(neg)
system.actors.add(ek)

print(ek.get_params())
print(pos.get_params())
print(neg.get_params())
print(pos[5, 5, 5].density)


ek_wall_left = ekboundaries.EKBoundary(
    shape=shapes.Wall(dist=1, normal=[1, 0, 0]), charge_density=-0.01)
ek_wall_right = ekboundaries.EKBoundary(
    shape=shapes.Wall(dist=-9, normal=[-1, 0, 0]), charge_density=0.01)
system.ekboundaries.add(ek_wall_left)
system.ekboundaries.add(ek_wall_right)


if not os.path.isdir("ek"):
    os.makedirs("ek")


n_int_cycles = 1000
for i in range(n_int_cycles):
    system.integrator.run(100)
    print("\rIntegrating: %03i" % i, end='', flush=True)

    pos.print_vtk_density("ek/pos_dens_%i.vtk" % i)
    neg.print_vtk_density("ek/neg_dens_%i.vtk" % i)
    pos.print_vtk_flux("ek/pos_flux_%i.vtk" % i)
    neg.print_vtk_flux("ek/neg_flux_%i.vtk" % i)
    ek.print_vtk_velocity("ek/ekv_%i.vtk" % i)
    ek.print_vtk_boundary("ek/ekb_%i.vtk" % i)
