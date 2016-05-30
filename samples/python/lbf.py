#
# Copyright (C) 2013,2014 The ESPResSo project
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
from __future__ import print_function
import espressomd._system as es
import espressomd
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import lb
import numpy as np

print("""
=======================================================
=                      lbf.py                         =
=======================================================

Program Information:""")
print(code_info.features())


system = espressomd.System()
system.time_step = 0.01
system.skin = 0.1
box_l = 50
system.box_l =[box_l, box_l, box_l]
system.periodic = [1,1,1]

system.part.add(id=0, pos=[box_l/2.0,box_l/2.0,box_l/2.0], fix=[1,1,1])
# system.part.add(id=0, pos=[box_l/2.0,box_l/2.0,box_l/2.0], ext_force=[0,0,1])


lbf = lb.LBFluid(agrid=1, fric=1, dens=1, visc=1, tau=0.01, ext_force=[0,0,-1.0/(box_l**3)])
system.actors.add(lbf)
print(system.actors)


f_list = []
for i in range(10):
    f_list.append(system.part[0].f)
    integrate.integrate(10)
    print(i)

f_list=np.array(f_list)

import matplotlib.pyplot as pp

fig1=pp.figure()
ax=fig1.add_subplot(111)
ax.plot(f_list[:,0],label="F_x")
ax.plot(f_list[:,1],label="F_y")
ax.plot(f_list[:,2],label="F_z")
ax.legend()
ax.set_xlabel("t")
ax.set_ylabel("F")

pp.show()
