################################################################################
#                                                                              #
# Copyright (C) 2010-2018 The ESPResSo project                                 #
#                                                                              #
# This file is part of ESPResSo.                                               #
#                                                                              #
# ESPResSo is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# ESPResSo is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.        #
#                                                                              #
################################################################################
#                                                                              #
#                  Active Matter: Swimmer Flow Field Tutorial                  #
#                                                                              #
##########################################################################

from __future__ import print_function

import numpy as np
import os
import sys

import espressomd
from espressomd import assert_features, lb


assert_features(["ENGINE", "LB_GPU", "MASS", "ROTATION", "ROTATIONAL_INERTIA"])

# Read in the hydrodynamic type (pusher/puller) and position

if len(sys.argv) != 3:
    print("Usage:", sys.argv[0], "<type> <pos>")
    exit()

mode = sys.argv[1]
pos = float(sys.argv[2])

##########################################################################

outdir = "./RESULTS_FLOW_FIELD/T_{}_P_{}/".format(mode, pos)
try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

# System parameters

length = 25.0
prod_steps = 1000
prod_length = 50
dt = 0.01

system = espressomd.System(box_l=[length, length, length])
system.cell_system.skin = 0.3
system.time_step = dt
system.min_global_cut = 1.0

##########################################################################

# Set the position of the particle

x0 = 0.5 * length
y0 = 0.5 * length
z0 = 0.5 * length + pos

# Sphere size, mass, and moment of inertia, dipole force

sph_size = 0.5
sph_mass = 4.8
Ixyz = 4.8
force = 0.1

# Setup the particle

system.part.add(
    pos=[x0, y0, z0], type=0, mass=sph_mass, rinertia=[Ixyz, Ixyz, Ixyz],
    swimming={'f_swim': force, 'mode': mode, 'dipole_length': sph_size + 0.5})

##########################################################################

# Setup the fluid (quiescent)

agrid = 1
vskin = 0.1
frict = 20.0
visco = 1.0
densi = 1.0

lbf = lb.LBFluidGPU(agrid=agrid, dens=densi,
                    visc=visco, tau=dt)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=frict, seed=42)

##########################################################################

# Output the coordinates

with open("{}/trajectory.dat".format(outdir), 'w') as outfile:
    print("####################################################", file=outfile)
    print("#        time        position       velocity       #", file=outfile)
    print("####################################################", file=outfile)

    # Production run

    for k in range(prod_steps):
        # Output quantities
        print("{time} {pos[0]} {pos[1]} {pos[2]} {vel[0]} {vel[1]} {vel[2]}"
              .format(time=system.time, pos=system.part[0].pos, vel=system.part[0].v),
              file=outfile)

        # Output 50 simulations
        if k % (prod_steps / 50) == 0:
            num = k // (prod_steps // 50)
            lbf.print_vtk_velocity("{}/lb_velocity_{}.vtk".format(outdir, num))
            system.part.writevtk(
                "{}/position_{}.vtk".format(outdir, num), types=[0])

        system.integrator.run(prod_length)
