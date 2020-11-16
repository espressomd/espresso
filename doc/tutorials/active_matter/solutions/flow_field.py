#
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
#
##########################################################################
#
#            Active Matter: Swimmer Flow Field Tutorial
#
##########################################################################

import os
import sys
import argparse

import espressomd
from espressomd import assert_features, lb


assert_features(["ENGINE", "CUDA", "MASS", "ROTATION", "ROTATIONAL_INERTIA"])


parser = argparse.ArgumentParser()
parser.add_argument('mode',
                    choices=['pusher', 'puller'],
                    help='Type of the swimmer')
parser.add_argument('z_position', type=float,
                    help="Z-position of the swimmer to start off from.")
args = parser.parse_args()

# Read in the hydrodynamic type (pusher/puller) and position
if len(sys.argv) != 3:
    print("Usage:", sys.argv[0], "<type> <pos>")
    exit()

mode = args.mode
pos = args.z_position

##########################################################################

outdir = "./RESULTS_FLOW_FIELD/T_{}_P_{}/".format(mode, pos)
os.makedirs(outdir, exist_ok=True)

# System parameters
LENGTH = 25.0
PROD_STEPS = 1000
PROD_LENGTH = 50
TIME_STEP = 0.01

system = espressomd.System(box_l=[LENGTH, LENGTH, LENGTH])
system.cell_system.skin = 0.3
system.time_step = TIME_STEP
system.min_global_cut = 1.0

##########################################################################

# Set the position of the particle
X0 = 0.5 * LENGTH
Y0 = 0.5 * LENGTH
Z0 = 0.5 * LENGTH + pos

# Sphere size, mass, and moment of inertia, dipole force
SPHERE_SIZE = 0.5
SPHERE_MASS = 4.8
I_XYZ = 4.8
FORCE = 0.1

# Setup the particle
system.part.add(
    pos=[X0, Y0, Z0], type=0, mass=SPHERE_MASS, rinertia=[I_XYZ, I_XYZ, I_XYZ],
    swimming={'f_swim': FORCE, 'mode': mode, 'dipole_length': SPHERE_SIZE + 0.5})

##########################################################################

# Setup the fluid (quiescent)
lbf = lb.LBFluidGPU(agrid=1.0, dens=1.0,
                    visc=1.0, tau=TIME_STEP)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=20.0, seed=42)

##########################################################################

# Output the coordinates
with open("{}/trajectory.dat".format(outdir), 'w') as outfile:
    print("####################################################", file=outfile)
    print("#        time        position       velocity       #", file=outfile)
    print("####################################################", file=outfile)

    # Production run
    for k in range(PROD_STEPS):
        if (k + 1) % 10 == 0:
            print('\rprogress: %.0f%%' % ((k + 1) * 100. / PROD_STEPS), end='')
            sys.stdout.flush()

        # Output quantities
        print("{time} {pos[0]} {pos[1]} {pos[2]} {vel[0]} {vel[1]} {vel[2]}"
              .format(time=system.time, pos=system.part[0].pos, vel=system.part[0].v),
              file=outfile)

        # Output 50 simulations
        if k % (PROD_STEPS / 50) == 0:
            num = k // (PROD_STEPS // 50)
            lbf.print_vtk_velocity("{}/lb_velocity_{}.vtk".format(outdir, num))
            system.part.writevtk(
                "{}/position_{}.vtk".format(outdir, num), types=[0])

        system.integrator.run(PROD_LENGTH)
print()
