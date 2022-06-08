#
# Copyright (C) 2022 The ESPResSo project
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

"""
Write ESPResSo trajectories in the H5MD format with a Lees-Edwards offset,
aperiodic boundaries and a fluctutating box size. Read the trajectory and
reconstruct the unfolded positions. See :ref:`Writing H5MD-files` for details.
"""

import espressomd
import espressomd.io.writer
import espressomd.lees_edwards

import os
import h5py
import numpy as np
import tempfile

# set particles outside the main box to get folded coordinates
system = espressomd.System(box_l=[6, 7, 8], periodicity=[True, True, False])
system.time_step = 0.1
system.cell_system.skin = 0.0
system.lees_edwards.set_boundary_conditions(
    shear_direction="x", shear_plane_normal="y",
    protocol=espressomd.lees_edwards.Off())
for i in range(12):
    p = system.part.add(pos=3 * [float(i - 4)], v=[0., 1., 0.])

# record particle positions in a H5MD file and in python lists
temp_directory = tempfile.TemporaryDirectory()
temp_file = os.path.join(temp_directory.name, 'test.h5')
h5 = espressomd.io.writer.h5md.H5md(file_path=temp_file)
xyz_folded = []
xyz_unfolded = []
# write coordinates with aperiodic boundary conditions
h5.write()
xyz_folded.append(system.part.all().pos_folded[:])
xyz_unfolded.append(system.part.all().pos[:])
# add Lees-Edwards offset
system.lees_edwards.protocol = espressomd.lees_edwards.LinearShear(
    shear_velocity=1., initial_pos_offset=0., time_0=0.)
system.integrator.run(20)
h5.write()
xyz_folded.append(system.part.all().pos_folded[:])
xyz_unfolded.append(system.part.all().pos[:])
# remove Lees-Edwards offset
system.lees_edwards.protocol = espressomd.lees_edwards.Off()
system.integrator.run(10)
h5.write()
xyz_folded.append(system.part.all().pos_folded[:])
xyz_unfolded.append(system.part.all().pos[:])
# resize box (simulates NpT)
system.box_l = system.box_l + 1.
system.integrator.run(10)
h5.write()
xyz_folded.append(system.part.all().pos_folded[:])
xyz_unfolded.append(system.part.all().pos[:])
h5.flush()
h5.close()
xyz_folded = np.array(xyz_folded)
xyz_unfolded = np.array(xyz_unfolded)

# read H5MD data
with h5py.File(temp_file, mode='r') as h5file:
    p_id = h5file['particles/atoms/id/value'][:]
    argsort = (np.tile(range(p_id.shape[0]), (p_id.shape[1], 1)).T,
               np.unravel_index(np.argsort(p_id, axis=-1), p_id.shape)[1])
    h5_pos = h5file['particles/atoms/position/value'][:][argsort]
    h5_img = h5file['particles/atoms/image/value'][:][argsort]
    h5_box = h5file['particles/atoms/box/edges/value'][:]

# apply unfolding
pos_folded = h5_pos
pos_unfolded = h5_img * h5_box[:, np.newaxis, :] + h5_pos

# report results from last frame
with np.printoptions(precision=1, suppress=True):
    print('unfolded coordinates:')
    print('  from ESPResSo')
    print(xyz_unfolded[-1])
    print('  from H5MD')
    print(pos_unfolded[-1])
    print('')
    print('folded coordinates:')
    print('  from ESPResSo')
    print(xyz_folded[-1])
    print('  from H5MD')
    print(pos_folded[-1])
