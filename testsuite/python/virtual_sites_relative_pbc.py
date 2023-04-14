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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.virtual_sites
import espressomd.lees_edwards
import numpy as np


@utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE", "ROTATION"])
class Test(ut.TestCase):
    """
    This test case checks the behavior of a virtual and real particle pair
    as they move towards, and eventually cross, a periodic boundary.
    Folded and unfolded coordinates are checked, as well as rotations,
    with and without Lees-Edwards boundary conditions.

    If this test fails, there must be a bug in the virtual site update method,
    the Lees-Edwards update method, or the particle position ghost communicator.
    """
    vs_dist = 1.
    system = espressomd.System(box_l=[20., 20., 20.])
    system.time_step = 0.01
    system.min_global_cut = 1.5 * vs_dist
    system.cell_system.skin = 0.1

    def setUp(self):
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

    def tearDown(self):
        self.system.part.clear()
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesOff()
        self.system.lees_edwards.protocol = None

    def check_pbc(self):
        system = self.system
        vs_dist = self.vs_dist
        box_l = system.box_l[0]
        director = np.array([1, 0, 0])
        unit_cell = np.array([0, 0, 0])
        eps = 1e-4
        v = (vs_dist / 2. + 1.5 * eps) / system.time_step
        le_offset = 0.
        if system.lees_edwards.protocol is not None:
            le_offset = system.lees_edwards.protocol.initial_pos_offset

        # make sure the periodic boundary (x-axis) uses the ghost communicator
        if np.prod(system.cell_system.node_grid) != 1:
            assert system.cell_system.node_grid[0] != 1
        # make sure the periodic boundary (x-axis) is the shear boundary
        if system.lees_edwards.protocol is not None:
            assert system.lees_edwards.shear_plane_normal == "x"
            assert system.lees_edwards.shear_direction == "y"

        def check_dist(real_part, virt_part, check_unfolded_dist=True):
            dist_folded = system.distance(real_part, virt_part)
            self.assertAlmostEqual(dist_folded, vs_dist)
            if check_unfolded_dist:
                dist_unfolded = np.linalg.norm(real_part.pos - virt_part.pos)
                self.assertAlmostEqual(dist_unfolded, vs_dist)

        def check_pos(p, ref_image_box, ref_pos):
            np.testing.assert_allclose(np.copy(p.pos), ref_pos, atol=1e-10)
            np.testing.assert_array_equal(np.copy(p.image_box), ref_image_box)

        for le_dir in (-1., +1.):
            if le_dir == -1:
                # set one particle near to the left periodic boundary and the
                # other one further away by one unit of vs distance
                start = 0.
            elif le_dir == +1:
                # set one particle near to the right periodic boundary and the
                # other one further away by one unit of vs distance
                start = box_l
            le_vec = le_dir * np.array([0., le_offset, 0.])
            image_box = le_dir * director
            # trajectory of the particle nearest to the boundary at each step
            pos_near = np.zeros((3, 3))
            pos_near[0] = [start - le_dir *
                           (eps * 1.0 + vs_dist * 0.0), 1., 1.]
            pos_near[1] = [start + le_dir *
                           (eps * 0.5 + vs_dist * 0.5), 1., 1.]
            pos_near[2] = [start + le_dir *
                           (eps * 2.0 + vs_dist * 1.0), 1. - le_dir * le_offset, 1.]
            pos_away = pos_near - [le_dir * vs_dist, 0., 0.]
            real_kwargs = {"v": le_dir * v * director, "director": director}
            virt_kwargs = {"virtual": True}
            # case 1: virtual site goes through boundary before real particle
            # case 2: virtual site goes through boundary after real particle
            # In the second time step, the virtual site and real particle are
            # in the same image box.
            real_part1 = system.part.add(pos=pos_away[0], **real_kwargs)
            virt_part1 = system.part.add(pos=pos_near[0], **virt_kwargs)
            real_part2 = system.part.add(pos=pos_near[0], **real_kwargs)
            virt_part2 = system.part.add(pos=pos_away[0], **virt_kwargs)
            virt_part1.vs_auto_relate_to(real_part1)
            virt_part2.vs_auto_relate_to(real_part2)
            system.integrator.run(0)
            check_dist(real_part1, virt_part1)
            check_dist(real_part2, virt_part2)
            check_pos(real_part1, unit_cell, pos_away[0])
            check_pos(virt_part1, unit_cell, pos_near[0])
            check_pos(real_part2, unit_cell, pos_near[0])
            check_pos(virt_part2, unit_cell, pos_away[0])
            system.integrator.run(1)
            check_dist(real_part1, virt_part1, check_unfolded_dist=False)
            check_dist(real_part2, virt_part2, check_unfolded_dist=False)
            check_pos(real_part1, unit_cell, pos_away[1] + 0 * le_vec)
            check_pos(virt_part1, image_box, pos_near[1] + 1 * le_vec)
            check_pos(real_part2, image_box, pos_near[1] - 1 * le_vec)
            check_pos(virt_part2, unit_cell, pos_away[1] - 2 * le_vec)
            system.integrator.run(1)
            check_dist(real_part1, virt_part1)
            check_dist(real_part2, virt_part2)
            check_pos(real_part1, image_box, pos_away[2])
            check_pos(virt_part1, image_box, pos_near[2])
            check_pos(real_part2, image_box, pos_near[2])
            check_pos(virt_part2, image_box, pos_away[2])

            system.part.clear()

            # case 1: virtual site goes through boundary before real particle
            # case 2: virtual site goes through boundary after real particle
            # Afterwards, the real particle director changes direction to bring
            # the virtual site in the same image box as the real particle.
            vs_offset = le_dir * director * vs_dist
            real_part1 = system.part.add(pos=pos_away[0], **real_kwargs)
            virt_part1 = system.part.add(pos=pos_near[0], **virt_kwargs)
            real_part2 = system.part.add(pos=pos_near[0], **real_kwargs)
            virt_part2 = system.part.add(pos=pos_away[0], **virt_kwargs)
            virt_part1.vs_auto_relate_to(real_part1)
            virt_part2.vs_auto_relate_to(real_part2)
            system.integrator.run(1)
            # flip director
            real_part1.director = -real_part1.director
            real_part2.director = -real_part2.director
            # all virtual sites rotate by pi around the real particle
            system.integrator.run(0)
            check_dist(real_part1, virt_part1)
            check_dist(real_part2, virt_part2)
            check_pos(real_part1, unit_cell, pos_away[1])
            check_pos(virt_part1, unit_cell, pos_away[1] - vs_offset)
            check_pos(real_part2, image_box, pos_near[1] - le_vec)
            check_pos(virt_part2, image_box, pos_near[1] - le_vec + vs_offset)

            system.part.clear()

    def test_pbc(self):
        system = self.system

        with self.subTest(msg='without Lees-Edwards boundary conditions'):
            self.check_pbc()

        system.part.clear()

        # add Lees-Edwards boundary conditions
        protocol = espressomd.lees_edwards.LinearShear(
            initial_pos_offset=0.1, time_0=0., shear_velocity=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="y", shear_plane_normal="x", protocol=protocol)

        with self.subTest(msg='with Lees-Edwards boundary conditions'):
            self.check_pbc()

    def test_image_box_update(self):
        system = self.system
        # virtual site image box is overriden by the real particle image box
        real_part = system.part.add(pos=[0., 0., +self.vs_dist / 2.])
        virt_part = system.part.add(
            pos=[0., 0., -self.vs_dist / 2. + 4. * system.box_l[2]], virtual=True)
        virt_part.vs_auto_relate_to(real_part)
        np.testing.assert_array_equal(np.copy(real_part.image_box), [0, 0, 0])
        np.testing.assert_array_equal(np.copy(virt_part.image_box), [0, 0, 3])
        system.integrator.run(0)
        np.testing.assert_array_equal(np.copy(real_part.image_box), [0, 0, 0])
        np.testing.assert_array_equal(np.copy(virt_part.image_box), [0, 0, -1])


if __name__ == "__main__":
    ut.main()
