#
# Copyright (C) 2021 The ESPResSo project
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
import espressomd
import espressomd.lb
import espressomd.lees_edwards

import unittest as ut
import unittest_decorators as utx
import numpy as np


system = espressomd.System(box_l=[17, 17, 1])
system.cell_system.skin = 0.1
system.time_step = 0.01


class LBContextManager:
    """
    Add an LB actor and remove it from the actor list at the end.
    """

    def __enter__(self):
        self.lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., viscosity=1., tau=system.time_step)
        system.actors.add(self.lbf)
        system.thermostat.set_lb(LB_fluid=self.lbf, gamma=1.0)
        return self.lbf

    def __exit__(self, exc_type, exc_value, exc_traceback):
        system.actors.remove(self.lbf)
        system.thermostat.turn_off()


class LEContextManager:
    """
    Add a Lees-Edwards linear shear boundary and remove it at the end.
    """

    def __init__(self, shear_direction, shear_plane_normal, offset):
        mapping = {'x': 0, 'y': 1, 'z': 2}
        self.shear_direction = mapping[shear_direction]
        self.shear_plane_normal = mapping[shear_plane_normal]
        self.offset = offset

    def __enter__(self):
        system.lees_edwards.protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=0., initial_pos_offset=self.offset, time_0=0.)
        system.lees_edwards.shear_direction = self.shear_direction
        system.lees_edwards.shear_plane_normal = self.shear_plane_normal

    def __exit__(self, exc_type, exc_value, exc_traceback):
        system.lees_edwards.protocol = espressomd.lees_edwards.Off()


@utx.skipIfMissingFeatures(['LB_WALBERLA'])
class LBLeesEdwards(ut.TestCase):

    """
    Check that velocities interpolated from spatially fixed particles wrap
    around shear boundaries with the correct offset. A two-dimensional LB
    grid is used for simplicity.

    """

    def tearDown(self):
        system.thermostat.turn_off()
        system.actors.clear()
        system.part.clear()
        system.lees_edwards.protocol = espressomd.lees_edwards.Off()

    def sample_lb_velocities(self, lbf):
        profiles = []
        for _ in range(5):
            system.integrator.run(2)
            vel_grid = lbf[:, :, :].velocity[:, :, 0, :]
            profiles.append(np.linalg.norm(vel_grid, axis=2))
        return profiles

    def check_profile(self, profile, stencil,
                      nodes_shifted, nodes_unshifted, tol):
        profile = np.copy(profile) / np.max(profile)
        for node in nodes_unshifted:
            self.assertAlmostEqual(profile[stencil[node]], 1.0, delta=tol)
        for node in nodes_shifted:
            ref = profile[stencil['C']]
            self.assertAlmostEqual(profile[stencil[node]], ref, delta=tol)
            node += '~'
            ref = 1.0 - ref
            self.assertAlmostEqual(profile[stencil[node]], ref, delta=tol)

    def test_velocity_shift_from_particles(self):
        """
        Place particles at the center and borders of a square (a cuboid LB
        grid with thickness 1 in the z-axis). Particles are fixed in space
        and apply a force on the fluid. The velocity vectors of particle
        pairs that are in contact across a periodic boundary are aligned,
        such that their contribution to the interpolation is constructive,
        i.e. at time = 0 the velocity of a LB cell containing a particle
        is 100%, and at time = 10 it is still 100% (70% from the particle,
        30% from the neighboring particle across the periodic boundary).

        Below is a diagram of the time evolution of such a system.
        The magnitude of the interpolated velocity is initially 5,
        and decreases to 4 at the next time step, with a magnitude
        of 1 for neighboring LB cells. LB cells at the boundaries
        remain at velocity = 5 because they gain 1 unit from the
        periodic images.

        .. code-block:: none

            /-------------\        /-------------\
            |      5      |        |     151     |
            |             |        |      1      |
            |             |        | 1    1    1 |
            | 5    5    5 |  --->  | 51  141  15 |
            |             |        | 1    1    1 |
            |             |        |      1      |
            |      5      |        |     151     |
            \-------------/        \-------------/


        When Lees-Edwards boundary conditions are present, contributions
        to the interpolation are no longer constructive across the shear
        boundary due to the shear offset.

        Below is a diagram of the time evolution of such a system,
        where the shear plane normal is the y-axis and the shear
        direction is the x-axis with an offset of 3 agrid:

        .. code-block:: none

            /-------------\        /-------------\
            |      5      |        |     141 1   |
            |             |        |      1      |
            |             |        | 1    1    1 |
            | 5    5    5 |  --->  | 51  141  15 |
            |             |        | 1    1    1 |
            |             |        |      1      |
            |      5      |        |   1 141     |
            \-------------/        \-------------/


        The interpolated velocity at the shear boundary is equal to
        the interpolated velocity of a particle moving diagonally.
        The central particle is moving diagonally and is used as a
        reference.

        """
        tol = 0.012

        # stencil for D2Q8
        stencil_D2Q8 = {'S': (8, 0), 'W': (0, 8), 'N': (8, 16), 'E': (16, 8),
                        'C': (8, 8)}

        # place particles at the square edges and at the center of the square
        for x, y in stencil_D2Q8.values():
            v = np.array([y == 8, x == 8, 0], dtype=float)
            v /= np.linalg.norm(v)
            system.part.add(pos=[x + 0.5, y + 0.5, 0.5], v=v, fix=3 * [True])

        # without Lees-Edwards, velocities remain unaffected
        with LBContextManager() as lbf:
            for profile in self.sample_lb_velocities(lbf):
                self.check_profile(profile, stencil_D2Q8, '', 'SNWE', tol)

        # with Lees-Edwards and no offset, velocities remain unaffected
        with LEContextManager('x', 'y', 0):
            with LBContextManager() as lbf:
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil_D2Q8, '', 'SNWE', tol)

        le_offset = 6

        # North and South are sheared horizontally
        with LEContextManager('x', 'y', le_offset):
            stencil = {'N~': (8 - le_offset, 0),
                       'S~': (8 + le_offset, 16),
                       **stencil_D2Q8}
            with LBContextManager() as lbf:
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil, 'SN', 'WE', tol)

        # East and West are sheared vertically
        with LEContextManager('y', 'x', le_offset):
            stencil = {'E~': (0, 8 - le_offset),
                       'W~': (16, 8 + le_offset),
                       **stencil_D2Q8}
            with LBContextManager() as lbf:
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil, 'WE', 'SN', tol)

    def test_velocity_shift_from_fluid_impulse(self):
        """
        Same test as ``test_velocity_shift_from_particles``, but the particle
        force on the fluid is simulated by manually changing the velocity of
        fluid nodes directly. The velocity is applied one agrid away from the
        shear boundary (at x=1), since velocities stored in the shear boundary
        at x=0 are copied to x=h without any offset.

        """
        tol = 0.08

        # stencil for D2Q8
        stencil_D2Q8 = {'S': (8, 1), 'W': (1, 8), 'N': (8, 15), 'E': (15, 8),
                        'C': (8, 8)}

        def create_impulse(lbf, stencil):
            # add velocities at the square edges and at the center
            for x, y in stencil.values():
                v = np.array([y == 8, x == 8, 0], dtype=float)
                v /= np.linalg.norm(v)
                lbf[x, y, 0].velocity = -0.05 * v

        # without Lees-Edwards, velocities remain unaffected
        with LBContextManager() as lbf:
            create_impulse(lbf, stencil_D2Q8)
            for profile in self.sample_lb_velocities(lbf):
                self.check_profile(profile, stencil_D2Q8, '', 'SNWE', tol)

        # with Lees-Edwards and no offset, velocities remain unaffected
        with LEContextManager('x', 'y', 0):
            with LBContextManager() as lbf:
                create_impulse(lbf, stencil_D2Q8)
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil_D2Q8, '', 'SNWE', tol)

        le_offset = 6

        # North and South are sheared horizontally
        with LEContextManager('x', 'y', le_offset):
            stencil = {'N~': (8 - le_offset, 1),
                       'S~': (8 + le_offset, 15),
                       **stencil_D2Q8}
            with LBContextManager() as lbf:
                create_impulse(lbf, stencil_D2Q8)
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil, 'SN', 'WE', tol)

        # East and West are sheared vertically
        with LEContextManager('y', 'x', le_offset):
            stencil = {'E~': (1, 8 - le_offset),
                       'W~': (15, 8 + le_offset),
                       **stencil_D2Q8}
            with LBContextManager() as lbf:
                create_impulse(lbf, stencil_D2Q8)
                for profile in self.sample_lb_velocities(lbf):
                    self.check_profile(profile, stencil, 'WE', 'SN', tol)


if __name__ == "__main__":
    ut.main()
