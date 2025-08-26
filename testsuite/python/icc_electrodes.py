#
# Copyright (C) 2023-2025 The ESPResSo project
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
import espressomd.electrostatics
import espressomd.electrostatic_extensions
import numpy as np


@utx.skipIfMissingFeatures(["ELECTROSTATICS", "EXTERNAL_FORCES"])
class Test(ut.TestCase):

    def test_electrodes(self):
        """
        This is a simplified version of the electrodes tutorial part 1.
        Two particles between the electrode plates are moving away from each
        other. Their force is compared against the analytical solution.
        """
        box_l_x = 80.
        box_l_y = 80.
        box_l_z = 5.
        elc_gap = 6. * box_l_z

        system = espressomd.System(box_l=[box_l_x, box_l_y, box_l_z + elc_gap])
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        pos1 = [box_l_x * 0.25, box_l_y / 2., box_l_z / 2.]
        pos2 = [box_l_x * 0.75, box_l_y / 2., box_l_z / 2.]
        p1 = system.part.add(pos=pos1, q=+1., fix=3 * [True])
        p2 = system.part.add(pos=pos2, q=-1., fix=3 * [True])

        xs = np.linspace(0., box_l_x, num=int(box_l_x), endpoint=False)
        ys = np.linspace(0., box_l_y, num=int(box_l_y), endpoint=False)
        surf_q = 1. / (len(xs) * len(ys))
        icc_particles_bot = []
        icc_particles_top = []
        lz = box_l_z
        for x in xs:
            for y in ys:
                icc_particles_bot.append(
                    system.part.add(pos=[x, y, 0.], q=-surf_q, fix=3 * [True]))
                icc_particles_top.append(
                    system.part.add(pos=[x, y, lz], q=+surf_q, fix=3 * [True]))

        n_icc_part = len(icc_particles_top) + len(icc_particles_bot)

        p3m = espressomd.electrostatics.P3M(
            prefactor=2., accuracy=1e-7, mesh=[100, 100, 96], cao=7, check_neutrality=False)
        elc = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=elc_gap, maxPWerror=1e-7, check_neutrality=False)
        system.electrostatics.solver = elc
        icc = espressomd.electrostatic_extensions.ICC(
            first_id=2, n_icc=n_icc_part, convergence=1e-3, eps_out=1.,
            ext_field=[0, 0, 0], max_iterations=100, relaxation=0.95,
            normals=np.tile([[0, 0, 1], [0, 0, -1]], [n_icc_part // 2, 1]),
            areas=np.ones(n_icc_part), sigmas=np.zeros(n_icc_part),
            epsilons=1e5 * np.ones(n_icc_part)
        )
        system.electrostatics.extension = icc

        n_points = 6
        r = np.logspace(0., box_l_z / 5., n_points)
        elc_forces_axial = np.empty((n_points, 2))
        p1.pos = [0., box_l_y / 2., box_l_z / 2.]
        for i in range(n_points):
            p2.pos = [r[i], box_l_y / 2., box_l_z / 2.]
            system.integrator.run(0, recalc_forces=True)
            elc_forces_axial[i, 0] = p1.f[0]
            elc_forces_axial[i, 1] = p2.f[0]
            # reset ICC charges
            for part in icc_particles_bot:
                part.q = +surf_q
            for part in icc_particles_bot:
                part.q = -surf_q

        def calc_sum(x):
            limit = 1000
            accumulator = 0.
            for n in range(-limit + 1, limit + 1):
                accumulator += (-1)**n * x / (x**2 + n**2)**(3. / 2.)
            return accumulator

        ref_force = calc_sum(r / box_l_z) * p3m.prefactor / box_l_z**2
        np.testing.assert_allclose(np.log(elc_forces_axial[:, 0]),
                                   np.log(ref_force), rtol=0.05)
        np.testing.assert_allclose(np.log(-elc_forces_axial[:, 1]),
                                   np.log(ref_force), rtol=0.05)


if __name__ == "__main__":
    ut.main()
