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
import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np


@utx.skipIfMissingFeatures(["ELECTROSTATICS", "EXTERNAL_FORCES"])
class test_icc(ut.TestCase):
    system = espressomd.System(box_l=[10, 10, 10])

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()

    def add_icc_particles(self, side_num_particles,
                          initial_charge, z_position):
        number = side_num_particles**2
        areas = self.system.box_l[0] * \
            self.system.box_l[1] / number * np.ones(number)
        normals = np.zeros((number, 3))
        normals[:, 2] = 1

        x_position = np.linspace(
            0,
            self.system.box_l[0],
            side_num_particles,
            endpoint=False)
        y_position = np.linspace(
            0,
            self.system.box_l[1],
            side_num_particles,
            endpoint=False)
        x_pos, y_pos = np.meshgrid(x_position, y_position)

        positions = np.stack((x_pos, y_pos, np.full_like(
            x_pos, z_position)), axis=-1).reshape(-1, 3)

        charges = np.full(number, initial_charge)
        fix = [(True, True, True)] * number

        return self.system.part.add(
            pos=positions, q=charges, fix=fix), normals, areas

    def common_setup(self, kwargs, error):
        from espressomd.electrostatic_extensions import ICC

        self.tearDown()
        part_slice, normals, areas = self.add_icc_particles(2, 0.01, 0)

        params = {"n_icc": len(part_slice),
                  "normals": normals,
                  "areas": areas,
                  "epsilons": np.ones_like(areas),
                  "first_id": part_slice.id[0],
                  "check_neutrality": False}

        params.update(kwargs)

        icc = ICC(**params)
        with self.assertRaisesRegex(Exception, error):
            self.system.actors.add(icc)

    def test_params(self):
        params = [({"n_icc": -1}, 'ICC: invalid number of particles'),
                  ({"first_id": -1}, 'ICC: invalid first_id'),
                  ({"max_iterations": -1}, 'ICC: invalid max_iterations'),
                  ({"convergence": -1}, 'ICC: invalid convergence value'),
                  ({"relaxation": -1}, 'ICC: invalid relaxation value'),
                  ({"relaxation": 2.1}, 'ICC: invalid relaxation value'),
                  ({"eps_out": -1}, 'ICC: invalid eps_out'),
                  ({"ext_field": 0}, 'A single value was given but 3 were expected'), ]

        for kwargs, error in params:
            self.common_setup(kwargs, error)

    def test_core_params(self):
        from espressomd.electrostatic_extensions import ICC

        self.tearDown()
        part_slice, normals, areas = self.add_icc_particles(5, 0.01, 0)

        params = {"n_icc": len(part_slice),
                  "normals": normals,
                  "areas": areas,
                  "epsilons": np.ones_like(areas),
                  "first_id": part_slice.id[0],
                  "check_neutrality": False}

        icc = ICC(**params)
        self.system.actors.add(icc)

        icc_params = icc.get_params()
        for key, value in params.items():
            np.testing.assert_allclose(value, np.copy(icc_params[key]))

    @utx.skipIfMissingFeatures(["P3M"])
    def test_dipole_system(self):
        from espressomd.electrostatics import P3M
        from espressomd.electrostatic_extensions import ICC

        BOX_L = 20.
        BOX_SPACE = 5.

        self.tearDown()
        self.system.box_l = [BOX_L, BOX_L, BOX_L + BOX_SPACE]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01

        N_ICC_SIDE_LENGTH = 10
        DIPOLE_DISTANCE = 5.0
        DIPOLE_CHARGE = 10.0

        part_slice_lower, normals_lower, areas_lower = self.add_icc_particles(
            N_ICC_SIDE_LENGTH, -0.0001, 0.)
        part_slice_upper, normals_upper, areas_upper = self.add_icc_particles(
            N_ICC_SIDE_LENGTH, 0.0001, BOX_L)

        assert (part_slice_upper.id[-1] - part_slice_lower.id[0] +
                1) == 2 * N_ICC_SIDE_LENGTH**2, "ICC particles not continuous"

        normals = np.vstack((normals_lower, -normals_upper))
        areas = np.hstack((areas_lower, areas_upper))
        epsilons = np.full_like(areas, 1e8)
        sigmas = np.zeros_like(areas)

        icc = ICC(n_icc=2 * N_ICC_SIDE_LENGTH**2,
                  normals=normals,
                  areas=areas,
                  epsilons=epsilons,
                  sigmas=sigmas,
                  convergence=1e-6,
                  max_iterations=100,
                  first_id=part_slice_lower.id[0],
                  eps_out=1.,
                  relaxation=0.75,
                  ext_field=[0, 0, 0])

        # Dipole in the center of the simulation box
        BOX_L_HALF = BOX_L / 2

        self.system.part.add(pos=[BOX_L_HALF, BOX_L_HALF, BOX_L_HALF - DIPOLE_DISTANCE / 2],
                             q=DIPOLE_CHARGE, fix=[True, True, True])
        self.system.part.add(pos=[BOX_L_HALF, BOX_L_HALF, BOX_L_HALF + DIPOLE_DISTANCE / 2],
                             q=-DIPOLE_CHARGE, fix=[True, True, True])

        p3m = P3M(prefactor=1, mesh=32, cao=7, accuracy=1e-5)

        self.system.actors.add(p3m)
        self.system.actors.add(icc)

        self.system.integrator.run(0)

        charge_lower = sum(part_slice_lower.q)
        charge_upper = sum(part_slice_upper.q)

        testcharge_dipole = DIPOLE_CHARGE * DIPOLE_DISTANCE
        induced_dipole = 0.5 * (abs(charge_lower) + abs(charge_upper)) * BOX_L

        self.assertAlmostEqual(1, induced_dipole / testcharge_dipole, places=4)


if __name__ == "__main__":
    ut.main()
