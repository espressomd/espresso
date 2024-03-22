#
# Copyright (C) 2010-2022 The ESPResSo project
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
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common
import scipy.optimize
import object_in_fluid as oif


@utx.skipIfMissingFeatures("MASS")
class Test(ut.TestCase):

    system = espressomd.System(box_l=(50., 50., 50.))
    system.time_step = 0.4
    system.cell_system.skin = 0.5

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def make_cell_type(self, **kwargs):
        return oif.OifCellType(
            nodes_file=str(tests_common.data_path("sphere393nodes.dat")),
            triangles_file=str(
                tests_common.data_path("sphere393triangles.dat")),
            system=self.system, **kwargs, kb=1.0, kal=1.0, kag=0.1, kv=0.1,
            check_orientation=False, resize=(3.0, 3.0, 3.0))

    def check_relaxation(self, **kwargs):
        half_box_l = np.copy(self.system.box_l) / 2.
        cell = oif.OifCell(cell_type=self.make_cell_type(**kwargs),
                           particle_type=0, origin=half_box_l)
        partcls = self.system.part.all()

        diameter_init = cell.diameter()
        self.assertAlmostEqual(diameter_init, 24., delta=1e-3)

        # OIF object is being stretched by factor 1.5
        partcls.pos = (partcls.pos - half_box_l) * 1.5 + half_box_l

        diameter_stretched = cell.diameter()
        self.assertAlmostEqual(diameter_stretched, 36., delta=1e-3)

        # Apply non-isotropic deformation
        partcls.pos = partcls.pos * np.array((0.96, 1.05, 1.02))

        # Test that restoring forces net to zero and don't produce a torque
        self.system.integrator.run(1)
        total_force = np.sum(partcls.f, axis=0)
        np.testing.assert_allclose(total_force, [0., 0., 0.], atol=1E-10)
        total_torque = np.zeros(3)
        for p in self.system.part:
            total_torque += np.cross(p.pos, p.f)
        np.testing.assert_allclose(total_torque, [0., 0., 0.], atol=1E-10)

        # Test relaxation to equilibrium diameter (overdamped system)
        self.system.thermostat.set_langevin(kT=0, gamma=0.7, seed=42)
        # warmup
        self.system.integrator.run(steps=50)
        # sampling
        xdata = []
        ydata = []
        for _ in range(20):
            self.system.integrator.run(steps=20)
            xdata.append(self.system.time)
            ydata.append(cell.diameter())
        # check exponential decay
        (prefactor, lam, _, diameter_final), _ = scipy.optimize.curve_fit(
            lambda x, a, b, c, d: a * np.exp(-b * x + c) + d, xdata, ydata,
            p0=[3., 0.03, 0., diameter_init],
            bounds=([-np.inf, 0., -np.inf, 0.], 4 * [np.inf]))
        self.assertGreater(prefactor, 0.)
        self.assertAlmostEqual(diameter_final, diameter_init, delta=0.005)
        self.assertAlmostEqual(lam, 0.0325, delta=0.0001)
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test_00_volume_relaxation(self):
        """
        Load a soft elastic sphere via ``object_in_fluid``, stretch it and
        check restoration of original volume due to elastic forces.
        """
        self.assertEqual(self.system.max_oif_objects, 0)
        with self.subTest(msg="linear stretching"):
            self.check_relaxation(kslin=1.)
        self.assertEqual(self.system.max_oif_objects, 1)
        with self.subTest(msg="neo-Hookean stretching"):
            self.check_relaxation(ks=1.)
        self.assertEqual(self.system.max_oif_objects, 2)

    def test_01_elastic_forces(self):
        half_box_l = np.copy(self.system.box_l) / 2.
        cell = oif.OifCell(cell_type=self.make_cell_type(),
                           particle_type=0, origin=half_box_l)
        # stretch cell
        partcls = self.system.part.all()
        partcls.pos = (partcls.pos - half_box_l) * 1.5 + half_box_l
        # reduce number of triangles to speed up calculation
        cell.mesh.triangles = cell.mesh.triangles[:20]
        # smoke test
        results = cell.elastic_forces(el_forces=6 * [0], f_metric=6 * [1])
        ref = [0., 1.36730815e-12, 22.4985704, 6838.5749, 7.3767594, 6816.6342]
        np.testing.assert_allclose(results, ref, atol=1e-10, rtol=1e-7)
        self.assertEqual(cell.elastic_forces(), 0)
        # check exceptions
        with self.assertRaises(Exception):
            cell.elastic_forces(el_forces=6 * [0], vtk_file="test")
        with self.assertRaises(Exception):
            cell.elastic_forces(el_forces=6 * [0], raw_data_file="test")


if __name__ == "__main__":
    ut.main()
