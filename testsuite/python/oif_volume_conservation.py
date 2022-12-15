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
class OifVolumeConservation(ut.TestCase):

    """Loads a soft elastic sphere via object_in_fluid, stretches it and checks
       restoration of original volume due to elastic forces."""

    system = espressomd.System(box_l=(50., 50., 50.))
    system.time_step = 0.4
    system.cell_system.skin = 0.5

    def check_relaxation(self, **kwargs):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        half_box_l = np.copy(self.system.box_l) / 2.

        # creating the template for OIF object
        cell_type = oif.OifCellType(
            nodes_file=str(tests_common.data_path("sphere393nodes.dat")),
            triangles_file=str(
                tests_common.data_path("sphere393triangles.dat")),
            system=self.system, **kwargs, kb=1.0, kal=1.0, kag=0.1, kv=0.1,
            check_orientation=False, resize=(3.0, 3.0, 3.0))

        # creating the OIF object
        cell0 = oif.OifCell(
            cell_type=cell_type, particle_type=0, origin=half_box_l)
        partcls = self.system.part.all()

        diameter_init = cell0.diameter()
        self.assertAlmostEqual(diameter_init, 24., delta=1e-3)

        # OIF object is being stretched by factor 1.5
        partcls.pos = (partcls.pos - half_box_l) * 1.5 + half_box_l

        diameter_stretched = cell0.diameter()
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
            ydata.append(cell0.diameter())
        # check exponential decay
        (prefactor, lam, _, diameter_final), _ = scipy.optimize.curve_fit(
            lambda x, a, b, c, d: a * np.exp(-b * x + c) + d, xdata, ydata,
            p0=[3., 0.03, 0., diameter_init],
            bounds=([-np.inf, 0., -np.inf, 0.], 4 * [np.inf]))
        self.assertGreater(prefactor, 0.)
        self.assertAlmostEqual(diameter_final, diameter_init, delta=0.005)
        self.assertAlmostEqual(lam, 0.0325, delta=0.0001)

    def test(self):
        self.assertEqual(self.system.max_oif_objects, 0)
        with self.subTest(msg='linear stretching'):
            self.check_relaxation(kslin=1.)
        self.assertEqual(self.system.max_oif_objects, 1)
        with self.subTest(msg='neo-Hookean stretching'):
            self.check_relaxation(ks=1.)
        self.assertEqual(self.system.max_oif_objects, 2)


if __name__ == "__main__":
    ut.main()
