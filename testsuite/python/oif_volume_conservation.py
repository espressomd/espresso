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
import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common


@utx.skipIfMissingFeatures("MASS")
class OifVolumeConservation(ut.TestCase):

    """Loads a soft elastic sphere via object_in_fluid, stretches it and checks
       restoration of original volume due to elastic forces."""

    def test(self):
        import object_in_fluid as oif

        system = espressomd.System(box_l=(10, 10, 10))
        self.assertEqual(system.max_oif_objects, 0)
        system.time_step = 0.4
        system.cell_system.skin = 0.5

        # creating the template for OIF object
        cell_type = oif.OifCellType(
            nodes_file=tests_common.abspath("data/sphere393nodes.dat"),
            triangles_file=tests_common.abspath("data/sphere393triangles.dat"),
            system=system, ks=1.0, kb=1.0, kal=1.0, kag=0.1, kv=0.1,
            check_orientation=False, resize=(3.0, 3.0, 3.0))

        # creating the OIF object
        cell0 = oif.OifCell(
            cell_type=cell_type, particle_type=0, origin=[5.0, 5.0, 5.0])
        self.assertEqual(system.max_oif_objects, 1)
        partcls = system.part.all()

        # fluid
        diameter_init = cell0.diameter()
        print(f"initial diameter = {diameter_init}")

        # OIF object is being stretched by factor 1.5
        partcls.pos = (partcls.pos - 5) * 1.5 + 5

        diameter_stretched = cell0.diameter()
        print(f"stretched diameter = {diameter_stretched}")

        # Apply non-isotropic deformation
        partcls.pos = partcls.pos * np.array((0.96, 1.05, 1.02))

        # Test that restoring forces net to zero and don't produce a torque
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.sum(
                partcls.f, axis=0), [
                0., 0., 0.], atol=1E-12)

        total_torque = np.zeros(3)
        for p in system.part:
            total_torque += np.cross(p.pos, p.f)
        np.testing.assert_allclose(total_torque, [0., 0., 0.], atol=2E-12)

        # main integration loop
        system.thermostat.set_langevin(kT=0, gamma=0.7, seed=42)
        # OIF object is let to relax into relaxed shape of the sphere
        for _ in range(2):
            system.integrator.run(steps=240)
            diameter_final = cell0.diameter()
            print(f"final diameter = {diameter_final}")
            self.assertAlmostEqual(
                diameter_final / diameter_init - 1, 0, delta=0.005)


if __name__ == "__main__":
    ut.main()
