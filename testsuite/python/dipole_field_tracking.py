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

import espressomd
import espressomd.magnetostatics
import espressomd.observables
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common


def dip_fld_kernel(dipx, dipy, dipz, mx, my, mz, x, y, z):
    dx, dy, dz = x - dipx, y - dipy, z - dipz
    dr = np.linalg.norm((dx, dy, dz))
    dr3 = dr**3
    dr5 = dr**5
    mr = dx * mx + dy * my + dz * mz
    Hx, Hy, Hz = 3. * dx * mr / dr5 - mx / dr3, 3. * dy * \
        mr / dr5 - my / dr3, 3. * dz * mr / dr5 - mz / dr3
    return Hx, Hy, Hz


def N2_loop(particle, slice_prop):
    storage = [0, 0, 0]
    particle_index, particle_x, particle_y, particle_z = particle.id, *particle.pos
    for part_id, part_pos, part_dip in slice_prop:
        if part_id != particle_index:
            dip_x, dip_y, dip_z = part_pos
            m_dip_x, m_dip_y, m_dip_z = part_dip
            dipole_field = dip_fld_kernel(
                dip_x, dip_y, dip_z, m_dip_x, m_dip_y, m_dip_z, particle_x, particle_y, particle_z)
            storage[0] = storage[0] + dipole_field[0]
            storage[1] = storage[1] + dipole_field[1]
            storage[2] = storage[2] + dipole_field[2]
    return storage


@utx.skipIfMissingFeatures(["DIPOLE_FIELD_TRACKING", "WCA"])
class Test(ut.TestCase):
    """
    Check the total dipole field for a magnetic LJ fluid (500 particles,
    density approx 0.002, mu^2=1, no PBC).
    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    data = np.loadtxt(tests_common.data_path("lj_system.dat"))
    pos = data[:, 1:4]

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.integrator.set_vv()

    def setUp(self):
        np.random.seed(42)
        system = self.system
        system.box_l = (50, 50, 50)
        system.cell_system.skin = 0.4
        system.time_step = 0.1
        system.periodicity = [True, True, True]

        n_part = 500
        system.part.add(type=n_part * [0],
                        pos=np.random.random((n_part, 3)) * system.box_l,
                        rotation=n_part * [(True, True, True)])

        # minimize system
        system.non_bonded_inter[0, 0].wca.set_params(epsilon=10., sigma=1.)
        system.integrator.set_steepest_descent(
            f_max=1., gamma=0.001, max_displacement=0.01)
        system.integrator.run(100)
        system.integrator.set_vv()
        system.periodicity = [False, False, False]

        orientor_list = np.random.standard_normal((n_part, 3))
        orientor_list_normalized = np.array(
            [(x / np.linalg.norm(x, axis=0)) for x in orientor_list])
        dip_mom = orientor_list_normalized
        system.part.all().dip = dip_mom

    def test_dds(self):
        solver = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.)
        self.system.actors.add(solver)
        self.system.integrator.run(steps=0)
        self.system.analysis.dipole_fields()
        slice_data = [(x.id, x.pos, x.dip) for x in self.system.part.all()]
        dip_fields_obs = espressomd.observables.ParticleDipoleFields(
            ids=self.system.part.all().id)
        dip_fields = dip_fields_obs.calculate()
        for val, p in zip(dip_fields, self.system.part.all()):
            np.testing.assert_allclose(val, N2_loop(p, slice_data))

        # check auto-update accumulator
        obs = espressomd.observables.ParticleDipoleFields(ids=[0, 1])
        acc = espressomd.accumulators.TimeSeries(obs=obs, delta_N=10)
        self.system.auto_update_accumulators.add(acc)
        self.system.integrator.run(steps=40)
        time_series = acc.time_series()
        rel_diff = 100. * (time_series[-1] - time_series[0]) / time_series[0]
        self.assertGreater(np.linalg.norm(rel_diff), 10)


if __name__ == "__main__":
    ut.main()
