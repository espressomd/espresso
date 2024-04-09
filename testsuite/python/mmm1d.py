#
# Copyright (C) 2013-2022 The ESPResSo project
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
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common
import espressomd.electrostatics


class ElectrostaticInteractionsTests:

    # Handle to espresso system
    system = espressomd.System(box_l=[10.0] * 3)
    system.periodicity = [0, 0, 1]
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.cell_system.set_n_square()
    system.thermostat.set_langevin(kT=0, gamma=1, seed=8)

    data = np.loadtxt(tests_common.data_path("mmm1d_data.txt"))
    p_pos = data[:, 1:4]
    p_q = data[:, 4]
    forces_target = data[:, 5:8]
    energy_target = -7.156365298205383

    def setUp(self):
        self.system.box_l = [10.0] * 3
        self.system.periodicity = [0, 0, 1]
        self.system.cell_system.set_n_square()

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def test_forces_and_energy(self):
        self.system.part.add(pos=self.p_pos, q=self.p_q)
        mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.system.integrator.run(steps=0)
        measured_f = np.copy(self.system.part.all().f)
        np.testing.assert_allclose(measured_f, self.forces_target,
                                   atol=self.allowed_error)
        measured_el_energy = self.system.analysis.energy()["coulomb"]
        self.assertAlmostEqual(
            measured_el_energy, self.energy_target, delta=self.allowed_error,
            msg="Measured energy deviates too much from stored result")

    def check_with_analytical_result(self, prefactor, accuracy):
        p = self.system.part.by_id(0)
        f_measured = p.f
        energy_measured = self.system.analysis.energy()["total"]
        target_energy_config = -1.00242505606 * prefactor
        target_force_z_config = 0.99510759 * prefactor

        self.assertAlmostEqual(
            f_measured[0], 0, delta=self.allowed_error,
            msg="Measured force in x deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[1], 0, delta=self.allowed_error,
            msg="Measured force in y deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[2], target_force_z_config, delta=accuracy,
            msg="Measured force in z deviates too much from analytical result")
        self.assertAlmostEqual(
            energy_measured, target_energy_config, delta=self.allowed_error,
            msg="Measured energy deviates too much from analytical result")

    def test_with_analytical_result(self):
        self.system.part.add(pos=[0, 0, 0], q=1)
        self.system.part.add(pos=[0, 0, 1], q=-1)
        mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.assertTrue(mmm1d.is_tuned)
        self.system.integrator.run(steps=0, recalc_forces=True)
        self.check_with_analytical_result(prefactor=1.0, accuracy=0.0004)

    def test_bjerrum_length_change(self):
        self.system.part.add(pos=[0, 0, 0], q=1)
        self.system.part.add(pos=[0, 0, 1], q=-1)
        mmm1d = self.MMM1D(prefactor=2.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.assertTrue(mmm1d.is_tuned)
        self.system.integrator.run(steps=0, recalc_forces=True)
        self.check_with_analytical_result(prefactor=2.0, accuracy=0.0017)

        # actor should remain in a valid state after a cell system reset
        forces1 = np.copy(self.system.part.all().f)
        self.system.change_volume_and_rescale_particles(
            self.system.box_l[0], "x")
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid
        self.system.integrator.run(steps=0, recalc_forces=True)
        forces2 = np.copy(self.system.part.all().f)
        np.testing.assert_allclose(forces1, forces2, atol=1e-12, rtol=0.)

    def test_infinite_wire(self):
        """
        For an infinite wire, the energy per ion is :math:`MC\\frac{q}{a}`
        with :math:`M = - 2\\ln{2}` the 1D Madelung constant, :math:`C`
        the electrostatics prefactor, :math:`q` the ion charge and
        :math:`a` the lattice constant. Likewise, the pressure for
        one ion can be derived as :math:`MC\\frac{q}{aV}` with
        :math:`V` the simulation box volume. For more details, see
        Orion Ciftja, "Equivalence of an infinite one-dimensional ionic
        crystal to a simple electrostatic model", Results in Physics,
        Volume 13, 2019, 102325, doi:10.1016/j.rinp.2019.102325
        """
        n_pairs = 128
        n_part = 2 * n_pairs
        self.system.box_l = 3 * [n_part]
        for i in range(n_pairs):
            self.system.part.add(pos=[0., 0., 2. * i + 0.], q=+1.)
            self.system.part.add(pos=[0., 0., 2. * i + 1.], q=-1.)
        self.system.actors.add(self.MMM1D(
            prefactor=1., maxPWerror=1e-20, far_switch_radius=n_pairs / 2.))
        energy = self.system.analysis.energy()["coulomb"]
        p_scalar = self.system.analysis.pressure()["coulomb"]
        p_tensor = self.system.analysis.pressure_tensor()["coulomb"]
        ref_energy = -np.log(2.) * n_part
        np.testing.assert_allclose(energy, ref_energy, atol=0., rtol=5e-7)
        np.testing.assert_allclose(p_scalar, 0., atol=1e-12)
        np.testing.assert_allclose(p_tensor, 0., atol=1e-12)


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class MMM1D_Test(ElectrostaticInteractionsTests, ut.TestCase):

    allowed_error = 2e-5
    MMM1D = espressomd.electrostatics.MMM1D


@utx.skipIfMissingFeatures(["ELECTROSTATICS", "MMM1D_GPU"])
@utx.skipIfMissingGPU()
class MMM1D_GPU_Test(ElectrostaticInteractionsTests, ut.TestCase):

    allowed_error = 1e-4
    MMM1D = espressomd.electrostatics.MMM1DGPU


if __name__ == "__main__":
    ut.main()
