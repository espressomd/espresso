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
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd.magnetostatics


@utx.skipIfMissingFeatures(["DP3M"])
class MagnetostaticsP3M(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])
    system.cell_system.skin = 0.

    def tearDown(self):
        self.system.part.clear()
        self.system.magnetostatics.clear()

    def test_dp3m(self):
        self.system.time_step = 0.01
        prefactor = 1.1
        box_vol = self.system.volume()
        self.system.part.add(pos=[1., 2., 3.])  # non-magnetic particle
        partcls = self.system.part.add(
            pos=[(4.0, 2.0, 2.0), (6.0, 2.0, 2.0)],
            dip=[(1.3, 2.1, -6.0), (7.3, 6.1, -4.0)])
        p1, p2 = partcls
        dip = np.copy(p1.dip + p2.dip)
        dp3m_params = {'accuracy': 1e-6,
                       'mesh': [49, 49, 49],
                       'cao': 7,
                       'r_cut': 4.739799499511719,
                       'alpha': 0.9056147262573242}
        mdlc_params = {'maxPWerror': 1e-5, 'gap_size': 5.}

        # reference values for energy and force calculated for prefactor = 1.1
        ref_dp3m_energy = 1.673333
        ref_dp3m_force = np.array([-3.54175042, -4.6761059, 9.96632774])
        ref_dp3m_torque1 = np.array([-3.29316117, -13.21245739, -5.33787892])
        ref_dp3m_torque2 = np.array([3.98103932, -7.47123148, -4.12823244])

        # check metallic case
        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=prefactor, epsilon='metallic', tune=False, **dp3m_params)
        self.system.magnetostatics.solver = dp3m
        self.system.integrator.run(0, recalc_forces=True)
        energy = self.system.analysis.energy()['dipolar']
        tol = 1e-5
        np.testing.assert_allclose(energy, ref_dp3m_energy, atol=tol)
        np.testing.assert_allclose(np.copy(p1.f), ref_dp3m_force, atol=tol)
        np.testing.assert_allclose(np.copy(p2.f), -ref_dp3m_force, atol=tol)
        np.testing.assert_allclose(
            np.copy(p1.convert_vector_space_to_body(p1.torque_lab)),
            ref_dp3m_torque1, atol=tol)
        np.testing.assert_allclose(
            np.copy(p2.convert_vector_space_to_body(p2.torque_lab)),
            ref_dp3m_torque2, atol=tol)

        # keep current values as reference to check for DP3M dipole correction
        ref_dp3m_energy_metallic = self.system.analysis.energy()['dipolar']
        ref_dp3m_forces_metallic = np.copy(partcls.f)
        ref_dp3m_torque_metallic = np.array([
            p1.convert_vector_space_to_body(p1.torque_lab),
            p2.convert_vector_space_to_body(p2.torque_lab)])

        # MDLC cancels out dipole correction
        self.system.magnetostatics.clear()
        mdlc = espressomd.magnetostatics.DLC(actor=dp3m, **mdlc_params)
        self.system.magnetostatics.solver = mdlc
        self.assertAlmostEqual(mdlc.prefactor, 1.1, delta=1e-12)

        # keep current values as reference to check for MDLC dipole correction
        self.system.integrator.run(0, recalc_forces=True)
        ref_mdlc_energy_metallic = self.system.analysis.energy()['dipolar']
        ref_mdlc_forces_metallic = np.copy(partcls.f)
        ref_mdlc_torque_metallic = np.copy(partcls.torque_lab)

        # solvers should remain in a valid state after a cell system reset
        self.system.change_volume_and_rescale_particles(
            self.system.box_l[0], "x")
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid
        self.system.integrator.run(0, recalc_forces=True)
        mdlc_forces = np.copy(partcls.f)
        mdlc_torque = np.copy(partcls.torque_lab)
        mdlc_energy = self.system.analysis.energy()['dipolar']
        np.testing.assert_allclose(mdlc_forces, ref_mdlc_forces_metallic,
                                   atol=1e-12)
        np.testing.assert_allclose(mdlc_torque, ref_mdlc_torque_metallic,
                                   atol=1e-12)
        np.testing.assert_allclose(mdlc_energy, ref_mdlc_energy_metallic,
                                   atol=1e-12)

        self.system.magnetostatics.clear()

        # check non-metallic case
        tol = 1e-10
        for epsilon in np.power(10., np.arange(-4, 5, 2)):
            dipole_correction = 4. * np.pi / box_vol / (1. + 2. * epsilon)
            e_correction = dipole_correction / 2. * np.linalg.norm(dip)**2
            t_correction = np.cross([p1.dip, p2.dip], dipole_correction * dip)
            ref_dp3m_energy = ref_dp3m_energy_metallic + prefactor * e_correction
            ref_dp3m_forces = ref_dp3m_forces_metallic
            ref_dp3m_torque = ref_dp3m_torque_metallic - prefactor * t_correction
            dp3m = espressomd.magnetostatics.DipolarP3M(
                prefactor=prefactor, epsilon=epsilon, tune=False, **dp3m_params)
            self.system.magnetostatics.solver = dp3m
            self.system.integrator.run(0, recalc_forces=True)
            dp3m_forces = np.copy(partcls.f)
            dp3m_torque = np.array([
                p1.convert_vector_space_to_body(p1.torque_lab),
                p2.convert_vector_space_to_body(p2.torque_lab)])
            dp3m_energy = self.system.analysis.energy()['dipolar']
            np.testing.assert_allclose(dp3m_forces, ref_dp3m_forces, atol=tol)
            np.testing.assert_allclose(dp3m_torque, ref_dp3m_torque, atol=tol)
            np.testing.assert_allclose(dp3m_energy, ref_dp3m_energy, atol=tol)

            # MDLC cancels out dipole correction
            ref_mdlc_energy = ref_mdlc_energy_metallic
            ref_mdlc_forces = ref_mdlc_forces_metallic
            ref_mdlc_torque = ref_mdlc_torque_metallic
            self.system.magnetostatics.clear()
            mdlc = espressomd.magnetostatics.DLC(actor=dp3m, **mdlc_params)
            self.system.magnetostatics.solver = mdlc
            self.system.integrator.run(0, recalc_forces=True)
            mdlc_forces = np.copy(partcls.f)
            mdlc_torque = np.copy(partcls.torque_lab)
            mdlc_energy = self.system.analysis.energy()['dipolar']
            np.testing.assert_allclose(mdlc_forces, ref_mdlc_forces, atol=tol)
            np.testing.assert_allclose(mdlc_torque, ref_mdlc_torque, atol=tol)
            np.testing.assert_allclose(mdlc_energy, ref_mdlc_energy, atol=tol)

            self.system.magnetostatics.clear()


if __name__ == "__main__":
    ut.main()
