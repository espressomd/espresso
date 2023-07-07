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

import espressomd
import espressomd.electrostatics
import tests_common


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class CoulombCloudWall(ut.TestCase):

    """
    Compare P3M CPU, P3M GPU and ScaFaCoS P2NFFT electrostatic forces
    and energy against stored data.

    """

    system = espressomd.System(box_l=[10., 10., 10.])
    data = np.genfromtxt(tests_common.data_path(
        "coulomb_cloud_wall_system.data"))

    tolerance = 1E-3
    p3m_params = {'r_cut': 1.001, 'accuracy': 1e-3,
                  'mesh': [64, 64, 64], 'cao': 7, 'alpha': 2.70746}

    # Reference energy from P3M
    reference_energy = 148.94229549

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4

        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        self.system.part.add(pos=self.data[:, 1:4], q=self.data[:, 4])
        self.reference_forces = self.data[:, 5:8]

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    def compare(self, method_name, prefactor, force_tol, energy_tol):
        # Compare forces and energy now in the system to reference data
        err_msg = f"difference too large for method {method_name}"

        # Force
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f) / prefactor,
            self.reference_forces, atol=force_tol, err_msg=f"Force {err_msg}")

        # Energy
        self.assertAlmostEqual(
            self.system.analysis.energy()["total"] / prefactor,
            self.reference_energy, delta=energy_tol,
            msg=f"Energy {err_msg}")

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_cpu(self):
        self.system.electrostatics.solver = espressomd.electrostatics.P3M(
            **self.p3m_params, prefactor=3., tune=False)
        self.system.integrator.run(0)
        self.compare("p3m", prefactor=3., force_tol=2e-3, energy_tol=1e-3)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_gpu(self):
        self.system.electrostatics.solver = espressomd.electrostatics.P3MGPU(
            **self.p3m_params, prefactor=2.2, tune=False)
        self.system.integrator.run(0)
        self.compare("p3m_gpu", prefactor=2.2, force_tol=2e-3, energy_tol=1e-3)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p2nfft")
    def test_scafacos_p2nfft(self):
        self.system.electrostatics.solver = espressomd.electrostatics.Scafacos(
            prefactor=2.8,
            method_name="p2nfft",
            method_params={"p2nfft_r_cut": 1.001, "tolerance_field": 1E-4})
        self.system.integrator.run(0)
        self.compare(
            "scafacos_p2nfft",
            prefactor=2.8,
            force_tol=1e-3,
            energy_tol=1e-3)

    def test_zz_deactivation(self):
        # Is the energy and force 0, if no methods active
        self.assertEqual(self.system.analysis.energy()["total"], 0.0)
        self.system.integrator.run(0, recalc_forces=True)
        for p in self.system.part:
            self.assertAlmostEqual(np.linalg.norm(p.f), 0, places=11)


if __name__ == "__main__":
    ut.main()
