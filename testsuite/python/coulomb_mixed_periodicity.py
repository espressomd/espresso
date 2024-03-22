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


@utx.skipIfMissingFeatures("ELECTROSTATICS")
class CoulombMixedPeriodicity(ut.TestCase):

    """Test mixed periodicity electrostatics"""

    system = espressomd.System(box_l=[10., 10., 10.])
    data = np.genfromtxt(tests_common.data_path(
        "coulomb_mixed_periodicity_system.data"))

    # Reference energy from MMM2D
    ref_energy = 216.640984711

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.

        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        self.system.part.add(pos=self.data[:, 1:4], q=self.data[:, 4])
        self.ref_forces = self.data[:, 5:8]

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    def compare(self, method_name, force_tol, energy_tol):
        self.system.integrator.run(0)
        forces_step1 = np.copy(self.system.part.all().f)
        energy_step1 = self.system.analysis.energy()["total"]

        err_msg = f"difference too large for method {method_name}"
        np.testing.assert_allclose(forces_step1, self.ref_forces, atol=force_tol,
                                   err_msg=f"Force {err_msg}")
        np.testing.assert_allclose(energy_step1, self.ref_energy, atol=energy_tol,
                                   err_msg=f"Energy {err_msg}")

        # triggering a solver re-initialization via a box resize
        # should not affect the forces nor the energies
        original_box_l = np.copy(self.system.box_l)
        self.system.box_l = original_box_l * 1.1
        self.system.box_l = original_box_l
        self.system.integrator.run(0)
        forces_step2 = np.copy(self.system.part.all().f)
        energy_step2 = self.system.analysis.energy()["total"]

        err_msg = f"method {method_name} deviates after cells reinitialization"
        np.testing.assert_allclose(forces_step1, forces_step2, atol=1e-12,
                                   err_msg=f"Force {err_msg}")
        np.testing.assert_allclose(energy_step2, energy_step1, rtol=1e-12,
                                   err_msg=f"Energy {err_msg}")

    def setup_elc_system(self):
        # Make sure, the data satisfies the gap
        for p in self.system.part:
            assert p.pos[2] >= 0. and p.pos[2] <= 9., f'particle {p.id} in gap'

        self.system.cell_system.set_regular_decomposition()
        self.system.cell_system.node_grid = sorted(
            self.system.cell_system.node_grid, key=lambda x: -x)
        self.system.periodicity = [True, True, True]

    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc_cpu(self):
        self.system.box_l = [10., 10., 12.]
        self.setup_elc_system()

        p3m = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-6, mesh=[42, 42, 50], r_cut=3.)
        elc = espressomd.electrostatics.ELC(
            actor=p3m, maxPWerror=1E-6, gap_size=3)

        self.system.electrostatics.solver = elc
        self.compare("elc", force_tol=1e-5, energy_tol=1e-4)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc_gpu(self):
        self.system.box_l = [10., 10., 12.]
        self.setup_elc_system()

        p3m = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-6, mesh=[42, 42, 50], r_cut=3.)
        elc = espressomd.electrostatics.ELC(
            actor=p3m, maxPWerror=1E-6, gap_size=3.)

        self.system.electrostatics.solver = elc
        self.compare("elc", force_tol=1e-5, energy_tol=1e-4)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p2nfft")
    def test_scafacos_p2nfft(self):
        self.system.box_l = [10., 10., 10.]
        self.system.periodicity = [True, True, False]
        self.system.cell_system.set_regular_decomposition()

        scafacos = espressomd.electrostatics.Scafacos(
            prefactor=1,
            method_name="p2nfft",
            method_params={
                "tolerance_field": 5E-5,
                "pnfft_n": "96,96,128",
                "pnfft_N": "96,96,128",
                "r_cut": 2.4,
                "pnfft_m": 3})
        self.system.electrostatics.solver = scafacos
        self.assertTrue(scafacos.call_method("get_near_field_delegation"))
        self.compare("scafacos_p2nfft", force_tol=1e-4, energy_tol=1.8e-3)

        # calculating near field in ScaFaCoS should yield the same result
        scafacos.call_method("set_near_field_delegation", delegate=False)
        self.assertFalse(scafacos.call_method("get_near_field_delegation"))
        self.compare("scafacos_p2nfft", force_tol=1e-4, energy_tol=1.8e-3)


if __name__ == "__main__":
    ut.main()
