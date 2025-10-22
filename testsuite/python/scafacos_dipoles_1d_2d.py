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


# This tests the scafacos p2nfft dipolar calculations by matching against
# reference data from direct summation. In 2d, reference data from the mdlc
# test case is used

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.magnetostatics
import tests_common


@utx.skipIfMissingFeatures(["SCAFACOS_DIPOLES"])
class Scafacos1d2d(ut.TestCase):

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.5
    system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.part.clear()
        self.system.magnetostatics.clear()
        self.system.periodicity = [True, True, True]

    def test_scafacos(self):
        system = self.system
        rho = 0.3

        # This is only for box size calculation. The actual particle number is
        # lower, because particles are removed from the mdlc gap region
        n_particle = 100

        particle_radius = 0.5

        box_l = np.cbrt(4 * n_particle * np.pi / (3 * rho)) * particle_radius

        for dim in (2, 1):
            system.box_l = 3 * [box_l]
            with self.subTest(f"{dim} dimensions"):
                # Read reference data
                if dim == 2:
                    file_prefix = "mdlc"
                    system.periodicity = [True, True, False]
                else:
                    system.periodicity = [True, False, False]
                    file_prefix = "scafacos_dipoles_1d"

                ref_e_path = tests_common.data_path(
                    f"{file_prefix}_reference_data_energy.dat")
                ref_e = float(np.genfromtxt(ref_e_path))

                # Particles
                data = np.genfromtxt(tests_common.data_path(
                    f"{file_prefix}_reference_data_forces_torques.dat"))
                system.part.add(pos=data[:, 1:4], dip=data[:, 4:7])
                system.part.all().rotation = 3 * [True]

                if dim == 2:
                    scafacos = espressomd.magnetostatics.Scafacos(
                        prefactor=1.,
                        method_name="p2nfft",
                        method_params={
                            "p2nfft_verbose_tuning": 0,
                            "pnfft_N": "80,80,160",
                            "pnfft_window_name": "bspline",
                            "pnfft_m": "4",
                            "p2nfft_ignore_tolerance": "1",
                            "pnfft_diff_ik": "0",
                            "p2nfft_r_cut": "6",
                            "p2nfft_alpha": "0.8",
                            "p2nfft_epsB": "0.05"})
                    # change box geometry in x,y direction to ensure that
                    # scafacos survives it
                    system.change_volume_and_rescale_particles(
                        1.3 * box_l, "z")
                    system.part.all().pos = data[:, 1:4]
                else:
                    # 1d periodic in x
                    scafacos = espressomd.magnetostatics.Scafacos(
                        prefactor=1.,
                        method_name="p2nfft",
                        method_params={
                            "p2nfft_verbose_tuning": 1,
                            "pnfft_N": "32,128,128",
                            "pnfft_direct": 0,
                            "p2nfft_r_cut": 2.855,
                            "p2nfft_alpha": "1.5",
                            "p2nfft_intpol_order": "-1",
                            "p2nfft_reg_kernel_name": "ewald",
                            "p2nfft_p": "16",
                            "p2nfft_ignore_tolerance": "1",
                            "pnfft_window_name": "bspline",
                            "pnfft_m": "8",
                            "pnfft_diff_ik": "1",
                            "p2nfft_epsB": "0.125"})
                    system.change_volume_and_rescale_particles(box_l, "z")

                system.magnetostatics.solver = scafacos
                system.integrator.run(0)

                fcs_f = np.copy(system.part.all().f)
                fcs_t = np.copy(system.part.all().torque_lab)
                fcs_e = system.analysis.energy()["dipolar"]
                ref_f = data[:, 7:10]
                ref_t = data[:, 10:13]

                tol_f = 1E-4
                tol_t = 1E-3
                tol_e = 1E-3

                np.testing.assert_allclose(fcs_e, ref_e, atol=tol_e,
                                           err_msg="Energy doesn't match")
                np.testing.assert_allclose(fcs_t, ref_t, atol=tol_t,
                                           err_msg="Torques don't match")
                np.testing.assert_allclose(fcs_f, ref_f, atol=tol_f,
                                           err_msg="Forces don't match")

                system.part.clear()
                system.magnetostatics.solver = None


if __name__ == "__main__":
    ut.main()
