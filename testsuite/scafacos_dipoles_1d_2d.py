# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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


# This tests the scafacos p2nfft dipolar calculations by matching against
# reference data from direct summation. In 2d, reference data from the mdlc
# test case is used

from __future__ import print_function
import os
import numpy as np
import unittest as ut
import espressomd
import espressomd.magnetostatics as magnetostatics
import tests_common


@ut.skipIf(not espressomd.has_features(["SCAFACOS_DIPOLES"]),
           "Features not available, skipping test!")
class Scafacos1d2d(ut.TestCase):

    def test_scafacos(self):
        rho = 0.3

        # This is only for box size calculation. The actual particle numbwe is
        # lower, because particles are removed from the mdlc gap region
        n_particle = 100

        particle_radius = 0.5
        dipole_lambda = 3.0

        #################################################

        box_l = pow(((4 * n_particle * 3.141592654) / (3 * rho)),
                    1.0 / 3.0) * particle_radius
        skin = 0.5

        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        # give Espresso some parameters
        s.time_step = 0.01
        s.cell_system.skin = skin
        s.box_l = box_l, box_l, box_l
        for dim in 2, 1:
            print("Dimension", dim)

            # Read reference data
            if dim == 2:
                file_prefix = "data/mdlc"
                s.periodicity = 1, 1, 0
            else:
                s.periodicity = 1, 0, 0
                file_prefix = "data/scafacos_dipoles_1d"

            f = open(tests_common.abspath(
                file_prefix + "_reference_data_energy.dat"))
            ref_E = float(f.readline())
            f.close()

            # Particles
            data = np.genfromtxt(tests_common.abspath(
                file_prefix + "_reference_data_forces_torques.dat"))
            for p in data[:, :]:
                s.part.add(id=int(p[0]), pos=p[1:4], dip=p[4:7],rotation=(1,1,1))

            if dim == 2:
                scafacos = magnetostatics.Scafacos(
                    prefactor=1,
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
                s.actors.add(scafacos)
                # change box geometry in x,y direction to ensure that
                # scafacos survives it
                s.box_l = np.array((1, 1, 1.3)) * box_l

            else:
                if dim == 1:
                    # 1d periodic in x
                    scafacos = magnetostatics.Scafacos(
                        prefactor=1,
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
                    s.box_l = np.array((1, 1, 1)) * box_l
                    s.actors.add(scafacos)
                else:
                    raise Exception("This shouldn't happen.")
            s.thermostat.turn_off()
            s.integrator.run(0)

            # Calculate errors

            err_f = np.sum(np.sqrt(
                np.sum((s.part[:].f - data[:, 7:10])**2, 1)), 0) / np.sqrt(data.shape[0])
            err_t = np.sum(np.sqrt(np.sum(
                (s.part[:].torque_lab - data[:, 10:13])**2, 1)), 0) / np.sqrt(data.shape[0])
            err_e = s.analysis.energy()["dipolar"] - ref_E
            print("Energy difference", err_e)
            print("Force difference", err_f)
            print("Torque difference", err_t)

            tol_f = 2E-3
            tol_t = 2E-3
            tol_e = 1E-3

            self.assertLessEqual(
                abs(err_e), tol_e, "Energy difference too large")
            self.assertLessEqual(
                abs(err_t), tol_t, "Torque difference too large")
            self.assertLessEqual(
                abs(err_f), tol_f, "Force difference too large")

            s.part.clear()
            del s.actors[0]


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
