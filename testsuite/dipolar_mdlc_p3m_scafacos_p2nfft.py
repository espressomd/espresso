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
import espressomd
import espressomd.magnetostatics as magnetostatics
import espressomd.magnetostatic_extensions as magnetostatic_extensions
import numpy as np
import unittest as ut
from tests_common import abspath


@ut.skipIf(not espressomd.has_features(["DIPOLES", "FFTW"]),
           "Features not available, skipping test!")
class Dipolar_p3m_mdlc_p2nfft(ut.TestCase):
    """Tests mdlc (2d)  as well as dipolar p3m and dipolar p2nfft (3d) against stored data.
       Validity of the stored data:
       2d: as long as this test AND the scafacos_dipolar_1d_2d test passes, we are save.
       3d: As long as the independently written p3m and p2nfft agree, we are save.
    """
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed  = s.cell_system.get_state()['n_nodes'] * [1234]
    s.time_step = 0.01
    s.cell_system.skin = .4
    s.periodicity = 1, 1, 1
    s.thermostat.turn_off()

    def test_mdlc(self):
        s = self.s
        s.part.clear()
        rho = 0.3

        # This is only for box size calculation. The actual particle numbwe is
        # lower, because particles are removed from the mdlc gap region
        n_particle = 100

        particle_radius = 0.5
        box_l = pow(((4 * n_particle * 3.141592654) / (3 * rho)),
                    1.0 / 3.0) * particle_radius
        s.box_l = box_l, box_l, box_l
        f = open(abspath("data/mdlc_reference_data_energy.dat"))
        ref_E = float(f.readline())
        f.close()

        # Particles
        data = np.genfromtxt(abspath("data/mdlc_reference_data_forces_torques.dat"))
        for p in data[:, :]:
            s.part.add(id=int(p[0]), pos=p[1:4], dip=p[4:7])
        s.part[:].rotation=(1,1,1)

        p3m = magnetostatics.DipolarP3M(prefactor=1, mesh=32, accuracy=1E-4)
        dlc = magnetostatic_extensions.DLC(maxPWerror=1E-5, gap_size=2.)
        s.actors.add(p3m)
        s.actors.add(dlc)
        s.thermostat.turn_off()
        s.integrator.run(0)
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

        self.assertLessEqual(abs(err_e), tol_e, "Energy difference too large")
        self.assertLessEqual(abs(err_t), tol_t, "Torque difference too large")
        self.assertLessEqual(abs(err_f), tol_f, "Force difference too large")

        s.part.clear()
        del s.actors[0]
        del s.actors[0]

    def test_p3m(self):
        s = self.s
        s.part.clear()
        rho = 0.09

        # This is only for box size calculation. The actual particle numbwe is
        # lower, because particles are removed from the mdlc gap region
        n_particle = 1000

        particle_radius = 1
        box_l = pow(((4 * n_particle * 3.141592654) / (3 * rho)),
                    1.0 / 3.0) * particle_radius
        s.box_l = box_l, box_l, box_l

        # Particles
        data = np.genfromtxt(abspath("data/p3m_magnetostatics_system.data"))
        for p in data[:, :]:
            s.part.add(id=int(p[0]), pos=p[1:4], dip=p[4:7])
        s.part[:].rotation=(1,1,1)

        p3m = magnetostatics.DipolarP3M(
            prefactor=1, mesh=32, accuracy=1E-6, epsilon="metallic")
        s.actors.add(p3m)
        s.integrator.run(0)
        expected = np.genfromtxt(abspath("data/p3m_magnetostatics_expected.data"))[:, 1:]
        err_f = np.sum(np.sqrt(
            np.sum((s.part[:].f - expected[:, 0:3])**2, 1)), 0) / np.sqrt(data.shape[0])
        err_t = np.sum(np.sqrt(np.sum(
            (s.part[:].torque_lab - expected[:, 3:6])**2, 1)), 0) / np.sqrt(data.shape[0])
        ref_E = 5.570
        err_e = s.analysis.energy()["dipolar"] - ref_E
        print("Energy difference", err_e)
        print("Force difference", err_f)
        print("Torque difference", err_t)

        tol_f = 2E-3
        tol_t = 2E-3
        tol_e = 1E-3

        self.assertLessEqual(abs(err_e), tol_e, "Energy difference too large")
        self.assertLessEqual(abs(err_t), tol_t, "Torque difference too large")
        self.assertLessEqual(abs(err_f), tol_f, "Force difference too large")

        s.part.clear()
        del s.actors[0]

    @ut.skipIf(not espressomd.has_features("SCAFACOS_DIPOLES"),
               "Skipped, because test requires SCAFACOS_DIPOLES")
    def test_scafacos_dipoles(self):
        s = self.s
        s.part.clear()
        rho = 0.09

        # This is only for box size calculation. The actual particle numbwe is
        # lower, because particles are removed from the mdlc gap region
        n_particle = 1000

        particle_radius = 1
        box_l = pow(((4 * n_particle * 3.141592654) / (3 * rho)),
                    1.0 / 3.0) * particle_radius
        s.box_l = box_l, box_l, box_l

        # Particles
        data = np.genfromtxt(abspath("data/p3m_magnetostatics_system.data"))
        for p in data[:, :]:
            s.part.add(id=int(p[0]), pos=p[1:4], dip=p[4:7],rotation=(1,1,1))

        scafacos = magnetostatics.Scafacos(
            prefactor=1,
            method_name="p2nfft",
            method_params={
                "p2nfft_verbose_tuning": 0,
                "pnfft_N": "32,32,32",
                "pnfft_n": "32,32,32",
                "pnfft_window_name": "bspline",
                "pnfft_m": "4",
                "p2nfft_ignore_tolerance": "1",
                "pnfft_diff_ik": "0",
                "p2nfft_r_cut": "11",
                "p2nfft_alpha": "0.31"})
        s.actors.add(scafacos)
        s.integrator.run(0)
        expected = np.genfromtxt(abspath("data/p3m_magnetostatics_expected.data"))[:, 1:]
        err_f = np.sum(np.sqrt(
            np.sum((s.part[:].f - expected[:, 0:3])**2, 1)), 0) / np.sqrt(data.shape[0])
        err_t = np.sum(np.sqrt(np.sum(
            (s.part[:].torque_lab - expected[:, 3:6])**2, 1)), 0) / np.sqrt(data.shape[0])
        ref_E = 5.570
        err_e = s.analysis.energy()["dipolar"] - ref_E
        print("Energy difference", err_e)
        print("Force difference", err_f)
        print("Torque difference", err_t)

        tol_f = 2E-3
        tol_t = 2E-3
        tol_e = 1E-3

        self.assertLessEqual(abs(err_e), tol_e, "Energy difference too large")
        self.assertLessEqual(abs(err_t), tol_t, "Torque difference too large")
        self.assertLessEqual(abs(err_f), tol_f, "Force difference too large")

        s.part.clear()
        del s.actors[0]


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
