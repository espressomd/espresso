#
# Copyright (C) 2010-2020 The ESPResSo project
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
import espressomd.magnetostatic_extensions
import numpy as np
import unittest as ut
import unittest_decorators as utx


@utx.skipIfMissingFeatures(["DIPOLES", "LENNARD_JONES"])
class dds(ut.TestCase):

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.1
    system.periodicity = [1, 1, 1]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.periodicity = [1, 1, 1]

    def dds_gpu_data(self):
        system = self.system

        system.periodicity = [0, 0, 0]
        dds_cpu = espressomd.magnetostatics.DipolarDirectSumGpu(prefactor=1.2)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(system.part[:].f)
        ref_t = np.copy(system.part[:].torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def dds_cpu_data(self):
        system = self.system

        system.periodicity = [0, 0, 0]
        dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.2)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(system.part[:].f)
        ref_t = np.copy(system.part[:].torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def dds_replica_data(self, n_replica):
        system = self.system

        system.periodicity = [1, 1, 1]
        dds_cpu = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1.2, n_replica=n_replica)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(system.part[:].f)
        ref_t = np.copy(system.part[:].torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def fcs_data(self, replicas):
        system = self.system

        system.periodicity = [1, 1, 1]
        scafacos_direct = espressomd.magnetostatics.Scafacos(
            prefactor=1.2,
            method_name="direct",
            method_params={"direct_periodic_images": replicas,
                           "direct_cutoff": 0.})
        system.actors.add(scafacos_direct)

        system.integrator.run(0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(system.part[:].f)
        ref_t = np.copy(system.part[:].torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 1,
               "Skipping test: only runs for n_nodes == 1")
    def test_dds(self):
        system = self.system
        system.box_l = [10., 10., 20.]

        # check for multiple configurations
        for i in range(4):
            np.random.seed(42 + i)

            # add particles
            N = 128
            dipole_modulus = 1.3
            part_pos = np.random.random((N, 3)) * system.box_l
            costheta = 2 * np.random.random(N) - 1
            sintheta = np.sin(np.arcsin(costheta))
            phi = 2 * np.pi * np.random.random(N)
            part_dip = np.array([sintheta * np.cos(phi) * dipole_modulus,
                                 sintheta * np.sin(phi) * dipole_modulus,
                                 costheta * dipole_modulus]).T
            system.part.add(pos=part_pos, dip=part_dip,
                            rotation=N * [(1, 1, 1)])

            # minimize system
            system.periodicity = [1, 1, 1]
            system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=1, cutoff=2**(1.0 / 6.0), shift="auto")
            system.integrator.set_steepest_descent(
                f_max=1, gamma=0.001, max_displacement=0.01)
            system.integrator.run(100)
            system.integrator.set_vv()

            # compute forces and energies for dawaanr
            dds_e, dds_f, dds_t = self.dds_cpu_data()

            # GPU and CPU implementations should agree within float precision
            if espressomd.has_features(
                    "DIPOLAR_DIRECT_SUM") and espressomd.gpu_available():
                gpu_e, gpu_f, gpu_t = self.dds_gpu_data()
                np.testing.assert_allclose(gpu_e, dds_e, atol=1e-5)
                np.testing.assert_allclose(gpu_f, dds_f, atol=1e-4)
                np.testing.assert_allclose(gpu_t, dds_t, atol=1e-4)

            # ESPResSo and ScaFaCoS should agree within double precision
            if espressomd.has_features("SCAFACOS_DIPOLES"):
                fcs_e, fcs_f, fcs_t = self.fcs_data("0,0,0")
                np.testing.assert_allclose(dds_e, fcs_e, rtol=1e-12)
                np.testing.assert_allclose(dds_f, fcs_f, rtol=1e-11)
                np.testing.assert_allclose(dds_t, fcs_t, rtol=1e-11)

            # with replica=0, should be identical to dawaanr
            dds_r0_e, dds_r0_f, dds_r0_t = self.dds_replica_data(0)
            np.testing.assert_allclose(dds_r0_e, dds_e, rtol=1e-12)
            np.testing.assert_allclose(dds_r0_f, dds_f, rtol=1e-11)
            np.testing.assert_allclose(dds_r0_t, dds_t, rtol=1e-11)

            system.part.clear()


if __name__ == "__main__":
    ut.main()
