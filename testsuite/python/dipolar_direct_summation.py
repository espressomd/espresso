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
import os
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common
OPEN_BOUNDARIES_REF_ENERGY = tests_common.abspath(
    "data/dipolar_open_boundaries_energy.npy")
OPEN_BOUNDARIES_REF_ARRAYS = tests_common.abspath(
    "data/dipolar_open_boundaries_arrays.npy")


@utx.skipIfMissingFeatures(["DIPOLES"])
class dds(ut.TestCase):

    system = espressomd.System(box_l=[3, 3, 3])

    system.time_step = 0.01
    system.cell_system.skin = 0.1
    system.periodicity = [False, False, False]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def dds_gpu_data(self):
        system = self.system

        dds_cpu = espressomd.magnetostatics.DipolarDirectSumGpu(prefactor=1.2)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(self.particles.f)
        ref_t = np.copy(self.particles.torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def dds_data(self):
        system = self.system

        dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.2)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(self.particles.f)
        ref_t = np.copy(self.particles.torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def dds_replica_data(self):
        system = self.system

        dds_cpu = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1.2, n_replica=0)
        system.actors.add(dds_cpu)

        system.integrator.run(steps=0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(self.particles.f)
        ref_t = np.copy(self.particles.torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    def fcs_data(self):
        system = self.system

        scafacos_direct = espressomd.magnetostatics.Scafacos(
            prefactor=1.2,
            method_name="direct",
            method_params={"direct_periodic_images": "0,0,0",
                           "direct_cutoff": 0})
        system.actors.add(scafacos_direct)

        system.integrator.run(0, recalc_forces=True)
        ref_e = system.analysis.energy()["dipolar"]
        ref_f = np.copy(self.particles.f)
        ref_t = np.copy(self.particles.torque_lab)

        system.actors.clear()

        return (ref_e, ref_f, ref_t)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 1,
               "Skipping test: only runs for n_nodes == 1")
    def test_gen_reference_data(self):
        filepaths = ('dipolar_direct_summation_energy.npy',
                     'dipolar_direct_summation_arrays.npy')
        for filepath in filepaths:
            if os.path.isfile(filepath):
                os.remove(filepath)

        self.gen_reference_data(filepaths[0], filepaths[1])
        for filepath in filepaths:
            self.assertTrue(os.path.isfile(filepath))

    def gen_reference_data(self, filepath_energy=OPEN_BOUNDARIES_REF_ENERGY,
                           filepath_arrays=OPEN_BOUNDARIES_REF_ARRAYS):
        system = self.system
        np.random.seed(42)

        # add particles
        N = 20
        dipole_modulus = 1.3
        part_pos = np.random.random((N, 3)) * system.box_l
        part_dip = dipole_modulus * tests_common.random_dipoles(N)
        self.particles = system.part.add(pos=part_pos, dip=part_dip,
                                         rotation=N * [(1, 1, 1)])

        # minimize system
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=10.0, sigma=1, cutoff=2**(1.0 / 6.0), shift="auto")
        system.integrator.set_steepest_descent(
            f_max=1, gamma=0.001, max_displacement=0.01)
        system.integrator.run(100)
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0.0, sigma=0, cutoff=0, shift=0)
        system.integrator.set_vv()
        assert system.analysis.energy()["total"] == 0

        # compute forces and energies for dawaanr
        ref_e, ref_f, ref_t = self.dds_data()
        np.save(
            filepath_energy,
            np.array([ref_e]),
            allow_pickle=False)
        np.save(
            filepath_arrays,
            np.hstack(
                (self.particles.pos_folded,
                 self.particles.dip,
                 ref_f,
                 ref_t)),
            allow_pickle=False)

    def check_open_bc(self, method_func=None, energy_tol=None,
                      force_tol=None, torque_tol=None):
        system = self.system

        ref_e = np.load(OPEN_BOUNDARIES_REF_ENERGY)[0]
        array_data = np.load(OPEN_BOUNDARIES_REF_ARRAYS)
        pos = array_data[:, :3]
        dip = array_data[:, 3:6]
        ref_f = array_data[:, 6:9]
        ref_t = array_data[:, 9:12]

        self.particles = system.part.add(
            pos=pos, dip=dip, rotation=[[1, 1, 1]] * len(pos))
        dds_e, dds_f, dds_t = method_func()
        self.assertAlmostEqual(dds_e, ref_e, delta=energy_tol)
        np.testing.assert_allclose(dds_f, ref_f, atol=force_tol)
        np.testing.assert_allclose(dds_t, ref_t, atol=torque_tol)

    def test_dds_cpu(self):
        self.check_open_bc(
            self.dds_data,
            energy_tol=1E-12,
            force_tol=1E-12,
            torque_tol=1E-12)

    def test_dds_cpu_replica_data(self):
        self.check_open_bc(
            self.dds_replica_data,
            energy_tol=1E-12,
            force_tol=1E-12,
            torque_tol=1E-12)

    @utx.skipIfMissingFeatures("DIPOLAR_DIRECT_SUM")
    @utx.skipIfMissingGPU()
    def test_dds_gpu(self):
        self.check_open_bc(
            self.dds_gpu_data,
            energy_tol=1E-5,
            force_tol=1E-4,
            torque_tol=1E-4)

    @ut.skipIf(not espressomd.has_features("SCAFACOS_DIPOLES") or
               "direct" not in espressomd.scafacos.available_methods(),
               "Skipping test: missing SCAFACOS_DIPOLES or 'direct' method")
    def test_dds_scafacos(self):
        self.check_open_bc(
            self.fcs_data,
            energy_tol=1E-12,
            force_tol=1E-12,
            torque_tol=1E-12)


if __name__ == "__main__":
    ut.main()
