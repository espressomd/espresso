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
import espressomd
import espressomd.magnetostatics
import pathlib
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common
OPEN_BOUNDARIES_REF_ENERGY = tests_common.data_path(
    "dipolar_open_boundaries_energy.npy")
OPEN_BOUNDARIES_REF_ARRAYS = tests_common.data_path(
    "dipolar_open_boundaries_arrays.npy")


@utx.skipIfMissingFeatures(["DIPOLES"])
class Test(ut.TestCase):

    system = espressomd.System(box_l=[3, 3, 3])

    system.time_step = 0.01
    system.cell_system.skin = 0.1
    system.periodicity = [False, False, False]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.periodicity = [False, False, False]

    def dds_gpu_data(self):
        system = self.system

        dds_cpu = espressomd.magnetostatics.DipolarDirectSumGpu(prefactor=1.2)
        system.actors.add(dds_cpu)
        # check MD cell reset has no impact
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid

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
        # check MD cell reset has no impact
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid

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

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_gen_reference_data(self):
        filepaths = ('dipolar_direct_summation_energy.npy',
                     'dipolar_direct_summation_arrays.npy')
        for filepath in filepaths:
            pathlib.Path(filepath).unlink(missing_ok=True)

        self.gen_reference_data(filepaths[0], filepaths[1])
        for filepath in filepaths:
            self.assertTrue(pathlib.Path(filepath).is_file())

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
                                         rotation=N * [(True, True, True)])

        # minimize system
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=10.0, sigma=1, cutoff=2**(1.0 / 6.0), shift="auto")
        system.integrator.set_steepest_descent(
            f_max=1, gamma=0.001, max_displacement=0.01)
        system.integrator.run(100)
        system.non_bonded_inter[0, 0].lennard_jones.deactivate()
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
            pos=pos, dip=dip, rotation=[[True, True, True]] * len(pos))
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

    @utx.skipIfMissingFeatures("DIPOLAR_DIRECT_SUM")
    @utx.skipIfMissingGPU()
    def test_dds_gpu(self):
        self.check_open_bc(
            self.dds_gpu_data,
            energy_tol=1E-5,
            force_tol=1E-4,
            torque_tol=1E-4)

    @utx.skipIfMissingFeatures(["SCAFACOS_DIPOLES"])
    @utx.skipIfMissingScafacosMethod("direct")
    def test_dds_scafacos(self):
        self.check_open_bc(
            self.fcs_data,
            energy_tol=1E-12,
            force_tol=1E-12,
            torque_tol=1E-12)

    def check_min_image_convention(self, solver, rtol):
        system = self.system
        tol = {"atol": 1e-10, "rtol": rtol}
        # place particles close to each other in space
        system.part.add(pos=[0.1, 0., 0.], dip=[0., 0., 1.],
                        rotation=[True, True, True])
        system.part.add(pos=[0.3, 0., 0.], dip=[0., 0.5, 0.5],
                        rotation=[True, True, True])
        ref_min_img_energy = 62.5
        ref_min_img_forces = 937.5 * np.array([[-1., 0., 0.], [+1., 0., 0.]])
        ref_min_img_torques = 62.5 * np.array([[+1., 0., 0.], [-1., 0., 0.]])
        system.actors.add(solver)

        # check min image against reference data
        system.periodicity = [False, False, False]
        system.integrator.run(steps=0, recalc_forces=True)
        min_img_energy = system.analysis.energy()["dipolar"]
        min_img_forces = np.copy(system.part.all().f)
        min_img_torques = np.copy(system.part.all().torque_lab)
        np.testing.assert_allclose(min_img_energy, ref_min_img_energy, **tol)
        np.testing.assert_allclose(min_img_forces, ref_min_img_forces, **tol)
        np.testing.assert_allclose(min_img_torques, ref_min_img_torques, **tol)

        # place particles far away from each other, but close across a boundary
        system.part.by_id(1).pos = [system.box_l[0] - 0.1, 0., 0.]

        # check that max image convention is applied for open boundaries
        system.periodicity = [False, False, False]
        system.integrator.run(steps=0, recalc_forces=True)
        max_img_energy = system.analysis.energy()["dipolar"]
        max_img_forces = np.copy(system.part.all().f)
        max_img_torques = np.copy(system.part.all().torque_lab)
        np.testing.assert_array_less(
            max_img_energy, min_img_energy / 1000.)
        np.testing.assert_array_less(
            np.abs(max_img_forces[:, 0]),
            np.abs(min_img_forces[:, 0]) / 1000.)
        np.testing.assert_array_less(
            np.abs(max_img_torques[:, 0]),
            np.abs(min_img_torques[:, 0]) / 1000.)
        np.testing.assert_allclose(
            max_img_torques,
            min_img_torques * max_img_energy / min_img_energy,
            **tol)

        # check that min image convention is applied for periodic boundaries
        system.periodicity = [True, True, True]
        system.integrator.run(steps=0, recalc_forces=True)
        min_img_energy = system.analysis.energy()["dipolar"]
        min_img_forces = np.copy(system.part.all().f)
        min_img_torques = np.copy(system.part.all().torque_lab)
        np.testing.assert_allclose(min_img_energy, ref_min_img_energy, **tol)
        np.testing.assert_allclose(min_img_forces, -ref_min_img_forces, **tol)
        np.testing.assert_allclose(min_img_torques, ref_min_img_torques, **tol)

    def test_min_image_convention_cpu(self):
        solver = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.)
        self.check_min_image_convention(solver, rtol=1e-10)

    @utx.skipIfMissingFeatures("DIPOLAR_DIRECT_SUM")
    @utx.skipIfMissingGPU()
    def test_min_image_convention_gpu(self):
        solver = espressomd.magnetostatics.DipolarDirectSumGpu(prefactor=1.)
        self.check_min_image_convention(solver, rtol=1e-5)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] == 1,
               "only runs for 2 or more MPI ranks")
    def test_inner_loop_consistency_cpu(self):
        system = self.system
        system.periodicity = [True, True, True]
        tol = {"atol": 1e-10, "rtol": 1e-10}
        p1 = system.part.add(pos=[0., 0., 0.], dip=[0., 0., 1.],
                             rotation=[True, True, True])
        p2 = system.part.add(pos=[1., 0., 0.], dip=[0., 0., 1.],
                             rotation=[True, True, True])
        for n_replicas in [0, 1]:
            system.actors.clear()
            solver = espressomd.magnetostatics.DipolarDirectSumCpu(
                prefactor=1., n_replicas=n_replicas)
            system.actors.add(solver)

            # intra-node calculation
            p1.pos = [system.box_l[0] / 2. - 0.1, 0., 2.]
            p2.pos = [system.box_l[0] / 2. + 0.1, 0., 0.]
            system.integrator.run(steps=0, recalc_forces=True)
            assert p1.node != p2.node
            node_01_energy = system.analysis.energy()["dipolar"]
            node_01_forces = np.copy(system.part.all().f)
            node_01_torques = np.copy(system.part.all().torque_lab)

            # inter-node calculation
            p1.pos = [0.1, 0., 2.]
            p2.pos = [0.3, 0., 0.]
            system.integrator.run(steps=0, recalc_forces=True)
            assert p1.node == p2.node
            node_00_energy = system.analysis.energy()["dipolar"]
            node_00_forces = np.copy(system.part.all().f)
            node_00_torques = np.copy(system.part.all().torque_lab)

            np.testing.assert_allclose(node_01_energy, node_00_energy, **tol)
            np.testing.assert_allclose(node_01_forces, node_00_forces, **tol)
            np.testing.assert_allclose(node_01_torques, node_00_torques, **tol)


if __name__ == "__main__":
    ut.main()
