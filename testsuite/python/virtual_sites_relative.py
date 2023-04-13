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
import espressomd
import espressomd.virtual_sites
import numpy as np
import tests_common


@utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE", "LENNARD_JONES"])
class VirtualSites(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    np.random.seed(42)

    def setUp(self):
        self.system.box_l = [10.0, 10.0, 10.0]
        self.system.cell_system.set_regular_decomposition(
            use_verlet_lists=True)

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesOff()

    def multiply_quaternions(self, a, b):
        return np.array(
            (a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
             a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
             a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3],
             a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1]))

    def director_from_quaternion(self, quat):
        return np.array((
            2 * (quat[1] * quat[3] + quat[0] * quat[2]),
            2 * (quat[2] * quat[3] - quat[0] * quat[1]),
            (quat[0] * quat[0] - quat[1] * quat[1]
             - quat[2] * quat[2] + quat[3] * quat[3])))

    def verify_vs(self, vs, verify_velocity=True):
        """Verify virtual site position and velocity."""
        self.assertTrue(vs.virtual)

        vs_r = vs.vs_relative

        # Get related particle
        rel = self.system.part.by_id(vs_r[0])

        # Distance
        d = self.system.distance(rel, vs)
        v_d = self.system.distance_vec(rel, vs)
        # Check distance
        self.assertAlmostEqual(d, vs_r[1], places=6)

        # check velocity
        if verify_velocity:
            self.assertLessEqual(np.linalg.norm(
                vs.v - rel.v - np.cross(rel.omega_lab, v_d)), 1E-6)

        # Check position
        self.assertLess(np.linalg.norm(
            v_d - vs_r[1] * self.director_from_quaternion(
                self.multiply_quaternions(rel.quat, vs_r[2]))), 1E-6)

    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertIsInstance(
            self.system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesOff)
        self.assertFalse(self.system.virtual_sites.have_quaternion)
        self.assertFalse(self.system.virtual_sites.override_cutoff_check)

        # Set properties
        self.system.virtual_sites.have_quaternion = True
        self.system.virtual_sites.override_cutoff_check = True
        self.assertTrue(self.system.virtual_sites.have_quaternion)
        self.assertTrue(self.system.virtual_sites.override_cutoff_check)

        # Switch implementation
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        self.assertIsInstance(
            self.system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesRelative)
        self.assertFalse(self.system.virtual_sites.have_quaternion)
        self.assertFalse(self.system.virtual_sites.override_cutoff_check)

    def test_vs_quat(self):
        self.system.time_step = 0.01
        self.system.min_global_cut = 0.23
        # First check that quaternion of virtual particle is unchanged if
        # have_quaternion is false.
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
            have_quaternion=False)
        self.assertFalse(self.system.virtual_sites.have_quaternion)
        p1 = self.system.part.add(pos=[1, 1, 1], rotation=3 * [True],
                                  omega_lab=[1, 1, 1])
        p2 = self.system.part.add(pos=[1, 1, 1], rotation=3 * [True])
        p2.vs_auto_relate_to(p1)
        np.testing.assert_array_equal(np.copy(p2.quat), [1, 0, 0, 0])
        self.system.integrator.run(1)
        np.testing.assert_array_equal(np.copy(p2.quat), [1, 0, 0, 0])
        # Now check that quaternion of the virtual particle gets updated.
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
            have_quaternion=True)
        self.system.integrator.run(1)
        self.assertRaises(AssertionError, np.testing.assert_array_equal,
                          np.copy(p2.quat), [1, 0, 0, 0])

        # co-aligned case
        p2.vs_quat = (1, 0, 0, 0)
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p2.director), np.copy(p1.director), atol=1E-12)

        # Construct a quaternion with perpendicular orientation.
        p0 = np.cos(np.pi / 4.0)
        p = np.array([0, np.sin(np.pi / 4.0), 0])
        q0 = 1.0
        q = np.array([0, 0, 0])
        r0 = p0 * q0
        r = -np.dot(q, p) + np.cross(q, p) + p0 * q + q0 * p
        p2.vs_quat = [r0, r[0], r[1], r[2]]
        self.system.integrator.run(1)
        # Check for orthogonality.
        self.assertAlmostEqual(
            np.dot(p1.director, p2.director), 0.0, delta=1E-12)
        # Check if still true after integration.
        self.system.integrator.run(1)
        self.assertAlmostEqual(np.dot(p1.director, p2.director), 0.0)

    def test_vs_exceptions(self):
        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        system.min_global_cut = 0.1
        p1 = system.part.add(pos=[0.0, 0.0, 0.0], rotation=3 * [True], id=1)
        p2 = system.part.add(pos=[1.0, 1.0, 1.0], rotation=3 * [True], id=2)
        p3 = system.part.add(pos=[1.0, 1.0, 1.0], rotation=3 * [True], id=3)
        p4 = system.part.add(pos=[1.0, 1.0, 1.0], rotation=3 * [True], id=4)
        # dangling virtual sites are not allowed
        with self.assertRaisesRegex(Exception, "Particle with id 4 is a dangling virtual site"):
            p4.virtual = True
            self.assertEqual(p4.vs_relative[0], -1)
            system.integrator.run(0, recalc_forces=True)
        p4.remove()
        # relating to anything else other than a particle or id is not allowed
        with self.assertRaisesRegex(ValueError, "Argument of 'vs_auto_relate_to' has to be of type ParticleHandle or int"):
            p2.vs_auto_relate_to('0')
        with self.assertRaisesRegex(ValueError, "Invalid particle id: -2"):
            p2.vs_auto_relate_to(-2)
        # relating to itself is not allowed
        with self.assertRaisesRegex(ValueError, "A virtual site cannot relate to itself"):
            p2.vs_auto_relate_to(p2)
        # relating to a deleted particle is not allowed
        with self.assertRaisesRegex(Exception, "No real particle with id 3 for virtual site with id 2"):
            p2.vs_auto_relate_to(p3)
            p3.remove()
            system.integrator.run(0, recalc_forces=True)
        if system.cell_system.get_state()["n_nodes"] > 1:
            with self.assertRaisesRegex(Exception, r"The distance between virtual and non-virtual particle \([0-9\.]+\) is larger than the minimum global cutoff"):
                p2.vs_auto_relate_to(p1)
            # If overridden this check should not raise an exception
            system.virtual_sites.override_cutoff_check = True
            p2.vs_auto_relate_to(p1)

    def test_pos_vel_forces(self):
        system = self.system
        system.cell_system.skin = 0.3
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        system.time_step = 0.004

        # Check setting of min_global_cut
        system.min_global_cut = 0.23
        self.assertEqual(system.min_global_cut, 0.23)

        # Place central particle + 3 vs
        p1 = system.part.add(rotation=3 * [True], pos=(0.5, 0.5, 0.5), id=1,
                             quat=(1, 0, 0, 0), omega_lab=(1, 2, 3))
        pos2 = (0.5, 0.4, 0.5)
        pos3 = (0.3, 0.5, 0.4)
        pos4 = (0.5, 0.5, 0.5)
        for pos in (pos2, pos3, pos4):
            p = system.part.add(rotation=3 * [True], pos=pos)
            p.vs_auto_relate_to(p1)
            # Was the particle made virtual
            self.assertTrue(p.virtual)
            # Are vs relative to id and
            vs_r = p.vs_relative
            # id
            self.assertEqual(vs_r[0], p1.id)
            # distance
            self.assertAlmostEqual(vs_r[1], system.distance(p1, p), places=6)

        # Move central particle and check vs placement
        p1.pos = (0.22, 0.22, 0.22)
        # Linear and rotational velocity on central particle
        p1.v = (0.45, 0.14, 0.447)
        p1.omega_lab = (0.45, 0.14, 0.447)
        system.integrator.run(0, recalc_forces=True)
        for p in system.part:
            if p.id != p1.id:
                self.verify_vs(p)

        # Check if still true, when non-virtual particle has rotated and a
        # linear motion
        p1.omega_lab = [-5., 3., 8.4]
        system.integrator.run(10)
        for p in system.part:
            if p.id != p1.id:
                self.verify_vs(p)

        if espressomd.has_features("EXTERNAL_FORCES"):
            # Test transfer of forces accumulating on virtual sites
            # to central particle
            f2 = np.array((3, 4, 5))
            f3 = np.array((-4, 5, 6))
            # Add forces to vs
            p2, p3 = system.part.by_ids([2, 3])
            p2.ext_force = f2
            p3.ext_force = f3
            system.integrator.run(0)
            # get force/torques on non-vs
            f = p1.f
            t = p1.torque_lab

            # Expected force = sum of the forces on the vs
            self.assertAlmostEqual(np.linalg.norm(f - f2 - f3), 0., delta=1E-6)

            # Expected torque
            # Radial components of forces on a rigid body add to the torque
            t_exp = np.cross(system.distance_vec(p1, p2), f2)
            t_exp += np.cross(system.distance_vec(p1, p3), f3)
            # Check
            self.assertAlmostEqual(np.linalg.norm(t_exp - t), 0., delta=1E-6)

    def run_test_lj(self):
        """
        This fills the system with vs-based dumbbells, adds a LJ potential,
        integrates and verifies forces. This is to make sure that no pairs
        get lost or are outdated in the short range loop.
        """
        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        # Parameters
        n = 40
        phi = 0.6
        sigma = 1.
        eps = .025
        cut = sigma * 2**(1. / 6.)

        kT = 2
        gamma = .5

        # box
        l = np.cbrt(n / 6. * np.pi * sigma**3 / phi)

        # Setup
        system.box_l = [l, l, l]
        system.min_global_cut = 0.501
        system.time_step = 0.01

        # Dumbbells consist of 2 virtual lj spheres + central particle
        # w/o interactions. For n spheres, n/2 dumbbells.
        for i in range(n // 2):
            # Type=1, i.e., no lj ia for the center of mass particles
            p3i = system.part.add(
                rotation=3 * [True], id=3 * i, pos=np.random.random(3) * l, type=1,
                omega_lab=0.3 * np.random.random(3), v=np.random.random(3))
            # lj spheres
            p3ip1 = system.part.add(rotation=3 * [True],
                                    id=3 * i + 1,
                                    pos=p3i.pos +
                                    p3i.director / 2.,
                                    type=0)
            p3ip2 = system.part.add(rotation=3 * [True],
                                    id=3 * i + 2,
                                    pos=p3i.pos -
                                    p3i.director / 2.,
                                    type=0)
            p3ip1.vs_auto_relate_to(p3i.id)
            self.verify_vs(p3ip1, verify_velocity=False)
            p3ip2.vs_auto_relate_to(p3i.id)
            self.verify_vs(p3ip2, verify_velocity=False)
        system.integrator.run(0, recalc_forces=True)
        # interactions
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
        # Remove overlap
        system.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        n_loops = 0
        n_max = 10
        while system.analysis.energy()["total"] > 10 * n and n_loops < n_max:
            system.integrator.run(20)
            n_loops += 1
        assert n_loops < n_max, "Steepest descent didn't converge"
        # Integrate
        system.integrator.set_vv()
        for i in range(10):
            # Langevin to maintain stability
            system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
            system.integrator.run(50)
            system.thermostat.turn_off()
            # Constant energy to get rid of thermostat forces in the
            # verification
            system.integrator.run(2)
            # Check the virtual sites config, pos and vel of the lj spheres
            for j in range(int(n / 2)):
                self.verify_vs(system.part.by_id(3 * j + 1))
                self.verify_vs(system.part.by_id(3 * j + 2))

            # Verify lj forces on the particles. The non-virtual particles are
            # skipped because the forces on them originate from the vss and not
            # the lj interaction
            tests_common.verify_lj_forces(system, 1E-10, 3 * np.arange(n // 2))

        # Test applying changes
        energy_pre_change = system.analysis.energy()['total']
        pressure_pre_change = system.analysis.pressure()['total']
        p0 = system.part.by_id(0)
        p0.pos = p0.pos + (2.2, -1.4, 4.2)
        energy_post_change = system.analysis.energy()['total']
        pressure_post_change = system.analysis.pressure()['total']
        self.assertNotAlmostEqual(energy_pre_change, energy_post_change)
        self.assertNotAlmostEqual(pressure_pre_change, pressure_post_change)

    def test_lj(self):
        """Run LJ fluid test for different cell systems."""

        self.system.cell_system.skin = 0.4
        with self.subTest(msg='N-square cell system with Verlet lists'):
            self.system.cell_system.set_n_square(use_verlet_lists=True)
            self.run_test_lj()
            self.tearDown()
        with self.subTest(msg='regular decomposition cell system with Verlet lists'):
            self.system.cell_system.set_regular_decomposition(
                use_verlet_lists=True)
            self.run_test_lj()
            self.tearDown()
        with self.subTest(msg='regular decomposition cell system without Verlet lists'):
            self.system.cell_system.set_regular_decomposition(
                use_verlet_lists=False)
            self.run_test_lj()
            self.tearDown()

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_zz_pressure_tensor(self):
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        system.min_global_cut = 0.2
        # Should not have a pressure
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesOff()
        pressure_tensor_vs = system.analysis.pressure_tensor()[
            "virtual_sites", 0]
        p_vs = system.analysis.pressure()["virtual_sites", 0]
        np.testing.assert_allclose(pressure_tensor_vs, 0., atol=1e-10)
        np.testing.assert_allclose(p_vs, 0., atol=1e-10)

        # vs relative contrib
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        p0 = system.part.add(id=0, pos=(0.0, 0.0, 0.0))
        p1 = system.part.add(id=1, pos=(0.1, 0.1, 0.1), ext_force=(1, 2, 3))
        p2 = system.part.add(id=2, pos=(0.1, 0.0, 0.0), ext_force=(-1, 0, 0))
        p1.vs_auto_relate_to(p0)
        p2.vs_auto_relate_to(p0)
        system.integrator.run(0, recalc_forces=True)
        sim_pressure = system.analysis.pressure()
        sim_pressure_tensor = system.analysis.pressure_tensor()
        pressure_tensor_total = sim_pressure_tensor["total"]
        pressure_tensor_vs_total = sim_pressure_tensor["virtual_sites"]
        pressure_tensor_vs = sim_pressure_tensor["virtual_sites", 0]
        p_total = sim_pressure["total"]
        p_vs_total = sim_pressure["virtual_sites"]
        p_vs = sim_pressure["virtual_sites", 0]

        # expected pressure_tensor
        s_expected = 1. / system.volume() * (
            np.outer(p1.ext_force, system.distance_vec(p1, p0))
            + np.outer(p2.ext_force, system.distance_vec(p2, p0)))
        np.testing.assert_allclose(
            pressure_tensor_total, s_expected, atol=1E-5)
        np.testing.assert_allclose(pressure_tensor_vs, s_expected, atol=1E-5)
        np.testing.assert_allclose(
            pressure_tensor_vs_total, s_expected, atol=1E-5)

        # Pressure
        self.assertAlmostEqual(p_total, np.sum(
            np.diag(s_expected)) / 3, places=5)
        self.assertAlmostEqual(p_vs_total, np.sum(
            np.diag(s_expected)) / 3, places=5)
        self.assertAlmostEqual(p_vs, np.sum(np.diag(s_expected)) / 3, places=5)


if __name__ == "__main__":
    ut.main(verbosity=2)
