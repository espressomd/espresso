#
# Copyright (C) 2013-2019 The ESPResSo project
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
if espressomd.has_features("VIRTUAL_SITES_RELATIVE"):
    from espressomd.virtual_sites import VirtualSitesRelative, VirtualSitesOff
import numpy as np

from tests_common import verify_lj_forces
from numpy import random


@utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
class VirtualSites(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    np.random.seed(42)

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
        """Verify vs position and (if compiled in) velocity."""
        self.assertTrue(vs.virtual)

        vs_r = vs.vs_relative

        # Get related particle
        rel = self.system.part[vs_r[0]]

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
        self.assertIsInstance(self.system.virtual_sites, VirtualSitesOff)

        # Switch implementation
        self.system.virtual_sites = VirtualSitesRelative()
        self.assertIsInstance(self.system.virtual_sites, VirtualSitesRelative)

    def test_vs_quat(self):
        self.system.part.clear()
        # First check that quaternion of virtual particle is unchanged if
        # have_quaterion is false.
        self.system.virtual_sites = VirtualSitesRelative(have_quaternion=False)
        self.assertFalse(self.system.virtual_sites.have_quaternion)
        self.system.part.add(id=0, pos=[1, 1, 1], rotation=[1, 1, 1],
                             omega_lab=[1, 1, 1])
        self.system.part.add(id=1, pos=[1, 1, 1], rotation=[1, 1, 1])
        self.system.part[1].vs_auto_relate_to(0)
        np.testing.assert_array_equal(
            np.copy(self.system.part[1].quat), [1, 0, 0, 0])
        self.system.integrator.run(1)
        np.testing.assert_array_equal(
            np.copy(self.system.part[1].quat), [1, 0, 0, 0])
        # Now check that quaternion of the virtual particle gets updated.
        self.system.virtual_sites = VirtualSitesRelative(have_quaternion=True)
        self.system.integrator.run(1)
        self.assertRaises(AssertionError, np.testing.assert_array_equal, np.copy(
            self.system.part[1].quat), [1, 0, 0, 0])

        # co-aligned case
        self.system.part[1].vs_quat = (1, 0, 0, 0)
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(self.system.part[1].director), np.copy(self.system.part[0].director), atol=1E-12)

        # Construct a quaternion with perpendicular orientation.
        p0 = np.cos(np.pi / 4.0)
        p = np.array([0, np.sin(np.pi / 4.0), 0])
        q0 = 1.0
        q = np.array([0, 0, 0])
        r0 = p0 * q0
        r = -np.dot(q, p) + np.cross(q, p) + p0 * q + q0 * p
        self.system.part[1].vs_quat = [r0, r[0], r[1], r[2]]
        self.system.integrator.run(1)
        # Check for orthogonality.
        self.assertAlmostEqual(
            np.dot(self.system.part[0].director, self.system.part[1].director), 0.0, delta=1E-12)
        # Check if still true after integration.
        self.system.integrator.run(1)
        self.assertAlmostEqual(
            np.dot(self.system.part[0].director, self.system.part[1].director), 0.0)

    def test_pos_vel_forces(self):
        system = self.system
        system.cell_system.skin = 0.3
        system.virtual_sites = VirtualSitesRelative()
        system.box_l = [10, 10, 10]
        system.part.clear()
        system.time_step = 0.004
        system.part.clear()
        system.thermostat.turn_off()
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0, shift=0)

        # Check setting of min_global_cut
        system.min_global_cut = 0.23
        self.assertEqual(system.min_global_cut, 0.23)

        # Place central particle + 3 vs
        system.part.add(rotation=(1, 1, 1), pos=(0.5, 0.5, 0.5), id=1,
                        quat=(1, 0, 0, 0), omega_lab=(1, 2, 3))
        pos2 = (0.5, 0.4, 0.5)
        pos3 = (0.3, 0.5, 0.4)
        pos4 = (0.5, 0.5, 0.5)
        cur_id = 2
        for pos in pos2, pos3, pos4:
            system.part.add(rotation=(1, 1, 1), pos=pos, id=cur_id)
            system.part[cur_id].vs_auto_relate_to(1)
            # Was the particle made virtual
            self.assertTrue(system.part[cur_id].virtual)
            # Are vs relative to id and
            vs_r = system.part[cur_id].vs_relative
            # id
            self.assertEqual(vs_r[0], 1)
            # distance
            self.assertAlmostEqual(vs_r[1], system.distance(
                system.part[1], system.part[cur_id]), places=6)
            cur_id += 1

        # Move central particle and check vs placement
        system.part[1].pos = (0.22, 0.22, 0.22)
        # Linear and rotational velocity on central particle
        system.part[1].v = (0.45, 0.14, 0.447)
        system.part[1].omega_lab = (0.45, 0.14, 0.447)
        system.integrator.run(0, recalc_forces=True)
        for i in [2, 3, 4]:
            self.verify_vs(system.part[i])

        # Check if still true, when non-virtual particle has rotated and a
        # linear motion
        system.part[1].omega_lab = [-5, 3, 8.4]
        system.integrator.run(10)
        for i in [2, 3, 4]:
            self.verify_vs(system.part[i])

        if espressomd.has_features("EXTERNAL_FORCES"):
            # Test transfer of forces accumulating on virtual sites
            # to central particle
            f2 = np.array((3, 4, 5))
            f3 = np.array((-4, 5, 6))
            # Add forces to vs
            system.part[2].ext_force = f2
            system.part[3].ext_force = f3
            system.integrator.run(0)
            # get force/torques on non-vs
            f = system.part[1].f
            t = system.part[1].torque_lab

            # Expected force = sum of the forces on the vs
            self.assertLess(np.linalg.norm(f - f2 - f3), 1E-6)

            # Expected torque
            # Radial components of forces on a rigid body add to the torque
            t_exp = np.cross(system.distance_vec(
                system.part[1], system.part[2]), f2)
            t_exp += np.cross(system.distance_vec(
                system.part[1], system.part[3]), f3)
            # Check
            self.assertLessEqual(np.linalg.norm(t_exp - t), 1E-6)

    def run_test_lj(self):
        """This fills the system with vs-based dumbells, adds a lj potential,
          integrates and verifies forces. This is to make sure that no pairs
          get lost or are outdated in the short range loop"""
        system = self.system
        system.virtual_sites = VirtualSitesRelative()
        # Parameters
        n = 90
        phi = 0.6
        sigma = 1.
        eps = .025
        cut = sigma * 2**(1. / 6.)

        kT = 2
        gamma = .5

        # box
        l = (n / 6. * np.pi * sigma**3 / phi)**(1. / 3.)

        # Setup
        system.box_l = l, l, l
        system.min_global_cut = 0.501
        system.part.clear()

        system.time_step = 0.01
        system.thermostat.turn_off()

        # Dumbells consist of 2 virtual lj spheres + central particle w/o interactions
        # For n spheres, n/2 dumbells.
        for i in range(int(n / 2)):
            # Type=1, i.e., no lj ia for the center of mass particles
            system.part.add(
                rotation=(1, 1, 1), id=3 * i, pos=random.random(3) * l, type=1,
                omega_lab=0.3 * random.random(3), v=random.random(3))
            # lj spheres
            system.part.add(rotation=(1, 1, 1), id=3 * i + 1,
                            pos=system.part[3 * i].pos +
                            system.part[3 * i].director / 2.,
                            type=0)
            system.part.add(rotation=(1, 1, 1), id=3 * i + 2,
                            pos=system.part[3 * i].pos -
                            system.part[3 * i].director / 2.,
                            type=0)
            system.part[3 * i + 1].vs_auto_relate_to(3 * i)
            self.verify_vs(system.part[3 * i + 1], verify_velocity=False)
            system.part[3 * i + 2].vs_auto_relate_to(3 * i)
            self.verify_vs(system.part[3 * i + 2], verify_velocity=False)
        system.integrator.run(0, recalc_forces=True)
        # interactions
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
        # Remove overlap
        system.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        while system.analysis.energy()["total"] > 10 * n:
            system.integrator.run(20)
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
            # Check the virtual sites config,pos and vel of the lj spheres
            for j in range(int(n / 2)):
                self.verify_vs(system.part[3 * j + 1])
                self.verify_vs(system.part[3 * j + 2])

            # Verify lj forces on the particles. The non-virtual particles are
            # skipped because the forces on them originate from the vss and not
            # the lj interaction
            verify_lj_forces(system, 1E-10, 3 *
                             np.arange(int(n / 2), dtype=int))

        # Test applying changes
        enegry_pre_change = system.analysis.energy()['total']
        pressure_pre_change = system.analysis.pressure()['total']
        system.part[0].pos = system.part[0].pos + (2.2, -1.4, 4.2)
        enegry_post_change = system.analysis.energy()['total']
        pressure_post_change = system.analysis.pressure()['total']
        self.assertNotAlmostEqual(enegry_pre_change, enegry_post_change)
        self.assertNotAlmostEqual(pressure_pre_change, pressure_post_change)

        # Turn off lj interaction
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0, shift=0)

    def test_lj(self):
        """Run LJ fluid test for different cell systems."""
        system = self.system

        system.cell_system.skin = 0.4
        system.cell_system.set_n_square(use_verlet_lists=True)
        self.run_test_lj()
        system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        self.run_test_lj()
        system.cell_system.set_domain_decomposition(use_verlet_lists=False)
        self.run_test_lj()

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_zz_stress_tensor(self):
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        system.min_global_cut = 0.2
        # Should not have a pressure
        system.virtual_sites = VirtualSitesOff()
        stress_vs = system.analysis.stress_tensor()["virtual_sites", 0]
        p_vs = system.analysis.pressure()["virtual_sites", 0]
        np.testing.assert_allclose(stress_vs, 0., atol=1e-10)
        np.testing.assert_allclose(p_vs, 0., atol=1e-10)

        # vs relative contrib
        system.virtual_sites = VirtualSitesRelative()
        system.part.clear()
        system.part.add(pos=(0, 0, 0), id=0)
        p = system.part.add(pos=(0.1, 0.1, 0.1), id=1, ext_force=(1, 2, 3))
        p.vs_auto_relate_to(0)
        p = system.part.add(pos=(0.1, 0, 0), id=2, ext_force=(-1, 0, 0))
        p.vs_auto_relate_to(0)
        system.integrator.run(0, recalc_forces=True)
        stress_total = system.analysis.stress_tensor()["total"]
        stress_vs_total = system.analysis.stress_tensor()["virtual_sites"]
        stress_vs = system.analysis.stress_tensor()["virtual_sites", 0]

        p_total = system.analysis.pressure()["total"]
        p_vs_total = system.analysis.pressure()["virtual_sites"]
        p_vs = system.analysis.pressure()["virtual_sites", 0]

        # expected stress
        s_expected = 1. / system.volume() * (
            np.outer(system.part[1].ext_force, system.distance_vec(
                system.part[1], system.part[0]))
            + np.outer(system.part[2].ext_force, system.distance_vec(system.part[2], system.part[0])))
        np.testing.assert_allclose(stress_total, s_expected, atol=1E-5)
        np.testing.assert_allclose(stress_vs, s_expected, atol=1E-5)
        np.testing.assert_allclose(stress_vs_total, s_expected, atol=1E-5)

        # Pressure
        self.assertAlmostEqual(p_total, np.sum(
            np.diag(s_expected)) / 3, places=5)
        self.assertAlmostEqual(p_vs_total, np.sum(
            np.diag(s_expected)) / 3, places=5)
        self.assertAlmostEqual(p_vs, np.sum(np.diag(s_expected)) / 3, places=5)


if __name__ == "__main__":
    ut.main()
