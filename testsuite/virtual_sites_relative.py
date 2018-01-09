
#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond
from espressomd.virtual_sites import VirtualSitesRelative, VirtualSitesOff

from tests_common import verify_lj_forces
from numpy import random


@ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),
           "Test requires VIRTUAL_SITES_RELATIVE")
class VirtualSites(ut.TestCase):
    s = espressomd.System()
    s.seed = range(s.cell_system.get_state()["n_nodes"])

    @classmethod
    def setUpClass(cls):
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
            (quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3])))

    def verify_vs(self, vs,verify_velocity=True):
        """Verify vs position and (if compiled in) velocity."""
        self.assertEqual(vs.virtual, 1)

        vs_r = vs.vs_relative

        # Get related particle
        rel = self.s.part[vs_r[0]]

        # Distance
        d = self.s.distance(rel, vs)
        v_d = self.s.distance_vec(rel, vs)
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
        self.assertTrue(isinstance(self.s.virtual_sites, VirtualSitesOff))

        # Switch implementation
        self.s.virtual_sites=VirtualSitesRelative(have_velocity=False)
        self.assertTrue(isinstance(self.s.virtual_sites, VirtualSitesRelative))
        self.assertEqual(self.s.virtual_sites.have_velocity,False)

    
    def test_pos_vel_forces(self):
        s = self.s
        s.cell_system.skin=0.3
        s.virtual_sites=VirtualSitesRelative(have_velocity=True)
        s.box_l = 10,10,10
        s.part.clear()
        s.time_step = 0.004
        s.part.clear()
        s.thermostat.turn_off()
        s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0, shift=0)

        # Check setting of min_global_cut
        s.min_global_cut = 0.23
        self.assertEqual(s.min_global_cut, 0.23)

        # Place central particle + 3 vs
        s.part.add(rotation=(1,1,1),pos=(0.5, 0.5, 0.5), id=1, quat=(
            1, 0, 0, 0), omega_lab=(1, 2, 3))
        pos2 = (0.5, 0.4, 0.5)
        pos3 = (0.3, 0.5, 0.4)
        pos4 = (0.5, 0.5, 0.5)
        cur_id = 2
        for pos in pos2, pos3, pos4:
            s.part.add(rotation=(1,1,1), pos=pos, id=cur_id)
            s.part[cur_id].vs_auto_relate_to(1)
            # Was the particle made virtual
            self.assertEqual(s.part[cur_id].virtual, 1)
            # Are vs relative to id and
            vs_r = s.part[cur_id].vs_relative
            # id
            self.assertEqual(vs_r[0], 1)
            # distance
            self.assertAlmostEqual(vs_r[1], s.distance(
                s.part[1], s.part[cur_id]), places=6)
            cur_id += 1

        # Move central particle and Check vs placement
        s.part[1].pos = (0.22, 0.22, 0.22)
        # linear and rotation velocity on central particle
        s.part[1].v = (0.45, 0.14, 0.447)
        s.part[1].omega_lab = (0.45, 0.14, 0.447)
        s.integrator.run(0, recalc_forces=True)
        # Ceck
        for i in 2, 3, 4:
            self.verify_vs(s.part[i])

        # Check if still true, when non-virtual particle has rotated and a
        # linear motion
        s.part[1].omega_lab = -5, 3, 8.4
        s.integrator.run(10)
        for i in 2, 3, 4:
            self.verify_vs(s.part[i])

        # Test transfer of forces accumulating on virtual sites
        # to central particle
        f2 = np.array((3, 4, 5))
        f3 = np.array((-4, 5, 6))
        # Add forces to vs
        s.part[2].ext_force = f2
        s.part[3].ext_force = f3
        s.integrator.run(0)
        # get force/torques on non-vs
        f = s.part[1].f
        t = s.part[1].torque_lab

        # Expected force = sum of the forces on the vs
        self.assertLess(np.linalg.norm(f - f2 - f3), 1E-6)

        # Expected torque
        # Radial components of forces on a rigid body add to the torque
        t_exp = np.cross(s.distance_vec(s.part[1], s.part[2]), f2)
        t_exp += np.cross(s.distance_vec(s.part[1], s.part[3]), f3)
        # Check
        self.assertLessEqual(np.linalg.norm(t_exp - t), 1E-6)


        # Check virtual sites without velocity
        s.virtual_sites.have_velocity=False

        v2=s.part[2].v
        s.part[1].v=17,-13.5,2
        s.integrator.run(0,recalc_forces=True)
        self.assertLess(np.linalg.norm(v2-s.part[2].v),1E-6)


        


    def run_test_lj(self):
        """This fills the system with vs-based dumbells, adds a lj potential
          integrates and verifies forces. This is to make sure, that no pairs
          get lost or are outdated in the short range loop"""
        s = self.s
        s.virtual_sites=VirtualSitesRelative(have_velocity=True)
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
        s.box_l = l, l, l
        s.min_global_cut = 0.501
        s.part.clear()

        s.time_step = 0.01
        s.thermostat.turn_off()


        # Dumbells consist of 2 virtual lj spheres + central particle w/o interactions
        # For n sphers n/2 dumbells.
        for i in range(int(n / 2)):
            # Type=1, i.e., no lj ia for the center of mass particles
            s.part.add(rotation=(1,1,1), id=3 * i, pos=random.random(3) * l, type=1,
                       omega_lab=0.3 * random.random(3), v=random.random(3))
            # lj spheres
            s.part.add(rotation=(1,1,1), id=3 * i + 1,
                       pos=s.part[3 * i].pos + s.part[3 * i].director / 2.,
                       type=0)
            s.part.add(rotation=(1,1,1), id=3 * i + 2,
                       pos=s.part[3 * i].pos - s.part[3 * i].director / 2.,
                       type=0)
            s.part[3 * i + 1].vs_auto_relate_to(3 * i)
            self.verify_vs(s.part[3*i+1],verify_velocity=False)
            s.part[3 * i + 2].vs_auto_relate_to(3 * i)
            self.verify_vs(s.part[3*i+2],verify_velocity=False)
        s.integrator.run(0,recalc_forces=True)
        # interactions
        s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
        # Remove overlap
        s.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        while s.analysis.energy()["total"] > 10 * n:
            s.integrator.run(20)
        # Integrate
        s.integrator.set_vv()
        for i in range(10):
            # Langevin to maintain stability
            s.thermostat.set_langevin(kT=kT, gamma=gamma)
            s.integrator.run(300)
            s.thermostat.turn_off()
            # Constant energy to get rid of thermostat forces in the
            # verification
            s.integrator.run(2)
            # Theck the virtual sites config,pos and vel of the lj spheres
            for j in range(int(n / 2)):
                self.verify_vs(s.part[3 * j + 1])
                self.verify_vs(s.part[3 * j + 2])

            # Verify lj forces on the particles. The non-virtual particles are skipeed
            # because the forces on them originate from the vss and not the lj
            # interaction
            verify_lj_forces(s, 1E-10, 3 * np.arange(int(n / 2), dtype=int))

        # Turn off lj interaction
        s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0, shift=0)

    @ut.skipIf(
        espressomd.has_features("VIRTUAL_SITES_THERMOSTAT"),
        "LJ fluid test only works when VIRTUAL_SITES_THERMOSTAT is not compiled in.")
    def test_lj(self):
        """Run LJ fluid test for different cell systems."""
        s = self.s

        s.cell_system.skin = 0.4 
        s.cell_system.set_n_square(use_verlet_lists=True)
        self.run_test_lj()
        s.cell_system.set_domain_decomposition(use_verlet_lists=True)
        self.run_test_lj()
        s.cell_system.set_domain_decomposition(use_verlet_lists=False)
        self.run_test_lj()

    def test_zz_stress_tensor(self):
        s=self.s
        s.time_step=0.01
        s.cell_system.skin=0.1
        s.min_global_cut=0.2
        # Should not have one if vs are turned off
        s.virtual_sites=VirtualSitesOff()
        self.assertTrue("virtual_sites" not in s.analysis.pressure())
        self.assertTrue("virtual_sites" not in s.analysis.stress_tensor())

        # vs relative contrib
        s.virtual_sites=VirtualSitesRelative()
        s.part.clear()
        s.part.add(pos=(0,0,0),id=0)
        p=s.part.add(pos=(0.1,0.1,0.1),id=1,ext_force=(1,2,3))
        p.vs_auto_relate_to(0)
        p=s.part.add(pos=(0.1,0,0),id=2,ext_force=(-1,0,0))
        p.vs_auto_relate_to(0)
        s.integrator.run(0,recalc_forces=True)
        stress_total=s.analysis.stress_tensor()["total"]
        stress_vs_total=s.analysis.stress_tensor()["virtual_sites"]
        stress_vs=s.analysis.stress_tensor()["virtual_sites",0]
        
        p_total=s.analysis.pressure()["total"]
        p_vs_total=s.analysis.pressure()["virtual_sites"]
        p_vs=s.analysis.pressure()["virtual_sites",0]

        # expected stress
        s_expected =1./s.volume() * (
           np.outer(s.part[1].ext_force,s.distance_vec(s.part[1],s.part[0]))
          +np.outer(s.part[2].ext_force,s.distance_vec(s.part[2],s.part[0])))
        np.testing.assert_allclose(stress_total,s_expected,atol=1E-5)
        np.testing.assert_allclose(stress_vs,s_expected,atol=1E-5)
        np.testing.assert_allclose(stress_vs_total,s_expected,atol=1E-5)
        
        # Pressure
        self.assertAlmostEqual(p_total,np.sum(np.diag(s_expected))/3,places=5)
        self.assertAlmostEqual(p_vs_total,np.sum(np.diag(s_expected))/3,places=5)
        self.assertAlmostEqual(p_vs,np.sum(np.diag(s_expected))/3,places=5)








if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
