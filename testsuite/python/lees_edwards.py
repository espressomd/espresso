import espressomd
from espressomd.interactions import HarmonicBond
import espressomd.lees_edwards as lees_edwards
from espressomd.virtual_sites import VirtualSitesRelative, VirtualSitesOff
from tests_common import verify_lj_forces

import unittest as ut
import numpy as np


class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=True)

    time_step = 0.5
    system.time_step = time_step

    def test_a_Protocol(self):
        """This test calculates if the offset and velocity of the Lees Edwards
        function are calculated correctly based on the input variables"""

        system = self.system

        # Protocol should be off by default
        self.assertTrue(self.system.lees_edwards.protocol is None)

        # Check if the set protocol is stored correctly
        params_lin = {
            'shear_velocity': 1.2,
            'initial_pos_offset': 0.1,
            'time_0': 0.1}
        system.time = 3.15
        lin_protocol = lees_edwards.LinearShear(**params_lin)
        system.lees_edwards.protocol = lin_protocol
        self.assertEqual(self.system.lees_edwards.protocol, lin_protocol)
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params_lin)
        expected_pos = params_lin['initial_pos_offset'] + \
            params_lin['shear_velocity'] * (system.time - params_lin['time_0'])
        expected_vel = params_lin['shear_velocity']
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity,
            expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check if the offset is determined correctly

        system.time = 0.0

        params_osc = {'amplitude': 2.3, 'omega': 2.51, 'time_0': -2.1}
        system.lees_edwards.protocol = lees_edwards.OscillatoryShear(
            **params_osc)
        for time in [0, 2.3]:

            system.time = time
            expected_pos = params_osc['amplitude'] * \
                np.sin(params_osc['omega'] * (time - params_osc['time_0']))
            expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
                np.cos(params_osc['omega'] * (time - params_osc['time_0']))
            self.assertAlmostEqual(
                system.lees_edwards.shear_velocity,
                expected_vel)
            self.assertAlmostEqual(
                system.lees_edwards.pos_offset, expected_pos)

        system.lees_edwards.protocol = None
        self.assertEqual(
            system.lees_edwards.shear_velocity, 0)
        self.assertEqual(
            system.lees_edwards.pos_offset, 0)


   def disable_test_b_BoundaryCrossing(self):
       """A particle crosses the upper and lower boundary to test if position
       and velocity are updated correctly."""
#
       system = self.system
       system.part.clear()
       system.time = 0.0
#
       # Set up a one particle system and check the position offset after crossing the boundary
       # Test for upper boundary
#
       velocity = 0.1
       dir = [0, 1, 2]
#
       for sheardir in dir:
#
   def disable_test_b_BoundaryCrossing(self):
       """A particle crosses the upper and lower boundary to test if position
       and velocity are updated correctly."""
#
       system = self.system
       system.part.clear()
       system.time = 0.0
#
       # Set up a one particle system and check the position offset after crossing the boundary
       # Test for upper boundary
#
       velocity = 0.1
       dir = [0, 1, 2]
#
       for sheardir in dir:
           for shear_plane_normal in dir:
               if sheardir != shear_plane_normal:
#
                   system.time = 0.0
                   system.lees_edwards.protocol = \
                       lees_edwards.LinearShear(
                           shear_velocity=velocity, 
                         shear_direction=sheardir, 
                         shear_plane_normal=shear_plane_normal)
#
                   pos = np.full([3], 2.5)
                   pos[shear_plane_normal] = system.box_l[
                       shear_plane_normal] - 0.05
                   vel = np.zeros([3])
                   vel[shear_plane_normal] = 0.2
                   system.part.add(pos=pos, v=vel, id=0, type=0)
#
                   pos_change = np.zeros([3])
                   pos_change[sheardir] = -system.time_step * 0.5 * velocity
                   pos_change[
                       shear_plane_normal] = vel[
                           shear_plane_normal] * system.time_step
                   vel_change = np.zeros([3])
                   vel_change[sheardir] = -velocity
#
                   expected_pos = system.part[0].pos + pos_change
                   expected_vel = system.part[0].v + vel_change
#
                   system.integrator.run(1)
#
                   np.testing.assert_almost_equal(
                       np.copy(system.part[0].pos), expected_pos)
                   np.testing.assert_almost_equal(
                       np.copy(system.part[0].v), expected_vel)
#
                   system.part.clear()
#
                   # Test for lower boundary
#
                   system.time = 0.0
                   pos = np.full([3], 2.5)
                   pos[shear_plane_normal] = 0.05
                   vel = np.zeros([3])
                   vel[shear_plane_normal] = -0.2
                   system.part.add(pos=pos, v=vel, id=0, type=0)
#
                   pos_change = np.zeros([3])
                   pos_change[sheardir] = system.time_step * 0.5 * velocity
                   pos_change[
                       shear_plane_normal] = vel[
                           shear_plane_normal] * system.time_step
                   vel_change = np.zeros([3])
                   vel_change[sheardir] = velocity
#
                   expected_pos = system.part[0].pos + pos_change
                   expected_vel = system.part[0].v + vel_change
#
                   system.integrator.run(1)
#
                   np.testing.assert_almost_equal(
                       np.copy(system.part[0].pos), expected_pos)
                   np.testing.assert_almost_equal(
                       np.copy(system.part[0].v), expected_vel)
#
                   system.part.clear()
#
   def disable_test_c_interactions(self):
       """We place two particles crossing a boundary and connect them with an unbonded
          and bonded interaction. We test with the resulting stress tensor if the offset
          is included properly"""
#
       system = self.system
       system.time = 0.0
       system.part.clear()
#
       dir = [0, 1, 2]
#
       for sheardir in dir:
           for shear_plane_normal in dir:
               if sheardir != shear_plane_normal:
#
                   offset = 1.0
                   system.lees_edwards.protocol = \
                       lees_edwards.LinearShear(
                           initial_pos_offset=offset,
                           shear_velocity=0,
                           shear_direction=sheardir,
                           shear_plane_normal=shear_plane_normal)
#
                   pos1 = np.full([3], 2.5)
                   pos1[shear_plane_normal] = 4.75
                   pos2 = np.full([3], 2.5)
                   pos2[shear_plane_normal] = 0.25
#
                   p1 = system.part.add(id=0, pos=pos1, fix=[1, 1, 1])
                   p2 = system.part.add(id=1, pos=pos2, fix=[1, 1, 1])
#
                   r = system.part[1].pos - system.part[0].pos
                   r[sheardir] += offset
                   r[shear_plane_normal] += system.box_l[0]
#
                   k = 1000
                   r_cut = 1.5
                   harm = HarmonicBond(k=k, r_0=0.0)
                   system.bonded_inter.add(harm)
                   system.part[0].add_bond((harm, 1))
                   system.non_bonded_inter[0, 0].soft_sphere.set_params(
                       a=k, n=-2, cutoff=r_cut)
#
                   system.integrator.run(1, recalc_forces=True)
                   np.testing.assert_allclose(
                     system.distance_vec(p1, p2), r)
                   np.testing.assert_allclose(
                     system.distance_vec(p2, p1), -r)
#
                   simulated_bondedstress = np.abs(
                       system.analysis.stress_tensor()['bonded'])
                   simulated_nonbondedstress = np.abs(
                       system.analysis.stress_tensor()['non_bonded'])
#
                   analytical_bondedstress = np.zeros([3, 3])
                   for i in range(3):
                       for j in range(3):
                           analytical_bondedstress[i, j] += k * r[i] * r[j]
                   analytical_bondedstress /= (system.box_l[0]**3.0)
#
                   analytical_nonbondedstress = np.zeros([3, 3])
                   for i in range(3):
                       for j in range(3):
                           analytical_nonbondedstress[
                               i, j] += 2 * k * r[i] * r[j]
                   analytical_nonbondedstress /= (system.box_l[0]**3.0)
#
                   np.testing.assert_array_equal(
                       simulated_bondedstress, analytical_bondedstress)
                   np.testing.assert_array_equal(
                       simulated_nonbondedstress, analytical_nonbondedstress)
#
                   system.part.clear()
                   
                   # Turn off soft sphere interaction
                   system.non_bonded_inter[0, 0].soft_sphere.set_params(
                       a=0, n=0, cutoff=0)
#
   def disable_test_d_vel_diff(self):
       """We place particles within a box and across a box boundary to check if
          the Lees Edwards shift is taken into account correctly."""
#
       system = self.system
       system.part.clear()
#
       system.lees_edwards.protocol = \
           lees_edwards.LinearShear(
               shear_velocity=1.5,
               shear_direction=0,
               shear_plane_normal=2)
       p1 = system.part.add(
           id=0, pos=[1, 1, 0.9 * system.box_l[2]], v=[2, -3.2, 3])
       p2 = system.part.add(
           id=1, pos=[1, 1, 0.1 * system.box_l[2]], v=[7, 6.1, -1])
#
       vel_diff = np.copy(system.velocity_difference(p1, p2))
       np.testing.assert_array_equal(
           vel_diff, p2.v - p1.v + np.array([system.lees_edwards.shear_velocity, 0, 0]))
#
       system.part.clear()
       system.lees_edwards.protocol = lees_edwards.Off()
#
       p1 = system.part.add(id=0, pos=[1, 1, 0.5], v=[0.4, 0.4, 0.4])
       p2 = system.part.add(id=1, pos=[1, 1, 1.0], v=[0.1, 0.1, 0.1])
       vel_diff = np.copy(system.velocity_difference(p1, p2))
       np.testing.assert_array_equal(
           vel_diff, p2.v - p1.v)
#
       system.part.clear()
#
   def disable_test_e_virt_sites_rotation(self):
       """A particle with virtual sites is plces on the boudary. We check if
          the forces yield the correct torque and if a rotation frequency is
          transmitted backl to the virtual sites."""
#
       system = self.system
       system.part.clear()
       system.min_global_cut = 2.5
       self.system.virtual_sites = VirtualSitesRelative(have_velocity=True)
#
       system.lees_edwards.protocol = \
           lees_edwards.LinearShear(
               shear_velocity=0.0,
               shear_direction=0,
               shear_plane_normal=1,
               initial_pos_offset=2.0
               )
#
       p1 = system.part.add(
           id=0, pos=[2.5, 5.0, 2.5], rotation=(1, 1, 1))
#
       p2 = system.part.add(
   pos=(
       2.5, 6.0, 2.5), ext_force=(
           1.0, 0., 0.), fix=(
               1, 1, 1))
       p2.vs_auto_relate_to(0)
       p3 = system.part.add(pos=(2.5, 4.0, 2.5),
                            ext_force=(-1.0, 0., 0.), fix=(1, 1, 1))
       p3.vs_auto_relate_to(0)
#
       system.integrator.run(0, recalc_forces=True)
#
       np.testing.assert_array_almost_equal(
           np.copy(p1.torque_lab), [0.0, 0.0, -2.0])
#
       p1.omega_lab = (0., 0., 2.5)
       system.integrator.run(0, recalc_forces=True)
#
       np.testing.assert_array_almost_equal(
           np.copy(p2.v), [-2.5, 0.0, 0.0])
       np.testing.assert_array_almost_equal(
           np.copy(p3.v), [2.5, -5.0, 0.0])
#
   def disable_test_f_virt_sites_interaction(self):
       """A virtual site interacts with a real particle via a DPD interaction to get
          a velocity dependent force. First we measure a force within the primary
          simulation box as reference. Then we compare it first with the situation of
          the interaction across the Lees Edwards boundary and second with the vector
          from the real particle to the virtual site crossing the boundary."""
#
       system = self.system
       system.time = 0.0
       system.part.clear()
       self.system.virtual_sites = VirtualSitesRelative(have_velocity=True)
#
       system.thermostat.set_dpd(kT=0.0, seed=1)
       system.non_bonded_inter[11, 11].dpd.set_params(
           weight_function=0, gamma=1.75, r_cut=2.,
           trans_weight_function=0, trans_gamma=1.5, trans_r_cut=2.0)
#
       system.lees_edwards.protocol = \
           lees_edwards.LinearShear(
               shear_velocity=2.0,
               shear_direction=0,
               shear_plane_normal=1,
               initial_pos_offset=0.0
               )
#
       p1 = system.part.add(
           id=0, pos=[2.5, 2.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
       p2 = system.part.add(pos=(2.5, 3.5, 2.5), type=11)
       p2.vs_auto_relate_to(0)
#
       p3 = system.part.add(pos=(2.5, 4.5, 2.5), type=11, v=(2.0, 1., 1.25))
#
       system.integrator.run(0, recalc_forces=True)
#
       f_p1 = np.copy(p1.f)
       f_p2 = np.copy(p2.f)
       f_p3 = np.copy(p3.f)
#
       system.part.clear()
#
       p1 = system.part.add(
           id=0, pos=[2.5, 3.75, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
       p2 = system.part.add(pos=(2.5, 4.75, 2.5), type=11)
       p2.vs_auto_relate_to(0)
#
       p3 = system.part.add(pos=(2.5, 5.75, 2.5), type=11, v=(0.0, 1., 1.25))
#
       system.integrator.run(0, recalc_forces=True)
#
       np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
       np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
       np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)
#
       system.part.clear()
#
       p1 = system.part.add(
           id=0, pos=[2.5, 4.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
       p2 = system.part.add(pos=(2.5, 5.5, 2.5), type=11)
       p2.vs_auto_relate_to(0)
#
       p3 = system.part.add(pos=(2.5, 6.5, 2.5), type=11, v=(0., 1., 1.25))
#
       system.integrator.run(0, recalc_forces=True)
#
       np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
       np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
       np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)
       
       system.part.clear()
#
       system.non_bonded_inter[11, 11].dpd.set_params(
           weight_function=0, gamma=0, r_cut=0,
           trans_weight_function=0, trans_gamma=0, trans_r_cut=0)
#
   def run_initial_pair_test(self, pos1, pos2, le_protocol):
       system = self.system
       system.part.clear()
       system.time = 0.0
       system.cell_system.set_n_square(use_verlet_lists=False)
       # Parameters
       sigma = 1.
       eps = 1E-7  # Weak force is intentiaonl, so particles aren't driven apart for a long time
       cut = sigma * 2**(1. / 6.)
#
       # box
       l = 10.0
#
       # Setup
       system.box_l = l, l, l
       system.time_step = 0.01
       system.thermostat.turn_off()
#
       system.lees_edwards.protocol = lees_edwards.Off() 
#
       p1 = system.part.add(pos=pos1)
       p2 = system.part.add(pos=pos2)
#
       # interactions
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
#
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
       return p1, p2
#
   const_offset_params = {
           'shear_velocity': 0.0,
           'shear_direction': 0,
           'shear_plane_normal': 1,
           'initial_pos_offset': 17.2}
#
   shear_params = {
           'shear_velocity': 0.1,
           'shear_direction': 0,
           'shear_plane_normal': 2,
           'initial_pos_offset': -np.sqrt(0.1)}
#
   def disable_test_z1_lj_interaction(self):
       """Takes two LJ particles across a boundary and checks with verify LJ forces
          if the interactin is right. Done in the columnar cell system and in n_square"""
#
       system = self.system
       for params in self.const_offset_params, self.shear_params:        
           p1, p2 = self.run_initial_pair_test(
               [1, 4, 5], [1, 4.99, 5], lees_edwards.LinearShear(**params))
           # Let them move at (nearly) constant distance
           p1.v = -1, 2, 3
           p2.v = p1.v
           for i in range(5000):
               system.integrator.run(5)
               assert np.linalg.norm(p1.f) > 0
               verify_lj_forces(system, 1E-10)
   
           #
           n_nodes = self.system.cell_system.get_state()['n_nodes']
           system.cell_system.node_grid = [1, 1, n_nodes]
           system.cell_system.set_domain_decomposition(
               fully_connected=[True, False, False])
   
           # Reset distance to 1
           p2.pos = p1.pos + np.array((0, 1, 0))
           # Let them move at (nearly) constant distance
           p1.v = -1, 2, 3
           p2.v = p1.v
           for i in range(5000):
               system.integrator.run(5)
               assert np.linalg.norm(p1.f) > 0
               verify_lj_forces(system, 1E-10)
   
       # Clean up 
       system.part.clear()
#
       # Turn off lj interaction
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=0, sigma=0, cutoff=0, shift=0)
   
   def disable_test_z2_lj_interaction(self):
       """Seets up moving particles across the LE boundary such that
         they will come into ia range during the simulation. Test that
         a force actuall occurs."""
#
       system = self.system
       for params in self.const_offset_params, self.shear_params:        
           print("Testing", params)
           dx = system.box_l[0] / 2
           p1, p2 = self.run_initial_pair_test(
               [0, 0.5, 0.5], [dx, 0.5, 0.5], lees_edwards.LinearShear(**params))
           # Let them move at (nearly) constant distance
           p1.v = [2, 0, 0]
           p2.v = [-.3, 0, 0]
#
           lj_cut = system.non_bonded_inter[p1.type, p2.type].lennard_jones.get_params()[
                                                                                       "cutoff"]
           distance = system.distance_vec(p1, p2)[0] - lj_cut
           vel_diff = system.velocity_difference(p1, p2)[0]
           encounter_time = np.abs(distance / vel_diff)
           print("expected encounter", encounter_time)
           before_steps = int(encounter_time / system.time_step - 1)
           for i in range(before_steps):
               system.integrator.run(1)
               np.testing.assert_equal(p1.f, [0, 0, 0])
               np.testing.assert_equal(p2.f, [0, 0, 0])
           system.integrator.run(2)
           print(system.distance_vec(p1, p2))
           self.assertGreater(np.linalg.norm(p1.f), 0)
       # Clean up 
       system.part.clear()
#
       # Turn off lj interaction
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=0, sigma=0, cutoff=0, shift=0)
   
   def disable_test_z3_moving_lj_particle_with_partner(self):
       """Same as before but this time with three interaction partners. One of them
          gets an artificial velocity with the same value as the shearing velocity."""
#
       system = self.system
       system.time = 0.0
       system.part.clear()
       system.cell_system.set_n_square(use_verlet_lists=False)
       # Parameters
       sigma = 1.
       eps = 1
       cut = sigma * 2**(1. / 6.)
#
       # box
       l = 10.0
#
       # Setup
       system.box_l = l, l, l
       system.time_step = 0.01
       system.thermostat.turn_off()
#
       system.lees_edwards.protocol = lees_edwards.Off() 
#
       system.part.add(pos=[2.5, 4.5, 2.5])
       system.part.add(pos=[1.5, 5.5, 2.5])
       system.part.add(pos=[1.5, 6.5, 2.5], v=[0.1, 0, 0])
#
       # interactions
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
#
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
#
       # Switch to constant offset protocol
       params = {
           'shear_velocity': 0.1,
           'shear_direction': 0,
           'shear_plane_normal': 1,
           'initial_pos_offset': 0.0}
       system.lees_edwards.protocol = lees_edwards.LinearShear(**params)
#
       for i in range(20):
           system.integrator.run(1, recalc_forces=True)
           verify_lj_forces(system, 1E-10)
#
       n_nodes = self.system.cell_system.get_state()['n_nodes']
       if n_nodes == 2:
           system.cell_system.node_grid = [1, 1, 2]
           system.cell_system.set_domain_decomposition(
               fully_connected=[True, False, False])
       
       if n_nodes == 4:
           system.cell_system.node_grid = [1, 2, 2]
           system.cell_system.set_domain_decomposition(
               fully_connected=[True, False, False])
#
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
#
       system.part.clear()
#
       # Turn off lj interaction
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=0, sigma=0, cutoff=0, shift=0)
#
   def setup_lj_liquid(self):
       system = self.system
       system.cell_system.set_n_square(use_verlet_lists=False)
       # Parameters
       n = 300 
       phi = 0.55
       sigma = 1.
       eps = 1
       cut = sigma * 2**(1. / 6.)
#
       # box
       l = (n / 6. * np.pi * sigma**3 / phi)**(1. / 3.)
#
       # Setup
       system.box_l = l, l, l
       system.part.clear()
#
       system.time_step = 0.01
       system.thermostat.turn_off()
#
       system.lees_edwards.protocol = lees_edwards.Off() 
       system.lees_edwards.pos_offset = 0
       system.lees_edwards.shear_velocity = 0
#
       system.part.add(pos=np.random.random((n, 3)) * l)
#
       # interactions
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
       # Remove overlap
       system.integrator.set_steepest_descent(
           f_max=0, gamma=0.05, max_displacement=0.05)
       while system.analysis.energy()["total"] > 8 * n:
           system.integrator.run(20)
           print(system.analysis.energy()["total"])
#
       system.integrator.set_vv()
       system.part[:].v = np.random.random((n, 3))
       for i in range(5):
           e_kin = 0.5 * np.sum(system.part[:].v**2)
           system.part[:].v = system.part[:].v / np.sqrt(e_kin)
           system.integrator.run(40)
       verify_lj_forces(system, 1E-10)
   
   def disable_test_z4_lj_fluid_constant_offset(self):
       """Simulates a static LJ liquid with a constant offset  and verifies forces.
          This is to make sure that the get_mi_works corrctly and now pairs get lost
          or are outdated in the short range loop. """
#
       system = self.system
       self.setup_lj_liquid()
       print(system.cell_system.get_state())
#
       # Switch to constant offset protocol
       system.lees_edwards.protocol = lees_edwards.LinearShear(**
           {
           'shear_velocity': 0.0,
           'shear_direction': 0,
           'shear_plane_normal': 1,
           'initial_pos_offset': -2 * system.box_l[0]})

       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
       print("Survived switch to const offset")
       
       n_nodes = self.system.cell_system.get_state()['n_nodes']
       if n_nodes == 2:
           system.cell_system.node_grid = [1, 1, 2]
       
       if n_nodes == 4:
           system.cell_system.node_grid = [1, 2, 2]
       system.cell_system.set_domain_decomposition(
               fully_connected=[True, False, False])
       
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
       print("Survived switch to columnar cs") 
       system.part.clear()
#
       # Turn off lj interaction
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=0, sigma=0, cutoff=0, shift=0)
#
   def disable_test_z5_lj(self):
       """Simulates an LJ liquid under linear shear and verifies forces. This is to make sure that no pairs
          get lost or are outdated in the short range loop.
          To have deterministic forces, velocity capping is used rather than a thermostat."""
       system = self.system
       self.setup_lj_liquid()
       system.time = 0
       params = {
           'shear_velocity': 0.1,
           'shear_direction': 2,
           'shear_plane_normal': 0,
           'initial_pos_offset': 0.}
       system.lees_edwards.protocol = lees_edwards.LinearShear(**params)
       system.integrator.run(0, recalc_forces=True)
       verify_lj_forces(system, 1E-10)
#
       system.part[:].v = np.random.random((len(system.part), 3))
       # Integrate
       for i in range(40):
           e_kin = 0.5 * np.sum(system.part[:].v**2)
           system.part[:].v = system.part[:].v / np.sqrt(e_kin)
           system.integrator.run(20)
       f1 = system.part[:].f
       p1 = system.part[:].pos_folded
#
       # Switch to constant offset protocol
       new_params = params.copy()
       new_offset = system.lees_edwards.pos_offset - \
           system.time_step * system.lees_edwards.shear_velocity
#
       new_params.update(shear_velocity=0,
           initial_pos_offset=new_offset)
       system.time = 0
       system.lees_edwards.protocol = lees_edwards.LinearShear(**new_params)
       system.integrator.run(0, recalc_forces=True)
       f2 = system.part[:].f
       np.testing.assert_allclose(f1, f2)
#
       p2 = system.part[:].pos_folded
       np.testing.assert_allclose(p1, p2)
#
       # Verify lj forces on the particles.
       verify_lj_forces(system, 1E-10)
       
       n_nodes = self.system.cell_system.get_state()['n_nodes']
       if n_nodes == 2:
           system.cell_system.node_grid = [1, 1, 2]
       
       if n_nodes == 4:
           system.cell_system.node_grid = [1, 2, 2]
       
       print(system.cell_system.get_state())
       system.cell_system.set_domain_decomposition(
           fully_connected=[True, False, False])
       
       system.integrator.run(0, recalc_forces=True)
       p3 = system.part[:].pos_folded
       np.testing.assert_allclose(p1, p3)
       
       f3 = system.part[:].f
       np.testing.assert_allclose(f1, f3)
#
       verify_lj_forces(system, 1E-10)
#
       system.part.clear()
#
       # Turn off lj interaction
       system.non_bonded_inter[0, 0].lennard_jones.set_params(
           epsilon=0, sigma=0, cutoff=0, shift=0)
#
if __name__ == "__main__":
    ut.main()
