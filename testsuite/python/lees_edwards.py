import espressomd
from espressomd.interactions import HarmonicBond
import espressomd.lees_edwards as lees_edwards
from espressomd.virtual_sites import VirtualSitesRelative, VirtualSitesOff
#from tests_common import verify_lj_forces

import unittest as ut
import numpy as np

params_lin = {'shear_velocity': 1.2, 'initial_pos_offset': 0.1, 'time_0': 0.1}
lin_protocol = lees_edwards.LinearShear(**params_lin)
params_osc = {'amplitude': 2.3, 'omega': 2.51, 'time_0': -2.1}


def axis(coord):
    """Return Cartesian axis from coordinate index"""
    res = np.zeros(3)
    res[coord] = 1
    return res


class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=True)

    time_step = 0.5
    system.time_step = time_step

    def test_00_is_off_by_default(self):

        system = self.system

        # Protocol should be off by default
        self.assertTrue(system.lees_edwards.protocol is None)
        self.assertEqual(
            system.lees_edwards.shear_velocity, 0)
        self.assertEqual(
            system.lees_edwards.pos_offset, 0)

    def test_protocols(self):
        """Tests shear velocity and pos offset vs protocol and time"""

        # Linear shear
        system = self.system
        system.time = 3.15
        system.lees_edwards.protocol = lin_protocol

        # check protocol assignment
        self.assertEqual(self.system.lees_edwards.protocol, lin_protocol)

        # check parameter setter/getter consistency
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params_lin)

        # Check pos offset and shear velocity
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

        # Oscillatory shear
        system.lees_edwards.protocol = lees_edwards.OscillatoryShear(
            **params_osc)

        # check parameter setter/getter consistency
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params_osc)

        # check pos offset and shear velocity at different times,
        # check that LE offsets are recalculated on simulation tmie change
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
        # Check that time change during integration updates offsets
        system.integrator.run(1)
        time = system.time
        expected_pos = params_osc['amplitude'] * \
            np.sin(params_osc['omega'] * (time - params_osc['time_0']))
        expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
            np.cos(params_osc['omega'] * (time - params_osc['time_0']))
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity,
            expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check that LE is diabled correctly
        system.lees_edwards.protocol = None
        self.assertEqual(
            system.lees_edwards.shear_velocity, 0)
        self.assertEqual(
            system.lees_edwards.pos_offset, 0)

    def test_boundary_crossing(self):
        """A particle crosses the upper and lower boundary to test if position
        and velocity are updated correctly."""
#
        system = self.system
        system.part.clear()
#
        # Set up a one particle system and check the position offset after crossing the boundary
        # Test for upper boundary
#
        dir = [0, 1, 2]
        system.lees_edwards.protocol = lin_protocol 
#
        for sheardir in dir:
            for shear_plane_normal in dir:
                if sheardir == shear_plane_normal: continue
                system.lees_edwards.shear_direction = sheardir
                system.lees_edwards.shear_plane_normal = shear_plane_normal

                shear_axis = axis(sheardir)
#
                pos = system.box_l - 0.01
                vel = np.array([1, 1, 1])
                p = system.part.add(pos=pos, v=vel)
#
                system.time = 0.0
                expected_pos = (pos + vel * system.time_step -
                                system.lees_edwards.pos_offset * shear_axis) 
                expected_delta_vel = -params_lin['shear_velocity'] * shear_axis
#
                system.integrator.run(1)
#
                np.testing.assert_almost_equal(p.v, vel + expected_delta_vel)
                np.testing.assert_almost_equal(p.pos, expected_pos)

                # Test the crossing lower boundary
                p.v = -vel
                expected_pos = p.pos + p.v * system.time_step + \
                    system.lees_edwards.pos_offset * shear_axis
                system.integrator.run(1)
                np.testing.assert_almost_equal(p.pos, expected_pos)
                np.testing.assert_almost_equal(p.v, -vel - expected_delta_vel)

    def test_distance_vel_diff(self):
        """check distance and velocity difference calculation across LE boundary
        """
        \
        system = self.system
        system.time = 0.0
        system.part.clear()
        system.lees_edwards.protocol = lin_protocol
        epsilon = 0.01
#
        dir = [0, 1, 2]
#
        for sheardir in dir:
            for shear_plane_normal in dir:
                if sheardir != shear_plane_normal:
                    #
                    system.lees_edwards.shear_direction = sheardir
                    system.lees_edwards.shear_plane_normal = shear_plane_normal
                    shear_axis = axis(sheardir)

                    system.lees_edwards.protocol = lin_protocol
                    p1 = system.part.add(
                        pos=[epsilon] * 3, v=np.random.random(3), fix=[1] * 3)
                    p2 = system.part.add(
                        pos=system.box_l - epsilon, v=np.random.random(3), fix=[1] * 3)
                    r_euclid = -2 * np.array([epsilon] * 3)

#
                    # check distance 
                    np.testing.assert_allclose(
                        system.distance_vec(p1, p2), r_euclid + system.lees_edwards.pos_offset * shear_axis)
                    np.testing.assert_allclose(
                        np.copy(system.distance_vec(p1, p2)),
                        -np.copy(system.distance_vec(p2, p1)))

                    # Ceck velocity difference
                    np.testing.assert_allclose(
                        system.velocity_difference(p1, p2),
                        p2.v - p1.v - system.lees_edwards.shear_velocity * shear_axis)

    def test_interactions(self):
        """We place two particles crossing a boundary and connect them with an unbonded
           and bonded interaction. We test with the resulting stress tensor if the offset
           is included properly"""
#
        system = self.system
        system.time = 0.0
        system.part.clear()
        system.lees_edwards.protocol = lin_protocol
        epsilon = 0.01
#
        dir = [0, 1, 2]
#
        for sheardir in dir:
            for shear_plane_normal in dir:
                if sheardir != shear_plane_normal:
                    #
                    system.lees_edwards.shear_direction = sheardir
                    system.lees_edwards.shear_plane_normal = shear_plane_normal

                    system.lees_edwards.protocol = lin_protocol
                    p1 = system.part.add(pos=[epsilon] * 3, fix=[1] * 3)
                    p2 = system.part.add(
                        pos=system.box_l - epsilon, fix=[1] * 3)

#
                    # check bonded interaction
                    k_bond = 2.1
                    r_cut = 1.5
                    harm = HarmonicBond(k=k_bond, r_0=0.0)
                    system.bonded_inter.add(harm)
                    p1.add_bond((harm, p2))

                    r_12 = system.distance_vec(p1, p2)
                    system.integrator.run(0)

                    np.testing.assert_allclose(
                        k_bond * r_12,
                        np.copy(p1.f))
                    np.testing.assert_allclose(np.copy(p1.f), -np.copy(p2.f))

                    np.testing.assert_allclose(
                        system.analysis.pressure_tensor()["bonded"],
                        np.outer(r_12, p2.f) / system.volume())

                    np.testing.assert_almost_equal(
                        system.analysis.energy()["bonded"],
                        1 / 2 * k_bond * system.distance(p1, p2)**2)
                    p1.bonds = [] 

                    # Check non bonded interaction
                    k_non_bonded = 3.2
                    # NOTE: The force is k*n *distance, hence the 1/2
                    system.non_bonded_inter[0, 0].soft_sphere.set_params(
                        a=k_non_bonded / 2, n=-2, cutoff=r_cut)
                    system.integrator.run(0)
                    r_12 = system.distance_vec(p1, p2)

                    np.testing.assert_allclose(
                        (k_non_bonded) * r_12,
                        np.copy(p1.f))
                    np.testing.assert_allclose(np.copy(p1.f), -np.copy(p2.f))

                    np.testing.assert_allclose(
                        system.analysis.pressure_tensor()["non_bonded"],
                        np.outer(r_12, p2.f) / system.volume())

                    np.testing.assert_almost_equal(
                        system.analysis.energy()["non_bonded"],
                        1 / 2 * k_non_bonded * system.distance(p1, p2)**2)

                    system.non_bonded_inter[0, 0].soft_sphere.set_params(
                        a=0, n=-2, cutoff=r_cut)
                    system.part.clear()

    def test_virt_sites_rotation(self):
        """A particle with virtual sites is plces on the boudary. We check if
           the forces yield the correct torque and if a rotation frequency is
           transmitted backl to the virtual sites."""

        system = self.system
        system.part.clear()
        system.min_global_cut = 2.5
        self.system.virtual_sites = VirtualSitesRelative()
        system.lees_edwards.protocol = lin_protocol 
        p1 = system.part.add(
            id=0, pos=[2.5, 5.0, 2.5], rotation=(1, 1, 1))

        p2 = system.part.add(
            pos=(
                2.5, 6.0, 2.5), ext_force=(
                1.0, 0., 0.), fix=(
                1, 1, 1))
        p2.vs_auto_relate_to(0)
        p3 = system.part.add(pos=(2.5, 4.0, 2.5),
                             ext_force=(-1.0, 0., 0.), fix=(1, 1, 1))
        p3.vs_auto_relate_to(0)

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(
            np.copy(p1.torque_lab), [0.0, 0.0, -2.0])

        p1.omega_lab = (0., 0., 2.5)
        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(
            np.copy(p2.v), [-2.5, 0.0, 0.0])
        np.testing.assert_array_almost_equal(
            np.copy(p3.v), [2.5, -5.0, 0.0])
        system.virtual_sites = VirtualSitesOff()

    def test_virt_sites_interaction(self):
        """A virtual site interacts with a real particle via a DPD interaction to get
           a velocity dependent force. First we measure a force within the primary
           simulation box as reference. Then we compare it first with the situation of
           the interaction across the Lees Edwards boundary and second with the vector
           from the real particle to the virtual site crossing the boundary."""

        system = self.system
        system.time = 0.0
        system.part.clear()
        self.system.virtual_sites = VirtualSitesRelative(have_velocity=True)

        system.thermostat.set_dpd(kT=0.0, seed=1)
        system.non_bonded_inter[11, 11].dpd.set_params(
            weight_function=0, gamma=1.75, r_cut=2.,
            trans_weight_function=0, trans_gamma=1.5, trans_r_cut=2.0)

        system.lees_edwards.protocol = \
            lees_edwards.LinearShear(
                shear_velocity=2.0,
                shear_direction=0,
                shear_plane_normal=1,
                initial_pos_offset=0.0
            )

        p1 = system.part.add(
            id=0, pos=[2.5, 2.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 3.5, 2.5), type=11)
        p2.vs_auto_relate_to(0)

        p3 = system.part.add(pos=(2.5, 4.5, 2.5), type=11, v=(2.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        f_p1 = np.copy(p1.f)
        f_p2 = np.copy(p2.f)
        f_p3 = np.copy(p3.f)

        system.part.clear()

        p1 = system.part.add(
            id=0, pos=[2.5, 3.75, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 4.75, 2.5), type=11)
        p2.vs_auto_relate_to(0)

        p3 = system.part.add(pos=(2.5, 5.75, 2.5), type=11, v=(0.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
        np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
        np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)

        system.part.clear()

        p1 = system.part.add(
            id=0, pos=[2.5, 4.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 5.5, 2.5), type=11)
        p2.vs_auto_relate_to(0)

        p3 = system.part.add(pos=(2.5, 6.5, 2.5), type=11, v=(0., 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
        np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
        np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)

        system.part.clear()

        system.non_bonded_inter[11, 11].dpd.set_params(
            weight_function=0, gamma=0, r_cut=0,
            trans_weight_function=0, trans_gamma=0, trans_r_cut=0)
        system.virtual_sites = VirtualSitesOff()


#     def run_initial_pair_test(self, pos1, pos2, le_protocol):
#         system = self.system
#         system.part.clear()
#         system.time = 0.0
#         system.cell_system.set_n_square(use_verlet_lists=False)
#         # Parameters
#         sigma = 1.
#         eps = 1E-7  # Weak force is intentiaonl, so particles aren't driven apart for a long time
#         cut = sigma * 2**(1. / 6.)
# #
#         # box
#         l = 10.0
# #
#         # Setup
#         system.box_l = l, l, l
#         system.time_step = 0.01
#         system.thermostat.turn_off()
# #
#         system.lees_edwards.protocol = lees_edwards.Off() 
# #
#         p1 = system.part.add(pos=pos1)
#         p2 = system.part.add(pos=pos2)
# #
#         # interactions
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
# #
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
#         return p1, p2
# #
#     const_offset_params = {
#         'shear_velocity': 0.0,
#         'shear_direction': 0,
#         'shear_plane_normal': 1,
#         'initial_pos_offset': 17.2}
# #
#     shear_params = {
#         'shear_velocity': 0.1,
#         'shear_direction': 0,
#         'shear_plane_normal': 2,
#         'initial_pos_offset': -np.sqrt(0.1)}
# #
# 
#     def disable_test_z1_lj_interaction(self):
#         """Takes two LJ particles across a boundary and checks with verify LJ forces
#            if the interactin is right. Done in the columnar cell system and in n_square"""
# #
#         system = self.system
#         for params in self.const_offset_params, self.shear_params:        
#             p1, p2 = self.run_initial_pair_test(
#                 [1, 4, 5], [1, 4.99, 5], lees_edwards.LinearShear(**params))
#             # Let them move at (nearly) constant distance
#             p1.v = -1, 2, 3
#             p2.v = p1.v
#             for i in range(5000):
#                 system.integrator.run(5)
#                 assert np.linalg.norm(p1.f) > 0
#                 verify_lj_forces(system, 1E-10)
# 
#             #
#             n_nodes = self.system.cell_system.get_state()['n_nodes']
#             system.cell_system.node_grid = [1, 1, n_nodes]
#             system.cell_system.set_domain_decomposition(
#                 fully_connected=[True, False, False])
# 
#             # Reset distance to 1
#             p2.pos = p1.pos + np.array((0, 1, 0))
#             # Let them move at (nearly) constant distance
#             p1.v = -1, 2, 3
#             p2.v = p1.v
#             for i in range(5000):
#                 system.integrator.run(5)
#                 assert np.linalg.norm(p1.f) > 0
#                 verify_lj_forces(system, 1E-10)
# 
#         # Clean up 
#         system.part.clear()
# #
#         # Turn off lj interaction
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=0, sigma=0, cutoff=0, shift=0)
# 
#     def disable_test_z2_lj_interaction(self):
#         """Seets up moving particles across the LE boundary such that
#           they will come into ia range during the simulation. Test that
#           a force actuall occurs."""
# #
#         system = self.system
#         for params in self.const_offset_params, self.shear_params:        
#             print("Testing", params)
#             dx = system.box_l[0] / 2
#             p1, p2 = self.run_initial_pair_test(
#                 [0, 0.5, 0.5], [dx, 0.5, 0.5], lees_edwards.LinearShear(**params))
#             # Let them move at (nearly) constant distance
#             p1.v = [2, 0, 0]
#             p2.v = [-.3, 0, 0]
# #
#             lj_cut = system.non_bonded_inter[p1.type, p2.type].lennard_jones.get_params()[
#                 "cutoff"]
#             distance = system.distance_vec(p1, p2)[0] - lj_cut
#             vel_diff = system.velocity_difference(p1, p2)[0]
#             encounter_time = np.abs(distance / vel_diff)
#             print("expected encounter", encounter_time)
#             before_steps = int(encounter_time / system.time_step - 1)
#             for i in range(before_steps):
#                 system.integrator.run(1)
#                 np.testing.assert_equal(p1.f, [0, 0, 0])
#                 np.testing.assert_equal(p2.f, [0, 0, 0])
#             system.integrator.run(2)
#             print(system.distance_vec(p1, p2))
#             self.assertGreater(np.linalg.norm(p1.f), 0)
#         # Clean up 
#         system.part.clear()
# #
#         # Turn off lj interaction
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=0, sigma=0, cutoff=0, shift=0)
# 
#     def disable_test_z3_moving_lj_particle_with_partner(self):
#         """Same as before but this time with three interaction partners. One of them
#            gets an artificial velocity with the same value as the shearing velocity."""
# #
#         system = self.system
#         system.time = 0.0
#         system.part.clear()
#         system.cell_system.set_n_square(use_verlet_lists=False)
#         # Parameters
#         sigma = 1.
#         eps = 1
#         cut = sigma * 2**(1. / 6.)
# #
#         # box
#         l = 10.0
# #
#         # Setup
#         system.box_l = l, l, l
#         system.time_step = 0.01
#         system.thermostat.turn_off()
# #
#         system.lees_edwards.protocol = lees_edwards.Off() 
# #
#         system.part.add(pos=[2.5, 4.5, 2.5])
#         system.part.add(pos=[1.5, 5.5, 2.5])
#         system.part.add(pos=[1.5, 6.5, 2.5], v=[0.1, 0, 0])
# #
#         # interactions
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
# #
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
# #
#         # Switch to constant offset protocol
#         params = {
#             'shear_velocity': 0.1,
#             'shear_direction': 0,
#             'shear_plane_normal': 1,
#             'initial_pos_offset': 0.0}
#         system.lees_edwards.protocol = lees_edwards.LinearShear(**params)
# #
#         for i in range(20):
#             system.integrator.run(1, recalc_forces=True)
#             verify_lj_forces(system, 1E-10)
# #
#         n_nodes = self.system.cell_system.get_state()['n_nodes']
#         if n_nodes == 2:
#             system.cell_system.node_grid = [1, 1, 2]
#             system.cell_system.set_domain_decomposition(
#                 fully_connected=[True, False, False])
# 
#         if n_nodes == 4:
#             system.cell_system.node_grid = [1, 2, 2]
#             system.cell_system.set_domain_decomposition(
#                 fully_connected=[True, False, False])
# #
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
# #
#         system.part.clear()
# #
#         # Turn off lj interaction
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=0, sigma=0, cutoff=0, shift=0)
# #
# 
#     def setup_lj_liquid(self):
#         system = self.system
#         system.cell_system.set_n_square(use_verlet_lists=False)
#         # Parameters
#         n = 300 
#         phi = 0.55
#         sigma = 1.
#         eps = 1
#         cut = sigma * 2**(1. / 6.)
# #
#         # box
#         l = (n / 6. * np.pi * sigma**3 / phi)**(1. / 3.)
# #
#         # Setup
#         system.box_l = l, l, l
#         system.part.clear()
# #
#         system.time_step = 0.01
#         system.thermostat.turn_off()
# #
#         system.lees_edwards.protocol = lees_edwards.Off() 
#         system.lees_edwards.pos_offset = 0
#         system.lees_edwards.shear_velocity = 0
# #
#         system.part.add(pos=np.random.random((n, 3)) * l)
# #
#         # interactions
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
#         # Remove overlap
#         system.integrator.set_steepest_descent(
#             f_max=0, gamma=0.05, max_displacement=0.05)
#         while system.analysis.energy()["total"] > 8 * n:
#             system.integrator.run(20)
#             print(system.analysis.energy()["total"])
# #
#         system.integrator.set_vv()
#         system.part[:].v = np.random.random((n, 3))
#         for i in range(5):
#             e_kin = 0.5 * np.sum(system.part[:].v**2)
#             system.part[:].v = system.part[:].v / np.sqrt(e_kin)
#             system.integrator.run(40)
#         verify_lj_forces(system, 1E-10)
# 
#     def disable_test_z4_lj_fluid_constant_offset(self):
#         """Simulates a static LJ liquid with a constant offset  and verifies forces.
#            This is to make sure that the get_mi_works corrctly and now pairs get lost
#            or are outdated in the short range loop. """
# #
#         system = self.system
#         self.setup_lj_liquid()
#         print(system.cell_system.get_state())
# #
#         # Switch to constant offset protocol
#         system.lees_edwards.protocol = lees_edwards.LinearShear(**
#                                                                 {
#                                                                     'shear_velocity': 0.0,
#                                                                     'shear_direction': 0,
#                                                                     'shear_plane_normal': 1,
#                                                                     'initial_pos_offset': -2 * system.box_l[0]})
# 
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
#         print("Survived switch to const offset")
# 
#         n_nodes = self.system.cell_system.get_state()['n_nodes']
#         if n_nodes == 2:
#             system.cell_system.node_grid = [1, 1, 2]
# 
#         if n_nodes == 4:
#             system.cell_system.node_grid = [1, 2, 2]
#         system.cell_system.set_domain_decomposition(
#             fully_connected=[True, False, False])
# 
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
#         print("Survived switch to columnar cs") 
#         system.part.clear()
# #
#         # Turn off lj interaction
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=0, sigma=0, cutoff=0, shift=0)
# #
# 
#     def disable_test_z5_lj(self):
#         """Simulates an LJ liquid under linear shear and verifies forces. This is to make sure that no pairs
#            get lost or are outdated in the short range loop.
#            To have deterministic forces, velocity capping is used rather than a thermostat."""
#         system = self.system
#         self.setup_lj_liquid()
#         system.time = 0
#         params = {
#             'shear_velocity': 0.1,
#             'shear_direction': 2,
#             'shear_plane_normal': 0,
#             'initial_pos_offset': 0.}
#         system.lees_edwards.protocol = lees_edwards.LinearShear(**params)
#         system.integrator.run(0, recalc_forces=True)
#         verify_lj_forces(system, 1E-10)
# #
#         system.part[:].v = np.random.random((len(system.part), 3))
#         # Integrate
#         for i in range(40):
#             e_kin = 0.5 * np.sum(system.part[:].v**2)
#             system.part[:].v = system.part[:].v / np.sqrt(e_kin)
#             system.integrator.run(20)
#         f1 = system.part[:].f
#         p1 = system.part[:].pos_folded
# #
#         # Switch to constant offset protocol
#         new_params = params.copy()
#         new_offset = system.lees_edwards.pos_offset - \
#             system.time_step * system.lees_edwards.shear_velocity
# #
#         new_params.update(shear_velocity=0,
#                           initial_pos_offset=new_offset)
#         system.time = 0
#         system.lees_edwards.protocol = lees_edwards.LinearShear(**new_params)
#         system.integrator.run(0, recalc_forces=True)
#         f2 = system.part[:].f
#         np.testing.assert_allclose(f1, f2)
# #
#         p2 = system.part[:].pos_folded
#         np.testing.assert_allclose(p1, p2)
# #
#         # Verify lj forces on the particles.
#         verify_lj_forces(system, 1E-10)
# 
#         n_nodes = self.system.cell_system.get_state()['n_nodes']
#         if n_nodes == 2:
#             system.cell_system.node_grid = [1, 1, 2]
# 
#         if n_nodes == 4:
#             system.cell_system.node_grid = [1, 2, 2]
# 
#         print(system.cell_system.get_state())
#         system.cell_system.set_domain_decomposition(
#             fully_connected=[True, False, False])
# 
#         system.integrator.run(0, recalc_forces=True)
#         p3 = system.part[:].pos_folded
#         np.testing.assert_allclose(p1, p3)
# 
#         f3 = system.part[:].f
#         np.testing.assert_allclose(f1, f3)
# #
#         verify_lj_forces(system, 1E-10)
# #
#         system.part.clear()
# #
#         # Turn off lj interaction
#         system.non_bonded_inter[0, 0].lennard_jones.set_params(
#             epsilon=0, sigma=0, cutoff=0, shift=0)
# 

#
if __name__ == "__main__":
    ut.main()
