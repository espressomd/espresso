# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import random, seed
import espressomd
import math


@ut.skipIf(not espressomd.has_features(["MASS",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "LANGEVIN_PER_PARTICLE"]),
           "Features not available, skipping test!")
class ThermoTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    system.cell_system.skin = 5.0

    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))
    gamma_global_rot = np.zeros((3))

    # Particle properties
    mass = 0.0
    J = 0.0, 0.0, 0.0

    # Per-particle type parameters.
    # 2 different langevin parameters for particles.
    kT_p = np.zeros((2))
    # gamma_tran/gamma_rot matrix: [2 kinds of particles] x [3 dimensions X Y Z]
    # These matrices are assigning per-particle in corresponding test cases.
    gamma_tran_p = np.zeros((2, 3))
    gamma_rot_p = np.zeros((2, 3))

    # These variables will take the values to compare with.
    # Depending on the test case following matrices either equals to the previous
    # or the global corresponding parameters. The corresponding setting effect is an essence of
    # all the test cases' differentiation here.
    halfkT_p_validate = np.zeros((2))
    gamma_tran_p_validate = np.zeros((2, 3))
    gamma_rot_p_validate = np.zeros((2, 3))
    # Diffusivity
    D_tran_p_validate = np.zeros((2, 3))

    @classmethod
    def setUpClass(cls):
        np.random.seed(15)

    def setUp(self):
        self.system.time = 0.0
        self.system.part.clear()
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.system.thermostat.turn_off()
    
    def set_initial_cond(self):
        """
        Set all the particles to zero coordinates and velocities; same for time.
        The quaternion is set to default value.

        """
        system = self.system
        system.time = 0.0
        system.part[:].pos = np.zeros((3))
        system.part[:].v = np.zeros((3))
        system.part[:].omega_body = np.zeros((3))
        system.part[:].quat = np.array([1., 0., 0., 0.])

    def generate_scalar_ranged_rnd(self, min_par, max_par):
        """
        Generate the scaled random scalar in the range between
        min_par*max_par and (min_par+1.0)*max_par.

        Parameters
        ----------
        min_par : :obj:`int`
                  Minimal value parameter.
        max_par : :obj:`int`
                  Maximal value parameter.

        """
        res = (min_par + np.random.random()) * max_par
        return res

    def generate_vec_ranged_rnd(self, min_par, max_par):
        """
        Generate the scaled random 3D vector with a magnitude
        in the range between sqrt(3)*min_par*max_par and
        sqrt(3)*(min_par+1.0)*max_par.

        Parameters
        ----------
        min_par : :obj:`int`
                  Minimal value parameter.
        max_par : :obj:`int`
                  Maximal value parameter.

        """
        res = (min_par + np.random.random(3)) * max_par
        return res

    def set_langevin_global_defaults(self):
        """
        Setup the expected NVT thermostat viscous friction parameters.

        """

        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k,:] = self.gamma_global[:]
            self.gamma_rot_p_validate[k,:] = self.gamma_global[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def set_langevin_global_defaults_rot_differ(self):
        """
        Setup the expected NVT thermostat viscous friction parameters
        with a rotation-specific gamma.

        """
        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k,:] = self.gamma_global[:]
            self.gamma_rot_p_validate[k,:] = self.gamma_global_rot[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def dissipation_param_setup(self):
        """
        Setup the parameters for the following dissipation
        test.

        """

        system = self.system
        # Time
        self.system.time_step = 0.007

        # Space
        box = 1.0
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        self.kT = 0.0
        # The translational gamma isotropy is required here.
        # Global gamma for tests without particle-specific gammas.
        #
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the gamma less than
        # some maximal value according to the value max_gamma_param.
        # Also, it cannot be very small (min_gamma_param), otherwise the thermalization will require
        # too much of the CPU time. Same: for all such gamma assignments throughout the test.
        #
        min_gamma_param = 0.5
        max_gamma_param = 2.0 / 3.0
        gamma_rnd = self.generate_scalar_ranged_rnd(
            min_gamma_param, max_gamma_param)
        self.gamma_global = gamma_rnd, gamma_rnd, gamma_rnd
        # Additional test case for the specific global rotational gamma set.
        self.gamma_global_rot = self.generate_vec_ranged_rnd(0.5, 2.0 / 3.0)
        # Per-paricle values:
        self.kT_p = 0.0, 0.0
        # Either translational friction isotropy is required
        # or both translational and rotational ones.
        # Otherwise these types of motion will interfere.
        # ..Let's test both cases depending on the particle index.
        self.gamma_tran_p[0, 0] = self.generate_scalar_ranged_rnd(0.5, 1.0)
        self.gamma_tran_p[0, 1] = self.gamma_tran_p[0, 0]
        self.gamma_tran_p[0, 2] = self.gamma_tran_p[0, 0]
        self.gamma_rot_p[0,:] = self.generate_vec_ranged_rnd(0.5, 2.0 / 3.0)
        self.gamma_tran_p[1, 0] = self.generate_scalar_ranged_rnd(0.5, 1.0)
        self.gamma_tran_p[1, 1] = self.gamma_tran_p[1, 0]
        self.gamma_tran_p[1, 2] = self.gamma_tran_p[1, 0]
        self.gamma_rot_p[1, 0] = self.generate_scalar_ranged_rnd(
            0.5, 2.0 / 3.0)
        self.gamma_rot_p[1, 1] = self.gamma_rot_p[1, 0]
        self.gamma_rot_p[1, 2] = self.gamma_rot_p[1, 0]

        # Particles
        self.mass = 12.74
        self.J = 10.0, 10.0, 10.0
        for i in range(2):
            system.part.add(rotation=(1, 1, 1), pos=(0.0, 0.0, 0.0), id=i)
            system.part[i].v = 1.0, 1.0, 1.0
            if "ROTATION" in espressomd.features():
                system.part[i].omega_body = 1.0, 1.0, 1.0
            system.part[i].mass = self.mass
            system.part[i].rinertia = self.J

    def dissipation_viscous_drag_setup_bd(self):
        """
        Setup the specific parameters for the following dissipation
        test of the viscous drag terminal velocity stationarity.
        It is used by the BD test cases only for the moment.

        """
        ## Time
        # Large time_step is OK for the BD by its definition & its benefits
        self.system.time_step = 17.0
        ## NVT thermostat
        # Isotropic reassignment is required here for the drag tests
        self.gamma_global_rot = np.zeros((3))
        self.gamma_global_rot[0] = (0.5 + np.random.random()) * 2.0 / 3.0
        self.gamma_global_rot[1] = self.gamma_global_rot[0]
        self.gamma_global_rot[2] = self.gamma_global_rot[0]
        # Isotropy is required here for the drag tests
        self.gamma_rot_p[0, 0] = self.generate_scalar_ranged_rnd(0.5, 2.0/3.0)
        self.gamma_rot_p[0, 1] = self.gamma_rot_p[0, 0]
        self.gamma_rot_p[0, 2] = self.gamma_rot_p[0, 0]

    def fluctuation_dissipation_param_setup(self, n):
        """
        Setup the parameters for the following fluctuation-dissipation
        test.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """
        # Time
        self.system.time_step = 0.03

        # Space
        box = 10.0
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = self.generate_scalar_ranged_rnd(0.3, 5)
        # See the above comment regarding the gamma assignments.
        # Note: here & hereinafter specific variations in these ranges are related to
        # the test execution duration to achieve the required statistical
        # averages faster.
        self.gamma_global = self.generate_vec_ranged_rnd(1.5, 2.0 / 3.0)
        self.gamma_global_rot = self.generate_vec_ranged_rnd(0.2, 20)
        # Per-particle parameters
        self.kT_p = 2.5, 2.0
        for k in range(2):
            self.gamma_tran_p[k,:] = self.generate_vec_ranged_rnd(0.4, 10.0)
            self.gamma_rot_p[k,:] = self.generate_vec_ranged_rnd(0.2, 20.0)

        # Particles
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the mass higher than
        # some minimal value according to the value min_mass_param.
        # Also, it is expected to test the large enough mass (max_mass_param).
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time.
        min_mass_param = 0.2
        max_mass_param = 7.0
        self.mass = self.generate_scalar_ranged_rnd(
            min_mass_param, max_mass_param)
        self.J = self.generate_vec_ranged_rnd(min_mass_param, max_mass_param)
        for i in range(n):
            for k in range(2):
                ind = i + k * n
                part_pos = np.random.random(3) * box
                part_v = 0.0, 0.0, 0.0
                part_omega_body = 0.0, 0.0, 0.0
                self.system.part.add(
                    rotation=(
                        1,
                        1,
                        1),
                    id=ind,
                    mass=self.mass,
                    rinertia=self.J,
                    pos=part_pos,
                    v=part_v)
                if "ROTATION" in espressomd.features():
                    self.system.part[ind].omega_body = part_omega_body

    def check_dissipation(self):
        """
        Check the dissipation relations: the simple viscous decelleration test.

        """

        system = self.system
        tol = 1.25E-4
        for i in range(100):
            for k in range(2):
                for j in range(3):
                    # Note: velocity is defined in the lab frame of reference
                    # while gamma_tr is defined in the body one.
                    # Hence, only isotropic gamma_tran_p_validate could be
                    # tested here.
                    self.assertLess(abs(
                        system.part[k].v[j] - math.exp(- self.gamma_tran_p_validate[k, j] * system.time / self.mass)), tol)
                    if "ROTATION" in espressomd.features():
                        self.assertLess(abs(
                            system.part[k].omega_body[j] - math.exp(- self.gamma_rot_p_validate[k, j] * system.time / self.J[j])), tol)

    # Note: the decelleration test is needed for the Langevin thermostat only. Brownian thermostat is defined
    # over a larger time-step by its concept.
    def check_dissipation_viscous_drag(self):
        """
        Check the dissipation relations: the drag terminal velocity tests,
        aka the drift in case of the electrostatics

        """
        system = self.system
        tol = 7E-3
        if "EXTERNAL_FORCES" in espressomd.features():
            # Just some random forces
            f0 = -1.2, 58.3578, 0.002
            f1 = -15.112, -2.0, 368.0
            system.part[0].ext_force = f0
            system.part[1].ext_force = f1
            if "ROTATION" in espressomd.features():
                # Just some random torques
                tor0 = 12, 0.022, 87
                tor1 = -0.03, -174, 368
                system.part[0].ext_torque = tor0
                system.part[1].ext_torque = tor1
                # Let's set the dipole perpendicular to the torque
                if "DIPOLES" in espressomd.features():
                    dip0 = 0.0, tor0[2], -tor0[1]
                    dip1 = -tor1[2], 0.0, tor1[0]
                    system.part[0].dip = dip0
                    system.part[1].dip = dip1
                    tmp_axis0 = np.cross(tor0, dip0) / (np.linalg.norm(tor0) * np.linalg.norm(dip0))
                    tmp_axis1 = np.cross(tor1, dip1) / (np.linalg.norm(tor1) * np.linalg.norm(dip1))
            # Small number of steps is enough for the terminal velocity within the BD by its definition.
            # A simulation of the original saturation of the velocity.
            system.integrator.run(7)
            system.time = 0.0
            for k in range(2):
                system.part[k].pos = np.zeros((3))
            if "DIPOLES" in espressomd.features():
                    system.part[0].dip = dip0
                    system.part[1].dip = dip1
            for i in range(3):
                # Small number of steps
                system.integrator.run(2)
                for k in range(3):
                    # Eq. (14.34) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010)
                    # First (deterministic) term of the eq. (14.34) of Schlick2010 taking into account eq. (14.35).
                    self.assertLess(
                        abs(system.part[0].v[k] - f0[k] / self.gamma_tran_p_validate[0, k]), tol)
                    self.assertLess(
                        abs(system.part[1].v[k] - f1[k] / self.gamma_tran_p_validate[1, k]), tol)
                    # Second (deterministic) term of the Eq. (14.39) of Schlick2010.
                    self.assertLess(
                        abs(system.part[0].pos[k] - system.time * f0[k] / self.gamma_tran_p_validate[0, k]), tol)
                    self.assertLess(
                        abs(system.part[1].pos[k] - system.time * f1[k] / self.gamma_tran_p_validate[1, k]), tol)
                    # Same, a rotational analogy.
                    if "ROTATION" in espressomd.features():
                        self.assertLess(abs(
                            system.part[0].omega_lab[k] - tor0[k] / self.gamma_rot_p_validate[0, k]), tol)
                        self.assertLess(abs(
                            system.part[1].omega_lab[k] - tor1[k] / self.gamma_rot_p_validate[1, k]), tol)
                if "ROTATION" in espressomd.features() and "DIPOLES" in espressomd.features():
                    # Same, a rotational analogy. One is implemented using a simple linear algebra;
                    # the polar angles with a sign control just for a correct inverse trigonometric functions application.
                    cos_alpha0 = np.dot(dip0, system.part[0].dip) / (np.linalg.norm(dip0) * system.part[0].dipm)
                    cos_alpha0_test = np.cos(system.time * np.linalg.norm(tor0) / self.gamma_rot_p_validate[0, 0])
                    sgn0 = np.sign(np.dot(system.part[0].dip, tmp_axis0))
                    sgn0_test = np.sign(np.sin(system.time * np.linalg.norm(tor0) / self.gamma_rot_p_validate[0, 0]))
                    
                    cos_alpha1 = np.dot(dip1, system.part[1].dip) / (np.linalg.norm(dip1) * system.part[1].dipm)
                    cos_alpha1_test = np.cos(system.time * np.linalg.norm(tor1) / self.gamma_rot_p_validate[1, 0])
                    sgn1 = np.sign(np.dot(system.part[1].dip, tmp_axis1))
                    sgn1_test = np.sign(np.sin(system.time * np.linalg.norm(tor1) / self.gamma_rot_p_validate[1, 0]))
                    
                    self.assertLess(abs(cos_alpha0 - cos_alpha0_test), tol)
                    self.assertLess(abs(cos_alpha1 - cos_alpha1_test), tol)
                    self.assertEqual(sgn0, sgn0_test)
                    self.assertEqual(sgn1, sgn1_test)

    def check_fluctuation_dissipation(self, n, therm_steps, loops):
        """
        Check the fluctuation-dissipation relations: thermalization
        and diffusion properties.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.
        therm_steps : :obj:`int`
            Number of thermalization steps.
        loops : :obj:`int`
            Number of iterations in the sampling.

        """
        system = self.system
        # The thermalization and diffusion test
        # Checks if every degree of freedom has 1/2 kT of energy, even when
        # mass and inertia tensor are active
        # Check the factual translational diffusion.
        #
        # matrices: [2 types of particless] x [3 dimensions X Y Z]
        # velocity^2, omega^2, position^2
        v2 = np.zeros((2, 3))
        o2 = np.zeros((2, 3))
        dr2 = np.zeros((2, 3))
        # Variance to compare with:
        sigma2_tr = np.zeros((2))
        # Comparable variance:
        dr_norm = np.zeros((2))

        pos0 = np.zeros((2 * n, 3))
        for p in range(n):
            for k in range(2):
                ind = p + k * n
                pos0[ind,:] = system.part[ind].pos
        dt0 = self.mass / self.gamma_tran_p_validate

        system.integrator.run(therm_steps)

        int_steps = 5
        for i in range(loops):
            system.integrator.run(int_steps)
            # Get kinetic energy in each degree of freedom for all particles
            for p in range(n):
                for k in range(2):
                    ind = p + k * n
                    v = system.part[ind].v
                    if "ROTATION" in espressomd.features():
                        o = system.part[ind].omega_body
                        o2[k,:] = o2[k,:] + np.power(o[:], 2)
                    pos = system.part[ind].pos
                    v2[k,:] = v2[k,:] + np.power(v[:], 2)
                    dr2[k,:] = np.power((pos[:] - pos0[ind,:]), 2)
                    dt = (int_steps * (i + 1) + therm_steps) * \
                        system.time_step
                    # translational diffusion variance: after a closed-form
                    # integration of the Langevin EOM;
                    # ref. the eq. (10.2.26) N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010)
                    # after simple transformations and the dimensional model
                    # matching (cf. eq. (10.1.1) there):
                    sigma2_tr[k] = 0.0
                    for j in range(3):
                        sigma2_tr[k] += self.D_tran_p_validate[k,
                                                               j] * (
                                                                   2.0 * dt + dt0[
                                                                       k,
                                                                       j] * (
                                                                           - 3.0 + 4.0 * math.exp(
                                                                               - dt / dt0[
                                                                                   k,
                                                                                                                            j]) - math.exp(- 2.0 * dt / dt0[k,
                                                                                                                                                            j])))
                    dr_norm[k] += (sum(dr2[k,:]) -
                                   sigma2_tr[k]) / sigma2_tr[k]

        tolerance = 0.15
        Ev = 0.5 * self.mass * v2 / (n * loops)
        Eo = 0.5 * self.J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k,:]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do[k] = sum(Eo[k,:]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do_vec[k,:] = Eo[k,:] / self.halfkT_p_validate[k] - 1.0
        dr_norm /= (n * loops)

        for k in range(2):
            self.assertLessEqual(
                abs(
                    dv[k]),
                tolerance,
                msg='Relative deviation in translational energy too large: {0}'.format(
                    dv[k]))
            if "ROTATION" in espressomd.features():
                self.assertLessEqual(
                    abs(
                        do[k]),
                    tolerance,
                    msg='Relative deviation in rotational energy too large: {0}'.format(
                        do[k]))
                self.assertLessEqual(abs(
                    do_vec[k, 0]), tolerance, msg='Relative deviation in rotational energy per the body axis X is too large: {0}'.format(do_vec[k, 0]))
                self.assertLessEqual(abs(
                    do_vec[k, 1]), tolerance, msg='Relative deviation in rotational energy per the body axis Y is too large: {0}'.format(do_vec[k, 1]))
                self.assertLessEqual(abs(
                    do_vec[k, 2]), tolerance, msg='Relative deviation in rotational energy per the body axis Z is too large: {0}'.format(do_vec[k, 2]))
            self.assertLessEqual(
                abs(
                    dr_norm[k]),
                tolerance,
                msg='Relative deviation in translational diffusion is too large: {0}'.format(
                    dr_norm[k]))

    def set_particle_specific_gamma(self, n):
        """
        Set the particle-specific gamma.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            # for the expected metrics calc
            self.gamma_tran_p_validate[k,:] = self.gamma_tran_p[k,:]
            self.gamma_rot_p_validate[k,:] = self.gamma_rot_p[k,:]
            # init
            for i in range(n):
                ind = i + k * n
                self.system.part[ind].gamma = self.gamma_tran_p[k,:]
                if "ROTATION" in espressomd.features():
                    self.system.part[ind].gamma_rot = self.gamma_rot_p[k,:]

    def set_particle_specific_temperature(self, n):
        """
        Set the particle-specific temperature.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            # expected
            self.halfkT_p_validate[k] = self.kT_p[k] / 2.0
            # init
            for i in range(n):
                ind = i + k * n
                self.system.part[ind].temp = self.kT_p[k]

    def set_diffusivity_tran(self):
        """
        Set the translational diffusivity to validate further.

        """

        for k in range(2):
            # Translational diffusivity for a validation
            self.D_tran_p_validate[k,:] = 2.0 * \
                self.halfkT_p_validate[k] / self.gamma_tran_p_validate[k,:]

    # Test case 0.0.0:
    # no particle specific values / dissipation only / LD only.
    # No meaning for the simple BD propagation cause
    # it has no inertial features.
    def test_case_000(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 0.0.1:
    # no particle specific values / dissipation viscous drag only / BD only.
    # LD will require too much computational time
    # (one is tested offline though).
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_001(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_param_setup()
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag()

    # Test case 0.1: no particle specific values / fluctuation & dissipation
    # Same particle and thermostat parameters for LD and BD are required in order
    # to test the BD consistency by means of NVT-ensemble.
    # Less number of steps and other specific integration parameters of BD
    # reflect its temporal scale advances.
    def test_case_01(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        therm_steps = 20
        loops = 250
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            # Large time-step is OK for BD.
            system.time_step = 42
            # Less number of loops are needed in case of BD because the velocity
            # distribution is already as required. It is not a result of a real dynamics.
            loops = 8
            # The BD does not require so the warmup. Only 1 step is enough.
            # More steps are taken just to be sure that they will not lead
            # to wrong results.
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 1.0.0: particle specific gamma but not temperature / dissipation
    # only / LD only
    def test_case_100(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 1.0.1: particle specific gamma but not temperature /
    # dissipation viscous drag only / BD only.
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_101(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_param_setup()
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            self.set_particle_specific_gamma(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag()

    # Test case 1.1: particle specific gamma but not temperature / fluctuation
    # & dissipation / LD and BD
    def test_case_11(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        therm_steps = 20
        loops = 250
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 42
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 2.0.0: particle specific temperature but not gamma / dissipation
    # only / LD only
    def test_case_200(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 2.0.1: particle specific temperature but not gamma / dissipation
    # viscous drag only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_201(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_param_setup()
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag()

    # Test case 2.1: particle specific temperature but not gamma / fluctuation
    # & dissipation / LD and BD
    def test_case_21(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        therm_steps = 20
        loops = 250
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 42
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.turn_off()
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 3.0.0: both particle specific gamma and temperature /
    # dissipation only / LD only
    def test_case_300(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 3.0.1: both particle specific gamma and temperature /
    # dissipation viscous drag only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_301(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_param_setup()
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            self.set_particle_specific_gamma(n)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag()

    # Test case 3.1: both particle specific gamma and temperature /
    # fluctuation & dissipation / LD and BD
    def test_case_31(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        therm_steps = 20
        loops = 250
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 42
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(kT=self.kT, gamma=self.gamma_global)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 4.0.0: no particle specific values / rotational specific global
    # thermostat / dissipation only / LD only
    def test_case_400(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT,
            gamma=self.gamma_global,
            gamma_rotation=self.gamma_global_rot)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 4.0.1: no particle specific values / rotational specific global
    # thermostat / dissipation only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_401(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_param_setup()
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults_rot_differ()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag()

    # Test case 4.1: no particle specific values / rotational specific global
    # thermostat / fluctuation & dissipation / LD and BD
    def test_case_41(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        therm_steps = 20
        loops = 250
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT,
            gamma=self.gamma_global,
            gamma_rotation=self.gamma_global_rot)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 42
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

if __name__ == '__main__':
    ut.main()
