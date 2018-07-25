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
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])
    es.seed = es.cell_system.get_state()['n_nodes'] * [1234]
    es.cell_system.skin = 5.0
    
    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))
    gamma_global_rot = np.zeros((3))
    
    # Particle properties
    mass = 0.0
    J = 0.0,0.0,0.0
    
    ## Per-particle type parameters.
    # 2 different langevin parameters for particles.
    kT_p = np.zeros((2))
    # gamma_tran/gamma_rot matrix: [2 kinds of particles] x [3 dimensions X Y Z]
    # These matrices are assigning per-particle in corresponding test cases.
    gamma_tran_p = np.zeros((2, 3))
    gamma_rot_p = np.zeros((2, 3))
    
    ## These variables will take the values to compare with.
    # Depending on the test case following matrices either equals to the previous
    # or the global corresponding parameters. The corresponding setting effect is an essence of
    # all the test cases' differentiation here.
    halfkT_p_validate = np.zeros((2))
    gamma_tran_p_validate = np.zeros((2, 3))
    gamma_rot_p_validate = np.zeros((2, 3))
    # Diffusivity
    D_tran_p_validate = np.zeros((2,3))
    
    @classmethod
    def setUpClass(cls):
        np.random.seed(15)

    def setUp(self):
        self.es.time = 0.0
        self.es.part.clear()
        print("\n")

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
        Setup the NVT thermostat viscous friction parameters.

        """

        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k, :] = self.gamma_global[:]
            self.gamma_rot_p_validate[k, :] = self.gamma_global[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def set_langevin_global_defaults_rot_differ(self):
        """
        Setup the NVT thermostat viscous friction parameters
        with a rotation-specific gamma.

        """

        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k, :] = self.gamma_global[:]
            self.gamma_rot_p_validate[k, :] = self.gamma_global_rot[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def dissipation_param_setup(self):
        """
        Setup the parameters for the following dissipation
        test.

        """

        ## Time
        self.es.time_step = 0.007
        
        ## Space
        box = 1.0
        self.es.box_l = box,box,box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0,0,0
        
        ## NVT thermostat
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
        max_gamma_param = 2.0/3.0
        gamma_rnd = self.generate_scalar_ranged_rnd(min_gamma_param, max_gamma_param)
        self.gamma_global = gamma_rnd, gamma_rnd, gamma_rnd
        # Additional test case for the specific global rotational gamma set.
        self.gamma_global_rot = self.generate_vec_ranged_rnd(0.5,2.0/3.0)
        # Per-paricle values:
        self.kT_p = 0.0,0.0
        # Either translational friction isotropy is required
        # or both translational and rotational ones.
        # Otherwise these types of motion will interfere.
        # ..Let's test both cases depending on the particle index.
        self.gamma_tran_p[0, 0] = self.generate_scalar_ranged_rnd(0.5,1.0)
        self.gamma_tran_p[0, 1] = self.gamma_tran_p[0, 0]
        self.gamma_tran_p[0, 2] = self.gamma_tran_p[0, 0]
        self.gamma_rot_p[0, :] = self.generate_vec_ranged_rnd(0.5,2.0/3.0)
        self.gamma_tran_p[1, 0] = self.generate_scalar_ranged_rnd(0.5,1.0)
        self.gamma_tran_p[1, 1] = self.gamma_tran_p[1, 0]
        self.gamma_tran_p[1, 2] = self.gamma_tran_p[1, 0]
        self.gamma_rot_p[1, 0] = self.generate_scalar_ranged_rnd(0.5,2.0/3.0)
        self.gamma_rot_p[1, 1] = self.gamma_rot_p[1, 0]
        self.gamma_rot_p[1, 2] = self.gamma_rot_p[1, 0]
        
        ## Particles
        self.mass = 12.74
        self.J = 10.0,10.0,10.0
        for i in range(2):
            self.es.part.add(rotation=(1,1,1), pos=(0.0,0.0,0.0), id=i)
            self.es.part[i].v = 1.0,1.0,1.0
            if "ROTATION" in espressomd.features():
                self.es.part[i].omega_body = 1.0,1.0,1.0
            self.es.part[i].mass = self.mass
            self.es.part[i].rinertia = self.J

    def fluctuation_dissipation_param_setup(self,n):
        """
        Setup the parameters for the following fluctuation-dissipation
        test.
    
        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        ## Time
        self.es.time_step = 0.03
        
        ## Space
        box = 10.0
        self.es.box_l = box,box,box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0,0,0
        
        ## NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = self.generate_scalar_ranged_rnd(0.3,5)
        # See the above comment regarding the gamma assignments.
        # Note: here & hereinafter specific variations in these ranges are related to
        # the test execution duration to achieve the required statistical averages faster.
        self.gamma_global = self.generate_vec_ranged_rnd(0.5,2.0/3.0)
        self.gamma_global_rot = self.generate_vec_ranged_rnd(0.2,20)
        # Per-particle parameters
        self.kT_p = 2.5,2.0
        for k in range(2):
            self.gamma_tran_p[k, :] = self.generate_vec_ranged_rnd(0.4,10.0)
            self.gamma_rot_p[k, :] = self.generate_vec_ranged_rnd(0.2,20.0)
        
        ## Particles
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the mass higher than
        # some minimal value according to the value min_mass_param.
        # Also, it is expected to test the large enough mass (max_mass_param).
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time.
        min_mass_param = 0.2
        max_mass_param = 7.0
        self.mass = self.generate_scalar_ranged_rnd(min_mass_param,max_mass_param)
        self.J = self.generate_vec_ranged_rnd(min_mass_param,max_mass_param)
        for i in range(n):
            for k in range(2):
                ind = i + k * n
                part_pos = np.random.random(3) * box
                part_v = 0.0,0.0,0.0
                part_omega_body = 0.0,0.0,0.0
                self.es.part.add(rotation=(1,1,1), id=ind, mass=self.mass, rinertia=self.J,
                                 pos=part_pos, v=part_v)
                if "ROTATION" in espressomd.features():
                    self.es.part[ind].omega_body = part_omega_body

    def check_dissipation(self):
        """
        Check the dissipation relations: the simple viscous decelleration test.

        """

        tol = 1.25E-4
        for i in range(100):
            for k in range(2):
                for j in range(3):
                    # Note: velocity is defined in the lab frame of reference
                    # while gamma_tr is defined in the body one.
                    # Hence, only isotropic gamma_tran_p_validate could be tested here.
                    self.assertLess(
                        abs(self.es.part[k].v[j] - math.exp(- self.gamma_tran_p_validate[k, j] * self.es.time / self.mass)), tol)
                    if "ROTATION" in espressomd.features():
                        self.assertLess(abs(
                            self.es.part[k].omega_body[j] - math.exp(- self.gamma_rot_p_validate[k, j] * self.es.time / self.J[j])), tol)

    def check_fluctuation_dissipation(self,n):
        """
        Check the fluctuation-dissipation relations: thermalization
        and diffusion properties.
    
        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        ## The thermalization and diffusion test
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
                pos0[ind, :] = self.es.part[ind].pos
        dt0 = self.mass / self.gamma_tran_p_validate

        loops = 250
        print("Thermalizing...")
        therm_steps = 20
        self.es.integrator.run(therm_steps)
        print("Measuring...")

        int_steps = 5
        for i in range(loops):
            self.es.integrator.run(int_steps)
            # Get kinetic energy in each degree of freedom for all particles
            for p in range(n):
                for k in range(2):
                    ind = p + k * n
                    v = self.es.part[ind].v
                    if "ROTATION" in espressomd.features():
                        o = self.es.part[ind].omega_body
                        o2[k, :] = o2[k, :] + np.power(o[:], 2)
                    pos = self.es.part[ind].pos
                    v2[k, :] = v2[k, :] + np.power(v[:], 2)
                    dr2[k, :] = np.power((pos[:] - pos0[ind, :]), 2)
                    dt = (int_steps * (i + 1) + therm_steps) * \
                        self.es.time_step
                    # translational diffusion variance: after a closed-form
                    # integration of the Langevin EOM;
                    # ref. the eq. (10.2.26) N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010)
                    # after simple transformations and the dimensional model matching (cf. eq. (10.1.1) there):
                    sigma2_tr[k] = 0.0
                    for j in range(3):
                        sigma2_tr[k] += self.D_tran_p_validate[k,
                                                           j] * (2.0 * dt + dt0[k,
                                                                              j] * (- 3.0 + 4.0 * math.exp(- dt / dt0[k,
                                                                                                                  j]) - math.exp(- 2.0 * dt / dt0[k,
                                                                                                                                                j])))
                    dr_norm[k] += (sum(dr2[k, :]) - sigma2_tr[k]) / sigma2_tr[k]

        tolerance = 0.15
        Ev = 0.5 * self.mass * v2 / (n * loops)
        Eo = 0.5 * self.J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k, :]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do[k] = sum(Eo[k, :]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do_vec[k, :] = Eo[k, :] / self.halfkT_p_validate[k] - 1.0
        dr_norm /= (n * loops)
        
        for k in range(2):
            print("\n")
            print("k = " + str(k))

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

    def set_particle_specific_gamma(self,n):
        """
        Set the particle-specific gamma.
    
        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            self.gamma_tran_p_validate[k, :] = self.gamma_tran_p[k, :]
            self.gamma_rot_p_validate[k, :] = self.gamma_rot_p[k, :]
            for i in range(n):
                ind = i + k * n
                self.es.part[ind].gamma = self.gamma_tran_p[k, :]
                if "ROTATION" in espressomd.features():
                    self.es.part[ind].gamma_rot = self.gamma_rot_p[k, :]

    def set_particle_specific_temperature(self,n):
        """
        Set the particle-specific temperature.
    
        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            self.halfkT_p_validate[k] = self.kT_p[k] / 2.0
            for i in range(n):
                ind = i + k * n
                self.es.part[ind].temp = self.kT_p[k]

    def set_diffusivity_tran(self):
        """
        Set the translational diffusivity to validate further.

        """

        for k in range(2):
            # Translational diffusivity for a validation
            self.D_tran_p_validate[k, :] = 2.0 * self.halfkT_p_validate[k] / self.gamma_tran_p_validate[k, :]

    # Test case 0.0: no particle specific values / dissipation only
    def test_case_00(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.check_dissipation()
    
    # Test case 0.1: no particle specific values / fluctuation & dissipation
    def test_case_01(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n)

    # Test case 1.0: particle specific gamma but not temperature / dissipation only
    def test_case_10(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 1.1: particle specific gamma but not temperature / fluctuation & dissipation
    def test_case_11(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n)

    # Test case 2.0: particle specific temperature but not gamma / dissipation only
    def test_case_20(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 2.1: particle specific temperature but not gamma / fluctuation & dissipation
    def test_case_21(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n)

    # Test case 3.0: both particle specific gamma and temperature / dissipation only
    def test_case_30(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation()

    # Test case 3.1: both particle specific gamma and temperature / fluctuation & dissipation
    def test_case_31(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n)

    # Test case 4.0: no particle specific values / rotational specific global thermostat / dissipation only
    def test_case_40(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup()
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.turn_off()
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global, gamma_rotation=self.gamma_global_rot)
        # Actual integration and validation run
        self.check_dissipation()
    
    # Test case 4.1: no particle specific values / rotational specific global thermostat / fluctuation & dissipation
    def test_case_41(self):
        # Each of 2 kind of particles will be represented by n instances:
        n = 200
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        self.es.thermostat.turn_off()
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global, gamma_rotation=self.gamma_global_rot)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
