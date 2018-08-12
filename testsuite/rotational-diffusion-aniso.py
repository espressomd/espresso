from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import random, seed
import espressomd
import math
import tests_common as tc

@ut.skipIf(not espressomd.has_features(["ROTATION",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "DIPOLES"]),
           "Features not available, skipping test!")
class RotDiffAniso(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System(box_l=[1.0,1.0,1.0])
    es.cell_system.skin = 5.0
    es.seed = 5
    
    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))
    
    # Particle properties
    J = 0.0,0.0,0.0
    
    @classmethod
    def setUpClass(cls):
        np.random.seed(4)

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

    def rot_diffusion_param_setup(self,n):
        """
        Setup the parameters for the rotational diffusion
        test check_rot_diffusion().
    
        Parameters
        ----------
        n : :obj:`int`
            Number of particles.

        """

        ## Time
        # The time step should be less than t0 ~ mass / gamma
        self.es.time_step = 3E-3
        
        ## Space
        box = 10.0
        self.es.box_l = box,box,box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0,0,0
        
        ## NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = self.generate_scalar_ranged_rnd(0.3,5)
        # Note: here & hereinafter specific variations in the random parameter ranges are related to
        # the test execution duration to achieve the required statistical averages faster.
        # The friction gamma_global should be large enough in order to have the small enough D ~ kT / gamma and
        # to observe the details of the original rotational diffusion: the Perrin1936 (see the reference below) tests are visible
        # only when the diffusive rotation is ~pi due to the exponential temporal dependencies (see the equations referred in the check_rot_diffusion()).
        # Also, t0 ~ J / gamma should be small enough
        # in order to remove the small-time-scale diffusion effects which do not fit the Perrin1936's
        # tests which are based on the partial differential equation (eq. (68), Perrin1934) leading only to the simple
        # classical Einstein-Smoluchowski equations of the diffusion
        # in a contrast of the eq. (10.2.26) [N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010)].
        self.gamma_global = 1E2 * self.generate_vec_ranged_rnd(0.5,0.7)
        
        ## Particles
        # As far as the problem characteristic time is t0 ~ J / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the moment of inertia higher than
        # some minimal value.
        # Also, it is expected to test the large enough J.
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time: the in silico time should clock over the t0.
        self.J = self.generate_vec_ranged_rnd(0.1,15.)
        for ind in range(n):
            part_pos = np.random.random(3) * box
            self.es.part.add(rotation=(1,1,1), id=ind, rinertia=self.J,
                                pos=part_pos)
            if "ROTATION" in espressomd.features():
                self.es.part[ind].omega_body = 0.0,0.0,0.0

    def check_rot_diffusion(self,n):
        """
        The rotational diffusion tests based on the reference work
        [Perrin, F. (1936) Journal de Physique et Le Radium, 7(1), 1-11. https://doi.org/10.1051/jphysrad:01936007010100]
        with a theoretical background of
        [Perrin, F. (1934) Journal de Physique et Le Radium, 5(10), 497-511. https://doi.org/10.1051/jphysrad:01934005010049700]
        
        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """
        # Global diffusivity tensor in the body frame:
        D = self.kT / self.gamma_global
        dt0 = self.J / self.gamma_global

        # Thermalizing...
        therm_steps = int(100)
        self.es.integrator.run(therm_steps)

        # Measuring...
        # Set the initial conditions according to the [Perrin1936], p.3.
        # The body angular velocity is rotated now, but there is
        # only the thermal velocity, hence, this does not impact the test and its physical context.
        for ind in range(n):
            self.es.part[ind].quat = 1.0, 0.0, 0.0, 0.0
        ## Average direction cosines
        # Diagonal ones:
        dcosjj_validate = np.zeros((3))
        dcosjj_dev = np.zeros((3))
        # same to the power of 2
        dcosjj2_validate = np.zeros((3))
        dcosjj2_dev = np.zeros((3))
        # The non-diagonal elements for 2 different tests: negative ("nn") and positive ("pp") ones.
        dcosijpp_validate = np.ones((3,3))
        dcosijpp_dev = np.zeros((3,3))
        dcosijnn_validate = np.ones((3,3))
        dcosijnn_dev = np.zeros((3,3))
        # The non-diagonal elements to the power of 2
        dcosij2_validate = np.ones((3,3))
        dcosij2_dev = np.zeros((3,3))
        
        self.es.time = 0.0
        int_steps = 10
        loops = 100
        for step in range(loops):
            self.es.integrator.run(steps = int_steps)
            dcosjj = np.zeros((3))
            dcosjj2 = np.zeros((3))
            dcosijpp = np.zeros((3,3))
            dcosijnn = np.zeros((3,3))
            dcosij2 = np.zeros((3,3))
            for ind in range(n):
                # Just a direction cosines functions averaging..
                dir_cos = tc.rotation_matrix_quat(self.es, ind)
                for j in range(3):
                    # the LHS of eq. (23) [Perrin1936].
                    dcosjj[j] += dir_cos[j,j]
                    # the LHS of eq. (32) [Perrin1936].
                    dcosjj2[j] += dir_cos[j,j]**2.0
                    for i in range(3):
                        if i != j:
                            # the LHS of eq. (24) [Perrin1936].
                            dcosijpp[i,j] += dir_cos[i,i]*dir_cos[j,j]+dir_cos[i,j]*dir_cos[j,i]
                            # the LHS of eq. (25) [Perrin1936].
                            dcosijnn[i,j] += dir_cos[i,i]*dir_cos[j,j]-dir_cos[i,j]*dir_cos[j,i]
                            # the LHS of eq. (33) [Perrin1936].
                            dcosij2[i,j] += dir_cos[i,j]**2.0
            dcosjj /= n
            dcosjj2 /= n
            dcosijpp /= n
            dcosijnn /= n
            dcosij2 /= n
            
            ### Actual comparison.
            
            tolerance = 0.195
            # Too small values of the direction cosines are out of interest compare to 0..1 range.
            min_value = 0.14
            
            # Eq. (23) [Perrin1936].
            dcosjj_validate[0] = np.exp(-(D[1] + D[2]) * self.es.time)
            dcosjj_validate[1] = np.exp(-(D[0] + D[2]) * self.es.time)
            dcosjj_validate[2] = np.exp(-(D[0] + D[1]) * self.es.time)
            dcosjj_dev = np.absolute(dcosjj - dcosjj_validate) / dcosjj_validate
            for j in range(3):
                if np.absolute(dcosjj_validate[j]) < min_value:
                    dcosjj_dev[j] = 0.0
            
            # Eq. (24) [Perrin1936].
            dcosijpp_validate[0,1] = np.exp(-(4 * D[2] + D[1] + D[0]) * self.es.time)
            dcosijpp_validate[1,0] = np.exp(-(4 * D[2] + D[1] + D[0]) * self.es.time)
            dcosijpp_validate[0,2] = np.exp(-(4 * D[1] + D[2] + D[0]) * self.es.time)
            dcosijpp_validate[2,0] = np.exp(-(4 * D[1] + D[2] + D[0]) * self.es.time)
            dcosijpp_validate[1,2] = np.exp(-(4 * D[0] + D[2] + D[1]) * self.es.time)
            dcosijpp_validate[2,1] = np.exp(-(4 * D[0] + D[2] + D[1]) * self.es.time)
            dcosijpp_dev = np.absolute(dcosijpp - dcosijpp_validate) / dcosijpp_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijpp_validate[i,j]) < min_value:
                        dcosijpp_dev[i,j] = 0.0
            
            # Eq. (25) [Perrin1936].
            dcosijnn_validate[0,1] = np.exp(-(D[1] + D[0]) * self.es.time)
            dcosijnn_validate[1,0] = np.exp(-(D[1] + D[0]) * self.es.time)
            dcosijnn_validate[0,2] = np.exp(-(D[2] + D[0]) * self.es.time)
            dcosijnn_validate[2,0] = np.exp(-(D[2] + D[0]) * self.es.time)
            dcosijnn_validate[1,2] = np.exp(-(D[2] + D[1]) * self.es.time)
            dcosijnn_validate[2,1] = np.exp(-(D[2] + D[1]) * self.es.time)
            dcosijnn_dev = np.absolute(dcosijnn - dcosijnn_validate) / dcosijnn_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijnn_validate[i,j]) < min_value:
                        dcosijnn_dev[i,j] = 0.0
            
            # Eq. (30) [Perrin1936].
            D0 = sum(D[:]) / 3.0
            D1D1 = 0.0
            for j in range(3):
                for i in range(3):
                    if i != j:
                        D1D1 += D[i]*D[j]
            D1D1 /= 6.0
            # Eq. (32) [Perrin1936].
            dcosjj2_validate = 1./3. + (1./3.) * (1. + (D-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0-np.sqrt(D0**2-D1D1))*self.es.time) \
                                     + (1./3.) * (1. - (D-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0+np.sqrt(D0**2-D1D1))*self.es.time)
            dcosjj2_dev = np.absolute(dcosjj2 - dcosjj2_validate) / dcosjj2_validate
            for j in range(3):
                if np.absolute(dcosjj2_validate[j]) < min_value:
                    dcosjj2_dev[j] = 0.0
            
            # Eq. (33) [Perrin1936].
            dcosij2_validate[0,1] = 1./3. - (1./6.) * (1. - (D[2]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0-np.sqrt(D0**2-D1D1))*self.es.time) \
                                          - (1./6.) * (1. + (D[2]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0+np.sqrt(D0**2-D1D1))*self.es.time)
            dcosij2_validate[1,0] = dcosij2_validate[0,1]
            
            dcosij2_validate[0,2] = 1./3. - (1./6.) * (1. - (D[1]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0-np.sqrt(D0**2-D1D1))*self.es.time) \
                                          - (1./6.) * (1. + (D[1]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0+np.sqrt(D0**2-D1D1))*self.es.time)
            dcosij2_validate[2,0] = dcosij2_validate[0,2]
            
            dcosij2_validate[1,2] = 1./3. - (1./6.) * (1. - (D[0]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0-np.sqrt(D0**2-D1D1))*self.es.time) \
                                          - (1./6.) * (1. + (D[0]-D0)/(2.*np.sqrt(D0**2-D1D1))) \
                               * np.exp(-6.*(D0+np.sqrt(D0**2-D1D1))*self.es.time)
            dcosij2_validate[2,1] = dcosij2_validate[1,2]
            dcosij2_dev = np.absolute(dcosij2 - dcosij2_validate) / dcosij2_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosij2_validate[i,j]) < min_value:
                        dcosij2_dev[i,j] = 0.0
            
            for j in range(3):
                self.assertLessEqual(
                        abs(
                            dcosjj_dev[j]),
                        tolerance,
                        msg='Relative deviation dcosjj_dev[{0}] in a rotational diffusion is too large: {1}'.format(
                            j,dcosjj_dev[j]))
                self.assertLessEqual(
                        abs(
                            dcosjj2_dev[j]),
                        tolerance,
                        msg='Relative deviation dcosjj2_dev[{0}] in a rotational diffusion is too large: {1}'.format(
                            j,dcosjj2_dev[j]))
                for i in range(3):
                    if i != j:
                        self.assertLessEqual(
                            abs(
                                dcosijpp_dev[i,j]),
                            tolerance,
                            msg='Relative deviation dcosijpp_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i,j,dcosijpp_dev[i,j]))
                        self.assertLessEqual(
                            abs(
                                dcosijnn_dev[i,j]),
                            tolerance,
                            msg='Relative deviation dcosijnn_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i,j,dcosijnn_dev[i,j]))
                        self.assertLessEqual(
                            abs(
                                dcosij2_dev[i,j]),
                            tolerance,
                            msg='Relative deviation dcosij2_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i,j,dcosij2_dev[i,j]))

    def test_case_00(self):
        n = 400
        self.rot_diffusion_param_setup(n)
        self.es.thermostat.set_langevin(kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.check_rot_diffusion(n)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
