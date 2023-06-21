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

import unittest as ut
import unittest_decorators as utx
import thermostats_common
import numpy as np

import espressomd.lb

"""
Check the lattice-Boltzmann thermostat with respect to the particle velocity
distribution.
"""

KT = 0.2
AGRID = 0.8
node_volume = AGRID**3
KIN_VISC = 0.3
DENS = 0.8
TIME_STEP = 0.005

LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'kinematic_viscosity': KIN_VISC,
             'tau': TIME_STEP,
             'kT': KT,
             'seed': 123}


class LBThermostatCommon(thermostats_common.ThermostatsCommon):

    """
    Check the LB thermostat.
    All temperature checks rely on https://en.wikipedia.org/wiki/Thermal_velocity
    with the root-mean-squared-velocity.
    The temperature check together with the friction test
    ensure that the noise amplitude is correct.
    """

    system = espressomd.System(box_l=[AGRID * 12] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    # use "small" particles such that radius=hydrodynamic radius
    global_gamma = 0.1 * AGRID * 6 * np.pi * KIN_VISC * DENS
    partcl_gamma = 0.05 * AGRID * 6 * np.pi * KIN_VISC * DENS
    # relaxation time 0.2 sim time units
    partcl_mass = 0.2 * partcl_gamma

    np.random.seed(41)

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf, seed=5, gamma=self.global_gamma)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def get_lb_kT(self, lbf):
        nodes_mass = lbf[:, :, :].density * node_volume
        nodes_vel_sq = np.sum(np.square(lbf[:, :, :].velocity), axis=3)
        return np.mean(nodes_mass * nodes_vel_sq) / 3.

    def get_lb_velocity(self, lbf):
        return np.mean(lbf[:, :, :].velocity, axis=(0, 1, 2))

    def test_fluid(self):
        self.system.integrator.run(100)
        fluid_kTs = []
        for _ in range(100):
            fluid_kTs.append(self.get_lb_kT(self.lbf))
            self.system.integrator.run(3)

        fluid_kT = np.mean(fluid_kTs)
        np.testing.assert_allclose(fluid_kT, KT, rtol=0.05)

    def check_partcl_temp(self, partcl_vel):
        partcl_vel_rms = np.sqrt(
            np.mean(np.linalg.norm(partcl_vel, axis=2)**2))

        np.testing.assert_allclose(
            np.mean(partcl_vel), 0, atol=0.05 * partcl_vel_rms)
        np.testing.assert_allclose(partcl_vel_rms**2 * self.partcl_mass / 3.,
                                   KT, rtol=0.05)

        vel_range = 2 * partcl_vel_rms
        n_bins = 7
        self.check_velocity_distribution(
            partcl_vel.reshape((-1, 3)), vel_range, n_bins, 0.02, KT,
            mass=self.partcl_mass)

    @utx.skipIfMissingFeatures(["MASS"])
    def test_temperature_with_particles(self):
        n_partcls_per_type = 50
        partcls_global_gamma = self.system.part.add(
            pos=np.random.random((n_partcls_per_type, 3)) * self.system.box_l,
            mass=n_partcls_per_type * [self.partcl_mass])
        partcls_per_part_gamma = self.system.part.add(
            pos=np.random.random((n_partcls_per_type, 3)) * self.system.box_l,
            mass=n_partcls_per_type * [self.partcl_mass],
            gamma=n_partcls_per_type * [self.partcl_gamma])
        self.system.integrator.run(50)
        n_samples = 400
        vel_global_gamma = np.zeros((n_samples, n_partcls_per_type, 3))
        vel_per_part_gamma = np.zeros_like(vel_global_gamma)
        fluid_kTs = []

        for i in range(n_samples):
            self.system.integrator.run(3)
            if i % 10 == 0:
                fluid_kTs.append(self.get_lb_kT(self.lbf))
            vel_global_gamma[i] = partcls_global_gamma.v
            vel_per_part_gamma[i] = partcls_per_part_gamma.v
        fluid_kT = np.mean(fluid_kTs)

        self.check_partcl_temp(vel_global_gamma)
        self.check_partcl_temp(vel_per_part_gamma)

        np.testing.assert_allclose(fluid_kT, KT, rtol=0.03)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_friction(self):
        """apply force and measure if the average velocity matches expectation"""

        # large force to get better signal  to noise ratio
        ext_force = np.array([1.2, 2.1, 1.1])
        n_partcls_each_type = 50
        partcls_global_gamma = self.system.part.add(
            pos=np.random.random((n_partcls_each_type, 3)) * self.system.box_l,
            ext_force=n_partcls_each_type * [ext_force],
            mass=n_partcls_each_type * [self.partcl_mass])
        partcls_per_partcl_gamma = self.system.part.add(
            pos=np.random.random((n_partcls_each_type, 3)) * self.system.box_l,
            ext_force=n_partcls_each_type * [ext_force],
            mass=n_partcls_each_type * [self.partcl_mass],
            gamma=n_partcls_each_type * [self.partcl_gamma])

        # add counterforce to fluid such that velocities cannot increase
        # indefinitely and we get a steady state relative velocity
        total_force_partcls = ext_force * 2 * n_partcls_each_type
        lb_volume = np.prod(self.system.box_l)
        self.lbf.ext_force_density = -total_force_partcls / lb_volume

        # let particles and fluid accelerate towards steady state before
        # measuring
        self.system.integrator.run(int(0.8 / self.system.time_step))
        startpos_global_gamma = np.copy(partcls_global_gamma.pos)
        startpos_per_p_gamma = np.copy(partcls_per_partcl_gamma.pos)

        run_time = 0.3
        n_steps = int(run_time / self.system.time_step)
        steps_per_sample = 10
        lb_vels = []
        for _ in range(int(n_steps / steps_per_sample)):
            self.system.integrator.run(steps_per_sample)
            lb_vels.append(self.get_lb_velocity(self.lbf))
        lb_vel = np.mean(lb_vels, axis=0)

        vel_global_gamma = np.mean(
            partcls_global_gamma.pos - startpos_global_gamma,
            axis=0) / run_time
        vel_per_part_gamma = np.mean(
            partcls_per_partcl_gamma.pos - startpos_per_p_gamma,
            axis=0) / run_time
        # get velocity relative to the fluid
        vel_global_gamma -= lb_vel
        vel_per_part_gamma -= lb_vel

        # average over 3 spatial directions (isotropic friction)
        global_gamma_measured = np.mean(ext_force / vel_global_gamma)
        per_partcl_gamma_measured = np.mean(ext_force / vel_per_part_gamma)

        np.testing.assert_allclose(
            global_gamma_measured,
            self.global_gamma,
            rtol=0.05)
        np.testing.assert_allclose(
            per_partcl_gamma_measured,
            self.partcl_gamma,
            rtol=0.05)

    @utx.skipIfMissingFeatures(["PARTICLE_ANISOTROPY",
                               "THERMOSTAT_PER_PARTICLE"])
    def test_exceptions(self):
        self.system.part.add(pos=[0., 0., 0.], gamma=[1., 2., 3.], id=2)
        with self.assertRaisesRegex(Exception, r"ERROR: anisotropic particle \(id 2\) coupled to LB"):
            self.system.integrator.run(1)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaThermostat(LBThermostatCommon, ut.TestCase):

    """Test for the CPU implementation of the LB."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaThermostatSinglePrecision(LBThermostatCommon, ut.TestCase):

    """Test for the CPU implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}


# @utx.skipIfMissingGPU()
# @utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
# class LBWalberlaGPUThermostat(LBThermostatCommon, ut.TestCase):

#    """Test for the GPU implementation of the LB."""

#    lb_class = espressomd.lb.LBFluidWalberlaGPU
#    lb_params = {"single_precision": True}


if __name__ == '__main__':
    ut.main()
