import unittest as ut
import numpy as np
import copy

import espressomd
import espressomd.lb
import espressomd.observables

"""
Tests for the LB fluid profile observables.

"""

BOX_L_X = 12.0
BOX_L_Y = 12.0
BOX_L_Z = 12.0
TIME_STEP = 0.1
AGRID = 0.5
VISC = .7
DENS = 1.7
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'fric': 1.0,
             'tau': TIME_STEP
             }

LB_VELOCITY_PROFILE_PARAMS = {
    'n_x_bins': int(BOX_L_X / AGRID),
    'n_y_bins': int(BOX_L_Y / AGRID),
    'n_z_bins': int(BOX_L_Z / AGRID),
    'min_x': 0.0,
    'min_y': 0.0,
    'min_z': 0.0,
    'max_x': BOX_L_X,
    'max_y': BOX_L_Y,
    'max_z': BOX_L_Z,
    'sampling_delta_x': AGRID,
    'sampling_delta_y': AGRID,
    'sampling_delta_z': AGRID,
    'sampling_offset_x': 0.5 * AGRID,
    'sampling_offset_y': 0.5 * AGRID,
    'sampling_offset_z': 0.5 * AGRID}

class ObservableProfileLBCommon(object):
    lbf = None
    system = espressomd.System(box_l=[12.0, 12.0, 12.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def set_fluid_velocities(self):
        """Set an x dependent fluid velocity."""
        for x in range(int(self.system.box_l[0] / AGRID)):
            for y in range(int(self.system.box_l[1] / AGRID)):
                for z in range(int(self.system.box_l[2] / AGRID)):
                    self.lbf[x, y, z].velocity = [float(x), 0.0, 0.0]

    def test_velocity_profile(self):
        self.set_fluid_velocities()
        obs = espressomd.observables.LBVelocityProfile(**LB_VELOCITY_PROFILE_PARAMS)
        obs_data = np.array(obs.calculate())
        obs_data = obs_data.reshape((LB_VELOCITY_PROFILE_PARAMS['n_x_bins'],
                                     LB_VELOCITY_PROFILE_PARAMS['n_y_bins'],
                                     LB_VELOCITY_PROFILE_PARAMS['n_z_bins'], 3))
        for x in range(obs_data.shape[0]):
            for y in range(obs_data.shape[1]):
                for z in range(obs_data.shape[2]):
                    self.assertAlmostEqual(obs_data[x, y, z, 0], float(x), places=5)
        self.assertEqual(obs.n_values(), LB_VELOCITY_PROFILE_PARAMS['n_x_bins'] * LB_VELOCITY_PROFILE_PARAMS['n_y_bins'] * LB_VELOCITY_PROFILE_PARAMS['n_z_bins'] * 3)

    def test_error_sampling_delta_of_0(self):
        lb_velocity_params_local = copy.copy(LB_VELOCITY_PROFILE_PARAMS)
        lb_velocity_params_local['sampling_delta_x'] = 0.0
        lb_velocity_params_local['sampling_delta_y'] = 0.0
        lb_velocity_params_local['sampling_delta_z'] = 0.0
        with self.assertRaises(RuntimeError):
            obs2 = espressomd.observables.LBVelocityProfile(**lb_velocity_params_local)

    def test_error_if_no_LB(self):
        self.system.actors.clear()
        obs = espressomd.observables.LBVelocityProfile(**LB_VELOCITY_PROFILE_PARAMS)
        with self.assertRaises(RuntimeError):
            obs.calculate()

    def test_error_if_empty_bin(self):
        lb_velocity_params_local = copy.copy(LB_VELOCITY_PROFILE_PARAMS)
        lb_velocity_params_local['sampling_delta_x'] = 3.0
        obs = espressomd.observables.LBVelocityProfile(**lb_velocity_params_local)
        with self.assertRaises(RuntimeError):
            obs.calculate()
        

@ut.skipIf(not espressomd.has_features(
    'LB') or espressomd.has_features('SHANCHEN'), "Skipping test due to missing features.")
class LBCPU(ut.TestCase, ObservableProfileLBCommon):
    """Test for the CPU implementation of the LB."""
    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)
        self.system.actors.clear()
        self.system.actors.add(self.lbf)

@ut.skipIf(not espressomd.has_features(
    'LB_GPU') or espressomd.has_features('SHANCHEN'), "Skipping test due to missing features.")
class LBGPU(ut.TestCase, ObservableProfileLBCommon):
    """Test for the GPU implementation of the LB."""
    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)
        self.system.actors.clear()
        self.system.actors.add(self.lbf)


if __name__ == "__main__":
    ut.main()
