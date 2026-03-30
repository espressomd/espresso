#
# Copyright (C) 2010-2026 The ESPResSo project
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

"""
Test that the stress tensor transforms correctly under rotations,
even for non-central forces (transverse DPD, angular bonds).

For a rotation matrix R, the stress tensor must satisfy:
    sigma'  =  R . sigma . R^T

where sigma is the stress tensor of the original system and sigma'
is the stress tensor of the rotated system.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np
from scipy.spatial.transform import Rotation

import espressomd
import espressomd.interactions


N_ROTATIONS = 5
BOX_L = 6.0

system = espressomd.System(box_l=[BOX_L] * 3)
system.time_step = 0.01
system.cell_system.skin = 0.4


def random_rotation_matrix(rng):
    """Return a random rotation matrix drawn uniformly from SO(3)."""
    return Rotation.random(random_state=rng).as_matrix()


def rotate_system(R, rotate_velocities=False):
    """Rotate all particle positions (and optionally velocities) by
    rotation matrix R around the center of the box."""
    center = np.array(system.box_l) / 2.0
    partcls = system.part.all()
    partcls.pos = (R @ (partcls.pos - center).T).T + center
    if rotate_velocities:
        partcls.v = (R @ partcls.v.T).T


def check_non_central_forces(test_case, particle_index, ref_direction):
    """Assert that the force on a particle has a non-zero component
    transverse to ref_direction, confirming non-central forces."""
    f = system.part.all().f[particle_index]
    d_hat = ref_direction / np.linalg.norm(ref_direction)
    f_transverse = f - np.dot(f, d_hat) * d_hat
    test_case.assertGreater(
        np.linalg.norm(f_transverse), 0.1,
        "Force is too close to central; test requires non-central forces")


def check_stress_tensor_rotation_covariance(
        test_case, label, rng_seed, rotate_velocities=False):
    """Verify that the stress tensor transforms as R.sigma.R^T under
    several random rotations."""
    system.integrator.run(0)
    stress_orig = np.copy(system.analysis.pressure_tensor()['total'])

    # Verify the stress tensor has non-zero off-diagonal elements,
    # otherwise the rotation covariance check would be trivially satisfied
    off_diag = stress_orig - np.diag(np.diag(stress_orig))
    test_case.assertGreater(
        np.linalg.norm(off_diag), 0.01,
        f"{label}: stress tensor has no significant off-diagonal elements; "
        "test would be trivial")

    rng = np.random.default_rng(seed=rng_seed)
    for i in range(N_ROTATIONS):
        R = random_rotation_matrix(rng)

        rotate_system(R, rotate_velocities=rotate_velocities)
        system.integrator.run(0)
        stress_rotated = np.copy(system.analysis.pressure_tensor()['total'])

        stress_expected = R @ stress_orig @ R.T
        np.testing.assert_allclose(
            stress_rotated, stress_expected, atol=1e-10,
            err_msg=f"{label}: stress tensor rotation covariance "
            f"failed for rotation {i}")

        # Rotate back to original for next iteration
        rotate_system(R.T, rotate_velocities=rotate_velocities)


@utx.skipIfMissingFeatures("DPD")
class StressTensorCovarianceDPD(ut.TestCase):
    """Stress tensor rotation covariance for transverse DPD forces."""

    system = system

    def setUp(self):
        self.system.part.clear()
        self.system.non_bonded_inter.reset()
        self.system.thermostat.turn_off()
        self.system.bonded_inter.clear()

    def tearDown(self):
        self.system.part.clear()
        self.system.non_bonded_inter.reset()
        self.system.thermostat.turn_off()
        self.system.bonded_inter.clear()

    def test_stress_tensor_rotation_covariance(self):
        center = np.array(system.box_l) / 2.0

        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=5.0, r_cut=1.5,
            trans_weight_function=0, trans_gamma=3.0, trans_r_cut=1.5)

        # DPD thermostat with kT=0 so forces are purely dissipative
        # (deterministic), no random kicks
        system.thermostat.set_dpd(kT=0., seed=42)

        # Place two particles close together with non-parallel velocities
        # so that the transverse DPD component is non-zero
        d = np.array([0.7, 0.3, 0.2])
        system.part.add(pos=center - d / 2, v=[1.0, -2.0, 0.5])
        system.part.add(pos=center + d / 2, v=[-0.5, 1.5, -1.0])

        system.integrator.run(0)
        check_non_central_forces(self, particle_index=0, ref_direction=d)
        check_stress_tensor_rotation_covariance(
            self, "DPD", rng_seed=87, rotate_velocities=True)


class StressTensorCovarianceAngleBond(ut.TestCase):
    """Stress tensor rotation covariance for angular bond forces."""

    system = system

    def setUp(self):
        self.system.part.clear()
        self.system.bonded_inter.clear()
        self.system.non_bonded_inter.reset()
        self.system.thermostat.turn_off()

    def tearDown(self):
        self.system.part.clear()
        self.system.bonded_inter.clear()
        self.system.non_bonded_inter.reset()
        self.system.thermostat.turn_off()

    def test_stress_tensor_rotation_covariance(self):
        center = np.array(system.box_l) / 2.0

        # Triangle of three particles with an angular bond.
        # Equilibrium angle differs from actual angle to produce forces.
        d01 = np.array([0.8, 0.3, 0.1])
        d02 = np.array([-0.2, 0.9, -0.4])

        angle_bond = espressomd.interactions.AngleHarmonic(
            bend=100.0, phi0=0.5 * np.pi)
        system.bonded_inter.add(angle_bond)

        p0 = system.part.add(pos=center)
        system.part.add(pos=center + d01)
        system.part.add(pos=center + d02)

        # The central particle holds the angle bond
        p0.add_bond((angle_bond, 1, 2))

        system.integrator.run(0)
        check_non_central_forces(self, particle_index=1, ref_direction=d01)
        check_stress_tensor_rotation_covariance(
            self, "Angle bond", rng_seed=42)


if __name__ == "__main__":
    ut.main()
