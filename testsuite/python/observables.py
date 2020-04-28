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
import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.observables


def calc_com_x(system, x, id_list):
    """Mass-weighted average, skipping virtual sites"""
    masses = system.part[id_list].mass

    # Filter out virtual particles by using mass=0 for them
    virtual = system.part[id_list].virtual
    for i in range(len(masses)):
        if virtual[i]:
            masses[i] = 0.

    com_x = np.average(
        getattr(system.part[id_list], x), weights=masses, axis=0)
    return com_x


class Observables(ut.TestCase):
    N_PART = 200
    # Handle for espresso system
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.part.add(
        id=np.arange(3, 3 + 2 * N_PART, 2),
        pos=np.random.random((N_PART, 3)) * system.box_l,
        v=np.random.random((N_PART, 3)) * 3.2 - 1,
        f=np.random.random((N_PART, 3)))

    if espressomd.has_features(["MASS"]):
        system.part[:].mass = np.random.random(N_PART)

    if espressomd.has_features(["DIPOLES"]):
        system.part[:].dip = np.random.random((N_PART, 3)) - .3

    if espressomd.has_features(["ROTATION"]):
        system.part[:].omega_body = np.random.random((N_PART, 3)) - .5
        system.part[:].torque_lab = np.random.random((N_PART, 3)) - .5
        system.part[:].quat = np.random.random((N_PART, 4))

    if espressomd.has_features("DIPOLES"):
        system.part[:].dipm = np.random.random(N_PART) + 2

    if espressomd.has_features("ELECTROSTATICS"):
        system.part[:].q = np.random.random(N_PART)

    if espressomd.has_features("VIRTUAL_SITES"):
        p = system.part[system.part[:].id[8]]
        p.virtual = True

    def generate_test_for_pid_observable(
            _obs_name, _pprop_name, _agg_type=None):
        """Generates test cases for observables working on particle id lists.

        """
        pprop_name = _pprop_name
        obs_name = _obs_name
        agg_type = _agg_type

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            # Randomly pick a subset of the particles
            id_list = sorted(
                np.random.choice(
                    self.system.part[:].id,
                    size=int(
                        self.N_PART * .9),
                    replace=False))
            for id in id_list:
                assert(self.system.part.exists(id))

            # Get data from particles
            part_data = getattr(self.system.part[id_list], pprop_name)

            # Reshape and aggregate to linear array
            if len(part_data.shape) > 1:
                if agg_type == "average":
                    part_data = np.average(part_data, 0)
                if agg_type == "sum":
                    part_data = np.sum(part_data, 0)
                if agg_type == 'com':
                    part_data = calc_com_x(self.system, pprop_name, id_list)

            # Data from observable
            observable = obs_name(ids=id_list)
            obs_data = observable.calculate()

            # Check
            self.assertEqual(obs_data.shape, part_data.shape)
            np.testing.assert_equal(id_list, observable.ids)

            np.testing.assert_array_almost_equal(
                obs_data,
                part_data, err_msg="Data did not agree for observable " +
                str(obs_name) +
                " and particle property " +
                pprop_name, decimal=11)

        return func

    test_pos = generate_test_for_pid_observable(
        espressomd.observables.ParticlePositions, "pos")
    test_v = generate_test_for_pid_observable(
        espressomd.observables.ParticleVelocities, "v")
    test_f = generate_test_for_pid_observable(
        espressomd.observables.ParticleForces, "f")
    test_com_position = generate_test_for_pid_observable(
        espressomd.observables.ComPosition, 'pos', 'com')
    test_com_velocity = generate_test_for_pid_observable(
        espressomd.observables.ComVelocity, 'v', 'com')

    if espressomd.has_features(["DIPOLES"]):
        test_mag_dip = generate_test_for_pid_observable(
            espressomd.observables.MagneticDipoleMoment, "dip", "sum")

    if espressomd.has_features(["ROTATION"]):
        test_body_angular_velocity = generate_test_for_pid_observable(
            espressomd.observables.ParticleBodyAngularVelocities, "omega_body")
        test_lab_angular_velocity = generate_test_for_pid_observable(
            espressomd.observables.ParticleAngularVelocities, "omega_lab")

    @utx.skipIfMissingFeatures(['ROTATION'])
    def test_particle_body_velocities(self):
        obs = espressomd.observables.ParticleBodyVelocities(
            ids=self.system.part[:].id)
        obs_data = obs.calculate()
        part_data = np.array([p.convert_vector_space_to_body(p.v)
                              for p in self.system.part])
        self.assertEqual(obs_data.shape, part_data.shape)
        np.testing.assert_array_almost_equal(part_data, obs_data,
                                             err_msg="Data did not agree for observable ParticleBodyVelocities and particle derived values.",
                                             decimal=9)

    def test_stress_tensor(self):
        s = self.system.analysis.stress_tensor()["total"]
        obs_data = espressomd.observables.StressTensor().calculate()
        self.assertEqual(obs_data.shape, s.shape)
        np.testing.assert_array_almost_equal(
            s,
            obs_data,
            err_msg="Stress tensor from analysis and observable did not agree",
            decimal=9)

    @utx.skipIfMissingFeatures('ELECTROSTATICS')
    def test_current(self):
        obs_data = espressomd.observables.Current(
            ids=self.system.part[:].id).calculate()
        part_data = self.system.part[:].q.dot(self.system.part[:].v)
        self.assertEqual(obs_data.shape, part_data.shape)
        np.testing.assert_array_almost_equal(
            obs_data, part_data, err_msg="Data did not agree for observable 'Current'", decimal=9)

    @utx.skipIfMissingFeatures('ELECTROSTATICS')
    def test_dipolemoment(self):
        obs = espressomd.observables.DipoleMoment(ids=self.system.part[:].id)
        obs_data = obs.calculate()
        part_data = self.system.part[:].q.dot(self.system.part[:].pos)
        self.assertEqual(obs_data.shape, part_data.shape)
        np.testing.assert_array_almost_equal(
            obs_data, part_data, err_msg="Data did not agree for observable 'DipoleMoment'", decimal=9)

    def test_com_force(self):
        id_list = sorted(
            np.random.choice(
                self.system.part[:].id,
                size=int(
                    self.N_PART * .9),
                replace=False))

        particles = self.system.part.select(
            lambda p: p.id in id_list and not p.virtual)

        np.testing.assert_allclose(
            np.sum(particles.f, axis=0),
            espressomd.observables.TotalForce(ids=id_list).calculate())


if __name__ == "__main__":
    ut.main()
