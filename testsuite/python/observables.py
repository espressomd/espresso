#
# Copyright (C) 2013-2022 The ESPResSo project
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
    partcls = system.part.by_ids(id_list)
    masses = partcls.mass

    # Filter out virtual particles by using mass=0 for them
    virtual = partcls.virtual
    masses[np.nonzero(virtual)] = 0.

    return np.average(getattr(partcls, x), weights=masses, axis=0)


class Observables(ut.TestCase):
    N_PART = 200
    # Handle for espresso system
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    partcls = system.part.add(
        id=np.arange(3, 3 + 2 * N_PART, 2),
        pos=np.random.random((N_PART, 3)) * system.box_l,
        v=np.random.random((N_PART, 3)) * 3.2 - 1,
        f=np.random.random((N_PART, 3)))

    if espressomd.has_features(["MASS"]):
        partcls.mass = np.random.random(N_PART)

    if espressomd.has_features(["DIPOLES"]):
        partcls.dip = np.random.random((N_PART, 3)) - .3

    if espressomd.has_features(["DIPOLE_FIELD_TRACKING"]):
        system.analysis.dipole_fields()

    if espressomd.has_features(["ROTATION"]):
        partcls.omega_body = np.random.random((N_PART, 3)) - .5
        partcls.torque_lab = np.random.random((N_PART, 3)) - .5
        direcs = np.random.random((N_PART, 3)) - 0.5
        direcs /= np.linalg.norm(direcs, axis=1)[:, None]
        partcls.director = direcs

    if espressomd.has_features("DIPOLES"):
        partcls.dipm = np.random.random(N_PART) + 2

    if espressomd.has_features("ELECTROSTATICS"):
        partcls.q = np.random.random(N_PART)

    if espressomd.has_features("VIRTUAL_SITES"):
        p = system.part.by_id(partcls.id[8])
        p.virtual = True

    def generate_test_for_pid_observable(
            _obs_class, _pprop_name, _agg_type=None):
        """Generates test cases for observables working on particle id lists.

        """
        pprop_name = _pprop_name
        obs_class = _obs_class
        agg_type = _agg_type

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            # Randomly pick a subset of the particles
            id_list = sorted(
                np.random.choice(
                    self.system.part.all().id,
                    size=int(
                        self.N_PART * .9),
                    replace=False))
            for id in id_list:
                self.assertTrue(self.system.part.exists(id))

            # Get data from particles
            if pprop_name == "f":
                for p_id in id_list:
                    if self.system.part.by_id(p_id).virtual:
                        id_list.remove(p_id)

            part_data = getattr(self.system.part.by_ids(id_list), pprop_name)

            # Reshape and aggregate to linear array
            if len(part_data.shape) > 1:
                if agg_type == "sum":
                    part_data = np.sum(part_data, 0)
                if agg_type == 'com':
                    part_data = calc_com_x(self.system, pprop_name, id_list)

            # Data from observable
            observable = obs_class(ids=id_list)
            obs_data = observable.calculate()

            # Check
            self.assertEqual(obs_data.shape, part_data.shape)
            np.testing.assert_equal(id_list, observable.ids)

            np.testing.assert_array_almost_equal(
                obs_data,
                part_data,
                err_msg=f"Data did not agree for observable {obs_class.__name__} and particle property {pprop_name}",
                decimal=11)

            # Test setters and getters
            self.assertEqual(observable.ids, id_list)
            with self.assertRaises(RuntimeError):
                observable.ids = [observable.ids[0]]

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

    if espressomd.has_features(["DIPOLE_FIELD_TRACKING"]):
        test_dip_fld = generate_test_for_pid_observable(
            espressomd.observables.ParticleDipoleFields, "dip_fld")

    if espressomd.has_features(["ROTATION"]):
        test_body_angular_velocity = generate_test_for_pid_observable(
            espressomd.observables.ParticleBodyAngularVelocities, "omega_body")
        test_lab_angular_velocity = generate_test_for_pid_observable(
            espressomd.observables.ParticleAngularVelocities, "omega_lab")
        test_director = generate_test_for_pid_observable(
            espressomd.observables.ParticleDirectors, "director")

    @ut.skipIf(espressomd.has_features(["DIPOLE_FIELD_TRACKING"]),
               "default dipole fields are needed")
    def test_director_no_dipole_fields(self):
        id_list = self.system.part.all().id
        observable = espressomd.observables.ParticleDipoleFields(ids=id_list)
        obs_data = observable.calculate()
        np.testing.assert_array_almost_equal(
            obs_data, self.N_PART * [[0., 0., 0.]], decimal=11)

    @ut.skipIf(espressomd.has_features(["ROTATION"]),
               "check default directors")
    def test_director_norotation(self):
        id_list = self.system.part.all().id
        observable = espressomd.observables.ParticleDirectors(ids=id_list)
        obs_data = observable.calculate()
        np.testing.assert_array_almost_equal(
            obs_data, self.N_PART * [[0., 0., 1.]], decimal=11)

    @utx.skipIfMissingFeatures(['ROTATION'])
    def test_particle_body_velocities(self):
        obs = espressomd.observables.ParticleBodyVelocities(
            ids=self.system.part.all().id)
        obs_data = obs.calculate()
        part_data = np.array([p.convert_vector_space_to_body(p.v)
                              for p in self.system.part])
        self.assertEqual(obs_data.shape, part_data.shape)
        np.testing.assert_array_almost_equal(part_data, obs_data,
                                             err_msg="Data did not agree for observable ParticleBodyVelocities and particle derived values.",
                                             decimal=9)

    def test_energy(self):
        s = self.system.analysis.energy()["total"]
        obs_data = espressomd.observables.Energy().calculate()
        self.assertEqual(obs_data.shape, (1,))
        np.testing.assert_array_almost_equal(
            obs_data,
            s,
            err_msg="Energy from analysis and observable did not agree",
            decimal=9)

    def test_pressure(self):
        s = self.system.analysis.pressure()["total"]
        obs_data = espressomd.observables.Pressure().calculate()
        self.assertEqual(obs_data.shape, (1,))
        np.testing.assert_array_almost_equal(
            obs_data,
            s,
            err_msg="Pressure from analysis and observable did not agree",
            decimal=9)

    def test_pressure_tensor(self):
        s = self.system.analysis.pressure_tensor()["total"]
        obs_data = espressomd.observables.PressureTensor().calculate()
        self.assertEqual(obs_data.shape, s.shape)
        np.testing.assert_array_almost_equal(
            obs_data,
            s,
            err_msg="Pressure tensor from analysis and observable did not agree",
            decimal=9)

    @utx.skipIfMissingFeatures('ELECTROSTATICS')
    def test_dipolemoment(self):
        obs = espressomd.observables.DipoleMoment(ids=self.partcls.id)
        obs_data = obs.calculate()
        part_data = self.partcls.q.dot(self.partcls.pos)
        self.assertEqual(obs_data.shape, part_data.shape)
        np.testing.assert_array_almost_equal(
            obs_data, part_data, err_msg="Data did not agree for observable 'DipoleMoment'", decimal=9)

    def test_com_force(self):
        id_list = sorted(
            np.random.choice(
                self.partcls.id,
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
