#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from numpy.random import random
from espressomd.interactions import FeneBond
import espressomd.observables

def calc_com_x(system, x):
    if espressomd.has_features(["MASS"]):
        com_x = np.average(getattr(system.part[:], x), weights=system.part[:].mass, axis=0)
    else:
        com_x = np.average(getattr(system.part[:], x), axis=0)
    return com_x

class Observables(ut.TestCase):
    N_PART = 1000
    # Handle for espresso system
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    
    def setUp(self):
        if not len(self.system.part):
            for i in range(self.N_PART):
                self.system.part.add(pos=random(3)*10, v=random(3), id=i)
                if espressomd.has_features(["MASS"]):
                    self.system.part[i].mass = random()
                if espressomd.has_features(["DIPOLES"]):
                    self.system.part[i].dip = random(3)
                if espressomd.has_features(["ROTATION"]):
                    self.system.part[i].omega_lab = random(3)
                if espressomd.has_features("ELECTROSTATICS"):
                    self.system.part[i].q = (1 if i % 2 == 0 else -1)

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
            # Get data from particles
            id_list = range(self.N_PART)
            part_data = getattr(self.system.part[id_list], pprop_name)
            # Reshape and aggregate to linear array
            if len(part_data.shape) > 1:
                if (agg_type == "average"):
                    part_data = average(part_data, 0)
                if (agg_type == "sum"):
                    part_data = sum(part_data, 0)
                if (agg_type == 'com'):
                    part_data = calc_com_x(self.system, pprop_name)
                part_data = part_data.flatten()

            # Data from observable
            observable = obs_name(ids=id_list)
            obs_data = observable.calculate()
            self.assertEqual(observable.n_values(), len(part_data))
            np.testing.assert_array_almost_equal(
                obs_data,
                part_data, err_msg="Data did not agree for observable " +
                str(obs_name) +
                " and particle property " +
                pprop_name, decimal=9)

        return func

    test_pos = generate_test_for_pid_observable(espressomd.observables.ParticlePositions, "pos")
    test_v = generate_test_for_pid_observable(espressomd.observables.ParticleVelocities, "v")
    test_f = generate_test_for_pid_observable(espressomd.observables.ParticleForces, "f")
    test_com_position = generate_test_for_pid_observable(espressomd.observables.ComPosition, 'pos', 'com')
    test_com_velocity = generate_test_for_pid_observable(espressomd.observables.ComVelocity, 'v', 'com')
    test_com_force = generate_test_for_pid_observable(espressomd.observables.ComForce, 'f', 'com')

    if espressomd.has_features(["DIPOLES"]):
        test_mag_dip = generate_test_for_pid_observable(
            espressomd.observables.MagneticDipoleMoment, "dip", "sum")

    if espressomd.has_features(["ROTATION"]):
        test_body_angular_velocity = generate_test_for_pid_observable(espressomd.observables.ParticleBodyAngularVelocities, "omega_body")
        test_lab_angular_velocity = generate_test_for_pid_observable(espressomd.observables.ParticleAngularVelocities, "omega_lab")

        def test_particle_body_velocities(self):
            obs = espressomd.observables.ParticleBodyVelocities(ids=range(self.N_PART))
            obs_data = obs.calculate()
            part_data = np.array([p.convert_vector_space_to_body(p.v) for p in self.system.part])
            np.testing.assert_array_almost_equal(part_data.flatten(), obs_data,
                err_msg="Data did not agree for observable ParticleBodyVelocities and particle derived values.",
                decimal=9)


    def test_stress_tensor(self):
        s = self.system.analysis.stress_tensor()["total"].reshape(9)
        obs_data = np.array(espressomd.observables.StressTensor().calculate())
        self.assertEqual(espressomd.observables.StressTensor().n_values(), len(s))
        np.testing.assert_array_almost_equal(
            s,
            obs_data,
            err_msg="Stress tensor from analysis and observable did not agree",
            decimal=9)

    @ut.skipIf(not espressomd.has_features('ELECTROSTATICS'), "Skipping test for Current observable due to missing features.")
    def test_current(self):
        obs_data = espressomd.observables.Current(ids=range(self.N_PART)).calculate()
        part_data = self.system.part[:].q.dot(self.system.part[:].v)
        self.assertEqual(espressomd.observables.Current(ids=range(self.N_PART)).n_values(), len(part_data.flatten()))
        np.testing.assert_array_almost_equal(obs_data, part_data, err_msg="Data did not agree for observable 'Current'", decimal=9)

    @ut.skipIf(not espressomd.has_features('ELECTROSTATICS'), "Skipping test for DipoleMoment observable due to missing features.")
    def test_dipolemoment(self):
        obs = espressomd.observables.DipoleMoment(ids=range(self.N_PART))
        obs_data = obs.calculate()
        part_data = self.system.part[:].q.dot(self.system.part[:].pos)
        self.assertEqual(obs.n_values(), len(part_data.flatten()))
        np.testing.assert_array_almost_equal(obs_data, part_data, err_msg="Data did not agree for observable 'DipoleMoment'", decimal=9)


if __name__ == "__main__":
    ut.main()
