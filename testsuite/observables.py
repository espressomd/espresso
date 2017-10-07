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
from espressomd.observables import *


class Observables(ut.TestCase):

    # Error tolerance when comparing arrays/tuples...
    tol = 1E-9

    # Handle for espresso system
    es = espressomd.System()

    def arraysNearlyEqual(self, a, b):
        """Test, if the magnitude of the difference between two arrays is smaller than the tolerance"""

        # Check length
        if len(a) != len(b):
            return False

        # We have to use a loop, since we can't be sure, we're getting numpy
        # arrays
        sum = 0.
        for i in range(len(a)):
            sum += abs(a[i] - b[i])

        if sum > self.tol:
            return False

        return True

    def setUp(self):
        if not len(self.es.part):
            for i in range(1000):
                self.es.part.add(pos=random(3), v=random(3), id=i)
                if espressomd.has_features(["MASS"]):
                    self.es.part[i].mass = random()
                if espressomd.has_features(["DIPOLES"]):
                    self.es.part[i].dip = random(3)
                if espressomd.has_features(["ROTATION"]):
                    self.es.part[i].omega_lab = random(3)

    def generate_test_for_pid_observable(
            _obs_name, _pprop_name, _agg_type=None):
        """Generates test cases for observables working on particle id lists"""

        pprop_name = _pprop_name
        obs_name = _obs_name
        agg_type = _agg_type

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            # Get data from particles
            id_list = range(100, 500, 2)
            part_data = getattr(self.es.part[id_list], pprop_name)
            # Reshape and aggregate to linear array
            if len(part_data.shape) > 1:
                if (agg_type == "average"):
                    part_data = average(part_data, 0)
                if (agg_type == "sum"):
                    part_data = sum(part_data, 0)

                part_data = part_data.reshape(part_data.size)

            # Data from observable
            obs_data = obs_name(ids=id_list).calculate()
            self.assertTrue(
                self.arraysNearlyEqual(
                    obs_data,
                    part_data),
                "Data did not agree for observable " +
                str(obs_name) +
                " and particle property " +
                pprop_name)

        return func

    test_pos = generate_test_for_pid_observable(ParticlePositions, "pos")
    test_v = generate_test_for_pid_observable(ParticleVelocities, "v")
    test_f = generate_test_for_pid_observable(ParticleForces, "f")

    com_force = generate_test_for_pid_observable(ComForce, "f", "sum")

    if espressomd.has_features(["DIPOLES"]):
        test_mag_dip = generate_test_for_pid_observable(
            MagneticDipoleMoment, "dip", "sum")

    # This is disabled as it does not currently work
    # if espressomd.has_features(["ROTATION"]):
    #    test_omega_body = generate_test_for_pid_observable(ParticleBodyVelocities,"omega_body")

    def test_stress_tensor(self):
        s = self.es.analysis.stress_tensor()["total"].reshape(9)
        obs_data = np.array(StressTensor().calculate())
        self.assertTrue(
            self.arraysNearlyEqual(
                s,
                obs_data),
            "Stress tensor from analysis and observable did not agree")

    def test_stress_tensor_acf(self):
        s = self.es.analysis.stress_tensor()["total"].reshape(9)
        s = np.array((s[1], s[5], s[6], s[0] - s[4], s[0] - s[8], s[4] - s[8]))
        obs_data = np.array(StressTensorAcf().calculate())
        self.assertTrue(
            self.arraysNearlyEqual(
                s,
                obs_data),
            "Stress tensor from analysis and observable StressTensorAcf did not agree")

    def test_com_position(self):
        if espressomd.has_features(["MASS"]):
            com = sum(
                (self.es.part[:].mass * self.es.part[:].pos.T).T, 0) / sum(self.es.part[:].mass)
        else:
            com = sum((self.es.part[:].pos.T).T, 0) / len(self.es.part)

        obs_data = ComPosition(ids=range(1000)).calculate()
        self.assertTrue(self.arraysNearlyEqual(com, obs_data),
                        "Center of mass observable wrong value")

    def test_com_velocity(self):
        if espressomd.has_features(["MASS"]):
            com_vel = sum(
                (self.es.part[:].mass * self.es.part[:].v.T).T, 0) / sum(self.es.part[:].mass)
        else:
            com_vel = sum((self.es.part[:].v.T).T, 0) / len(self.es.part)
        obs_data = ComVelocity(ids=range(1000)).calculate()
        self.assertTrue(self.arraysNearlyEqual(com_vel, obs_data),
                        "Center of mass velocity observable wrong value")


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
