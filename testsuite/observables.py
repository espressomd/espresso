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
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])

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
            np.testing.assert_array_almost_equal(
                obs_data,
                part_data, err_msg="Data did not agree for observable " +
                str(obs_name) +
                " and particle property " +
                pprop_name, decimal=9)

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
        np.testing.assert_array_almost_equal(
            s,
            obs_data,
            err_msg="Stress tensor from analysis and observable did not agree",
            decimal=9)

    def test_com_position(self):
        if espressomd.has_features(["MASS"]):
            com = sum(
                (self.es.part[:].mass * self.es.part[:].pos.T).T, 0) / sum(self.es.part[:].mass)
        else:
            com = sum((self.es.part[:].pos.T).T, 0) / len(self.es.part)

        obs_data = ComPosition(ids=range(1000)).calculate()
        np.testing.assert_array_almost_equal(
            com, obs_data, err_msg="Center of mass observable wrong value", decimal=9)

    def test_com_velocity(self):
        if espressomd.has_features(["MASS"]):
            com_vel = sum(
                (self.es.part[:].mass * self.es.part[:].v.T).T, 0) / sum(self.es.part[:].mass)
        else:
            com_vel = sum((self.es.part[:].v.T).T, 0) / len(self.es.part)
        obs_data = ComVelocity(ids=range(1000)).calculate()
        np.testing.assert_array_almost_equal(
            com_vel,
            obs_data,
            err_msg="Center of mass velocity observable wrong value",
            decimal=9)


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
