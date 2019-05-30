#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import numpy as np
import numpy.testing as npt

import espressomd
from espressomd import electrostatics


class pressureViaVolumeScaling(object):

    def __init__(self, system, kbT):
        self.system = system
        self.kbT = kbT
        self.old_box_lengths = np.copy(system.box_l)
        self.old_volume = np.prod(self.old_box_lengths)
        dV_div_old_volume = 0.001
        self.dV = -dV_div_old_volume * self.old_volume
        self.new_volume = self.old_volume + self.dV
        self.new_box_l = (self.new_volume)**(1. / 3.)

        self.list_of_previous_values = []
    
    def measure_pressure_via_volume_scaling(
            self):
        # taken from "Efficient pressure estimation in molecular simulations without evaluating the virial"
        # only works so far for isotropic volume changes, i.e. the isotropic
        # pressure
        energy = self.system.analysis.energy()
        Epot_old = energy["total"] - energy["kinetic"]
        self.system.change_volume_and_rescale_particles(self.new_box_l, "xyz")
        self.system.integrator.run(0)
        energy = self.system.analysis.energy()
        Epot_new = energy["total"] - energy["kinetic"]
        self.system.change_volume_and_rescale_particles(
            self.old_box_lengths[0], "xyz")
        self.system.integrator.run(0)
        DeltaEpot = Epot_new - Epot_old
        particle_number = len(self.system.part[:].id)
        current_value = (self.new_volume / self.old_volume)**particle_number * \
            np.exp(-DeltaEpot / self.kbT)
        self.list_of_previous_values.append(current_value)
    
    def get_result(self):
        average_value = np.mean(self.list_of_previous_values)

        pressure = self.kbT / self.dV * np.log(average_value)
        return pressure


@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class VirialPressureConsistency(ut.TestCase):

    """Test the consistency of the core implementation of the virial pressure with an analytical relation which allows
       for the calculation of the pressure as a volume derivative of a function of the potential energy change on infinitesimal volume changes.
       The relation and its derivation can be found in the paper with the name "Efficient pressure estimation in molecular simulations without evaluating the virial"
       by Harismiadis, V. I., J. Vorholz, and A. Z. Panagiotopoulos. 1996

    """
    # Handle to espresso system
    system = espressomd.System(box_l=[50, 50, 50])

    def setUp(self):
        np.random.seed(seed=1)
        self.system.seed = range(
            self.system.cell_system.get_state()["n_nodes"])
        self.system.time_step = 0.01
        self.kT = 0.5
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1.0 / 6.0), shift="auto")
        num_part = 40
        mass = 1
        
        for i in range(num_part):
            self.system.part.add(pos=np.random.random(
                3) * self.system.box_l, q=1, v=np.sqrt(self.kT / mass) * np.random.normal(loc=[0, 0, 0]))
            self.system.part.add(pos=np.random.random(
                3) * self.system.box_l, q=-1, v=np.sqrt(self.kT / mass) * np.random.normal(loc=[0, 0, 0]))

        #############################################################
        #  Warmup Integration                                       #
        #############################################################

        self.system.integrator.set_steepest_descent(
            f_max=0,
            gamma=0.001,
            max_displacement=0.01)

        # warmup
        while self.system.analysis.energy()["total"] > 10 * num_part:
            print("minimization: {:.1f}".format(
                self.system.analysis.energy()["total"]))
            self.system.integrator.run(10)
        self.system.integrator.set_vv()
        self.system.thermostat.set_langevin(kT=self.kT, gamma=1.0, seed=41)

    def test_p3m_pressure(self):
        pressures_via_virial = []
        pressures_via_volume_scaling = []
        p3m = electrostatics.P3M(
            prefactor=2.0,
            accuracy=1e-3,
            mesh=16,
            cao=6,
            r_cut=1.4941e-01 * self.system.box_l[0])
        self.system.actors.add(p3m)
        print("Tune skin: {}".format(self.system.cell_system.tune_skin(
                                     min_skin=0.0, max_skin=2.5, tol=0.05, int_steps=100)))

        num_samples = 100
        pressure_via_volume_scaling = pressureViaVolumeScaling(
            self.system, self.kT)
        for i in range(num_samples):
            self.system.integrator.run(100)
            pressures_via_virial.append(
                self.system.analysis.pressure()['total'])
            pressure_via_volume_scaling.measure_pressure_via_volume_scaling()
        pressure_virial = np.mean(pressures_via_virial)
        abs_deviation_in_percent = abs(
            pressure_virial / pressure_via_volume_scaling.get_result() - 1.0) * 100.0  # should be 0% ideally
        npt.assert_array_less(
            abs_deviation_in_percent,
            5.0)  # devation should be below 5%


if __name__ == "__main__":
    ut.main()
