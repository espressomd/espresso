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


@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS"]),
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
        self.system.time_step = 0.01
        self.kT = 0.5
        self.system.thermostat.set_langevin(kT=self.kT, gamma=1.0)
        self.system.seed = range(
            self.system.cell_system.get_state()["n_nodes"])
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1.0 / 6.0), shift="auto")
        num_part = 40
        mass = 1
        np.random.seed(seed=1)

        for i in range(num_part):
            self.system.part.add(pos=np.random.random(
                3) * self.system.box_l, q=1, v=np.sqrt(self.kT / mass) * np.random.normal(loc=[0, 0, 0]))
            self.system.part.add(pos=np.random.random(
                3) * self.system.box_l, q=-1, v=np.sqrt(self.kT / mass) * np.random.normal(loc=[0, 0, 0]))

    def pressure_via_volume_scaling(
        self,
        system,
     kbT,
     list_of_previous_values):
        # taken from "Efficient pressure estimation in molecular simulations without evaluating the virial"
        # only works so far for isotropic volume changes, i.e. the isotropic
        # pressure
        energy = system.analysis.energy()
        Epot_old = energy["total"] - energy["kinetic"]
        old_box_lengths = system.box_l
        old_volume = np.prod(old_box_lengths)
        dV_div_old_volume = 0.001
        dV = -dV_div_old_volume * old_volume
        new_volume = old_volume + dV
        new_box_l = (new_volume)**(1. / 3.)
        system.change_volume_and_rescale_particles(new_box_l, "xyz")
        system.integrator.run(0)
        energy = system.analysis.energy()
        Epot_new = energy["total"] - energy["kinetic"]
        system.change_volume_and_rescale_particles(old_box_lengths[0], "xyz")
        system.integrator.run(0)
        DeltaEpot = Epot_new - Epot_old
        particle_number = len(system.part[:].id)
        current_value = (new_volume / old_volume)**particle_number * \
            np.exp(-DeltaEpot / kbT)
        list_of_previous_values.append(current_value)
        average_value = np.mean(list_of_previous_values)

        pressure = kbT / dV * np.log(average_value)
        return pressure

    def test_p3m_pressure(self):
        pressures_via_virial = []
        pressures_via_volume_scaling = []
        p3m = electrostatics.P3M(prefactor=2.0, accuracy=1e-3)
        self.system.actors.add(p3m)
        num_samples = 100
        pressure_via_volume_scaling = np.nan
        for i in range(num_samples):
            self.system.integrator.run(100)
            pressures_via_virial.append(
                self.system.analysis.pressure()['total'])
            pressure_via_volume_scaling = self.pressure_via_volume_scaling(
                self.system, self.kT, pressures_via_volume_scaling)
        pressure_virial = np.mean(pressures_via_virial)
        ratio_pressure_virial_pressure_volume_scaling = pressure_via_volume_scaling / \
            pressure_virial  # should be 1 ideally
        self.assertAlmostEqual(
            ratio_pressure_virial_pressure_volume_scaling-1.0, 0.0, places=2,
                               msg="Difference to between isotropic virial pressure and pressure via volume derivative of potential energy is too big. The result must be self-consistent.")


if __name__ == "__main__":
    ut.main()
