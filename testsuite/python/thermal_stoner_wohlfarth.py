#
# Copyright (C) 2025 The ESPResSo project
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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.propagation
Propagation = espressomd.propagation.Propagation


@utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH", "EXTERNAL_FORCES"])
class Test(ut.TestCase):
    # parameters taken from thermal_stoner_wohlfarth_fluid.py
    system = espressomd.System(box_l=(29.69314567, 29.69314567, 29.69314567))
    skin = 0.4
    seed = 42
    np.random.seed(seed)
    time_step = 0.001
    temperature = 1.
    # ani_energy = K1 * V, where kT_KVm_inv was previously ani_param = ani_energy/kT
    # So ani_energy = kT_KVm_inv * kT
    kT_KVm_inv = 5.  # old ani_param value
    ani_energy = temperature * kT_KVm_inv
    dt_incr = 0.001 * 3.437060795580368e-08
    HK_inv = 0.17501031139401407
    dip_reduced = 1.7501031139401464
    gamma_T = 74.86576383782938
    gamma_R = 24.955254612609792
    tau0_inv = 735412234.8230474
    error = 0.045
    default_magnetodynamics = {
        "is_enabled": True, "anisotropy_field_inv": HK_inv,
        "sat_mag": dip_reduced, "anisotropy_energy": ani_energy,
        "sw_dt_incr": dt_incr, "sw_tau0_inv": tau0_inv}

    def setUp(self):
        system = self.system
        system.cell_system.skin = 0.4
        system.min_global_cut = 1.
        system.time_step = 0.001
        system.periodicity = [True, True, True]
        system.thermostat.set_langevin(kT=self.temperature, gamma=self.gamma_T,
                                       gamma_rotation=self.gamma_R, seed=self.seed)

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        for x in self.system.constraints:
            self.system.constraints.remove(x)

    def _init_virtual_site_pair(self):
        self.system.part.clear()
        p1 = self.system.part.add(pos=[0, 0, 0], director=[1, 0, 0],
                                  rotation=[False, False, False],
                                  fix=[True, True, True])
        p2 = self.system.part.add(
            pos=p1.pos, dip=[1, 2, 3], rotation=[False, False, False],
            magnetodynamics=self.default_magnetodynamics)
        p2.vs_auto_relate_to(p1)
        p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT
        return p1, p2

    def _find_phi_minima(self, p2, max_iterations=1000):
        # for h=0.5 and psi=90 degrees, there are two minima at 60 and 300
        # degrees
        found_min1, found_min2 = False, False
        count = 0
        while not (found_min1 and found_min2) and count < max_iterations:
            self.system.integrator.run(100)
            phi0_deg = np.degrees(p2.magnetodynamics["sw_phi_0"])
            if np.isclose(phi0_deg, 60., atol=1e-06):
                found_min1 = True
            if np.isclose(phi0_deg, 300., atol=1e-06):
                found_min2 = True
            count += 1
        return found_min1, found_min2

    def _check_zero_field_flips(self, p2, phi0_start, max_iterations=10000):
        phi_no_flip, phi_yes_flip = False, False
        count = 0
        while not (phi_no_flip and phi_yes_flip) and count < max_iterations:
            old_dip = np.copy(p2.dip)
            self.system.integrator.run(1)
            new_phi = p2.magnetodynamics["sw_phi_0"]
            new_dip = np.copy(p2.dip)
            if phi0_start == new_phi and not phi_no_flip:
                np.testing.assert_allclose(old_dip, new_dip, atol=1e-06)
                phi_no_flip = True
            elif phi0_start != new_phi and not phi_yes_flip:
                np.testing.assert_allclose(-1 * old_dip, new_dip, atol=1e-06)
                phi_yes_flip = True

        return phi_no_flip, phi_yes_flip

    def _apply_single_field_z_axis(self, h_reduced):
        for x in self.system.constraints:
            self.system.constraints.remove(x)
        ExtH = espressomd.constraints.HomogeneousMagneticField(
            H=(0., 0., h_reduced))
        self.system.constraints.add(ExtH)

    def test_minimal_no_field(self):
        p1, p2 = self._init_virtual_site_pair()
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p1.director), np.copy(p2.director), atol=1e-06)
        test_flags = self._check_zero_field_flips(p2, 0.)
        self.assertEqual(all(test_flags), True)
        test_flags = self._check_zero_field_flips(p2, np.pi)
        self.assertEqual(all(test_flags), True)

    def test_minimal_field(self):
        _, p2 = self._init_virtual_site_pair()
        self.system.integrator.run(1)
        # critical reduced field i.e. h=1
        self._apply_single_field_z_axis(6.)
        self.system.integrator.run(0, recalc_forces=True)
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p2.director), np.array([0, 0, 1]), atol=1e-06)
        # reduced field h=0.5
        self._apply_single_field_z_axis(2.8569745177717567)
        found_min1, found_min2 = self._find_phi_minima(p2)
        self.assertEqual(found_min1, True)
        self.assertEqual(found_min2, True)


if __name__ == "__main__":
    ut.main()
