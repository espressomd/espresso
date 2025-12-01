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
import espressomd.polymer
import espressomd.propagation
import espressomd.observables
Propagation = espressomd.propagation.Propagation


def generate_random_unit_vectors(N_PART):
    z = np.random.uniform(-1., 1., N_PART)
    r = np.sqrt(1. - z**2)
    phi = np.random.uniform(0., 2. * np.pi, N_PART)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.column_stack((x, y, z))


@utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH", "EXTERNAL_FORCES"])
class Test(ut.TestCase):
    """
    Check the total dipole field for a magnetic LJ fluid (500 particles,
    density approx 0.002, mu^2=1, no PBC).
    """
    # analytical solution for a ferrofluid in the thermal Stoner-Wohlfarth
    # model. Obtained from Eq.17 in :cite:`mostarac25a`.
    res_dict_fluid = {3.4283694213261087: 0.91,
                      1.1427898071087026: 0.62, 0.28569745177717565: 0.2}
    # analytical solution for a solid superparamagnet in the thermal
    # Stoner-Wohlfarth model. Obtained from Eq.15 in :cite:`mostarac25a`.
    res_dict_solid = {3.4283694213261087: 0.8,
                      1.1427898071087026: 0.54, 0.28569745177717565: 0.19}
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
    n_part = 100
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

    def _init_particles(self):
        system = self.system
        self.system.part.clear()
        orientor_list = generate_random_unit_vectors(N_PART=self.n_part)
        dip_mom_list = self.dip_reduced * orientor_list
        positions = espressomd.polymer.linear_polymer_positions(
            n_polymers=self.n_part, beads_per_chain=1, min_distance=1., bond_length=1., seed=self.seed)
        positions = np.reshape(positions, (-1, 3))
        particles = system.part.add(pos=positions, director=orientor_list)
        particles.rotation = (True, True, True)
        for p1, dipm_el in zip(list(particles), dip_mom_list):
            p2 = system.part.add(
                pos=p1.pos, dip=dipm_el, rotation=[False, False, False],
                magnetodynamics=self.default_magnetodynamics)
            p2.vs_auto_relate_to(p1)
            p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT

    def _apply_single_field_z_axis(self, h_reduced):
        for x in self.system.constraints:
            self.system.constraints.remove(x)
        ExtH = espressomd.constraints.HomogeneousMagneticField(
            H=(0., 0., h_reduced))
        self.system.constraints.add(ExtH)

    def _measure_dipole_moment(self, steps):
        dipm_tot = espressomd.observables.MagneticDipoleMoment(
            ids=self.system.part.select(lambda p: p.magnetodynamics["is_enabled"]).id)
        norm = 1 / (self.dip_reduced * self.n_part)
        self.system.integrator.run(steps)
        mag_el = dipm_tot.calculate() * norm
        return mag_el[-1]

    def test_tSW_fluid(self):
        STEPS = 12477
        self.n_part = 100
        for h_reduced, res in self.res_dict_fluid.items():
            self._init_particles()
            self._apply_single_field_z_axis(h_reduced)
            self.assertAlmostEqual(
                self._measure_dipole_moment(STEPS), res, delta=self.error)

    def test_tSW_solid(self):
        STEPS = 3447
        self.n_part = 500
        system = self.system
        for h_reduced, res in self.res_dict_solid.items():
            self._init_particles()
            part_slice = system.part.select(lambda p: p.is_virtual() == False)
            part_slice.rotation = [False, False, False]
            part_slice.fix = [True, True, True]
            self._apply_single_field_z_axis(h_reduced)
            self.assertAlmostEqual(
                self._measure_dipole_moment(STEPS), res, delta=self.error)


if __name__ == "__main__":
    ut.main()
