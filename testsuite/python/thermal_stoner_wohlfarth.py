#
# Copyright (C) 2023 The ESPResSo project
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
import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd.polymer
import espressomd.propagation
Propagation = espressomd.propagation.Propagation
from espressomd.observables import MagneticDipoleMoment


def generate_random_unit_vectors(N_PART):
    z = np.random.uniform(-1, 1, N_PART)
    r = np.sqrt(1 - z * z)
    phi = np.random.uniform(0, 2 * np.pi, N_PART)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.column_stack((x, y, z))


@utx.skipIfMissingFeatures(["NLOPT"])
class Test(ut.TestCase):
    """
    Check the total dipole field for a magnetic LJ fluid (500 particles,
    density approx 0.002, mu^2=1, no PBC).
    """
    # Values coorespond to analytical solution for a ferrofluid in the thermal Stoner-Wohlfarth model. Obtained from Eq.17 in https://doi.org/10.1103/PhysRevB.111.014438.
    res_dict_fluid = {3.4283694213261087: 0.91,
                      1.1427898071087026: 0.62, 0.28569745177717565: 0.2}
    # Values coorespond to analytical solution for a solid superparamagnet  in the thermal Stoner-Wohlfarth model. Obtained from Eq.15 in https://doi.org/10.1103/PhysRevB.111.014438.
    res_dict_solid = {3.4283694213261087: 0.8,
                      1.1427898071087026: 0.54, 0.28569745177717565: 0.19}
    system = espressomd.System(box_l=(29.69314567, 29.69314567, 29.69314567))
    skin = 0.4
    seed = 42
    np.random.seed(seed)
    time_step = 0.001
    temperature = 1
    # ani_energy = K1 * V, where kT_KVm_inv was previously ani_param = ani_energy/kT
    # So ani_energy = kT_KVm_inv * kT
    kT_KVm_inv = 5  # old ani_param value
    ani_energy = temperature * kT_KVm_inv
    dt_incr = 0.001 * 3.437060795580368e-08
    HK_inv = 0.17501031139401407
    dip_reduced = 1.7501031139401464
    gamma_T = 74.86576383782938
    gamma_R = 24.955254612609792
    tau0_inv = 735412234.8230474
    n_part = 100
    error = 0.035

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

    def _init_virtual_site_pair(self):
        self.system.part.clear()
        p1 = self.system.part.add(pos=[0, 0, 0], director=[1, 0, 0])
        p1.rotation = (False, False, False)
        p1.fix = (True, True, True)
        p2 = self.system.part.add(
            pos=p1.pos, dip=[1, 2, 3], rotation=[False, False, False], magnetodynamics={'is_enabled': True, 'anisotropy_field_inv': self.HK_inv, 'sat_mag': self.dip_reduced, 'anisotropy_energy': self.ani_energy, 'sw_dt_incr': self.dt_incr, 'sw_tau0_inv': self.tau0_inv})
        p2.vs_auto_relate_to(p1)
        p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT
        return p1, p2
    # for h=0.5 and psi=90 degrees, there are two minima, at 60 and 300 degrees respectively.

    def _find_phi_minima(self, p2, max_iterations=1000):
        found_min1, found_min2 = False, False
        count = 0
        while (not found_min1 or not found_min2) and count < max_iterations:
            self.system.integrator.run(100)
            phi0_deg = np.degrees(p2.magnetodynamics['sw_phi_0'])
            if np.isclose(phi0_deg, 60., atol=1e-06):
                found_min1 = True
            if np.isclose(phi0_deg, 300., atol=1e-06):
                found_min2 = True
            count += 1
        return found_min1, found_min2

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
                pos=p1.pos, dip=dipm_el, rotation=[False, False, False], magnetodynamics={'is_enabled': True, 'anisotropy_field_inv': self.HK_inv, 'sat_mag': self.dip_reduced, 'anisotropy_energy': self.ani_energy, 'sw_dt_incr': self.dt_incr, 'sw_tau0_inv': self.tau0_inv})
            p2.vs_auto_relate_to(p1)
            p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT

    def _apply_single_field_z_axis(self, h_reduced):
        for x in self.system.constraints:
            self.system.constraints.remove(x)
        ExtH = espressomd.constraints.HomogeneousMagneticField(
            H=(0, 0, h_reduced))
        self.system.constraints.add(ExtH)

    def _measure_dipole_moment(self, steps):
        dipm_tot = MagneticDipoleMoment(
            ids=self.system.part.select(lambda p: p.magnetodynamics['is_enabled'] == True).id)
        norm = 1 / (self.dip_reduced * self.n_part)
        self.system.integrator.run(steps)
        mag_el = dipm_tot.calculate() * norm
        return mag_el[-1]

    @utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH"])
    def test_minimal(self):
        p1, p2 = self._init_virtual_site_pair()
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p1.director), np.copy(p2.director), atol=1e-06)
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

    @utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH"])
    def test_tSW_fluid(self):

        STEPS = 12477
        self.n_part = 100
        for h_reduced, res in self.res_dict_fluid.items():
            self._init_particles()
            self._apply_single_field_z_axis(h_reduced)
            self.assertAlmostEqual(
                self._measure_dipole_moment(STEPS), res, delta=self.error)

    @utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH"])
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
