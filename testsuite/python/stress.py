# Copyright (C) 2010-2019 The ESPResSo project
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
from espressomd.interactions import HarmonicBond
from espressomd.interactions import FeneBond
from espressomd.observables import StressTensor

from tests_common import fene_force2

import numpy as np

# allowed deviation from analytical results
tol = 1.0e-13


def stress_kinetic(vel):
    '''Analytical result for convective stress'''
    return np.einsum('ij,ik->jk', vel, vel) / np.prod(np.copy(system.box_l))


def stress_bonded(pos):
    '''Analytical result for stress originating from bonded forces'''
    stress = np.zeros([3, 3])
    for p1, p2 in zip(pos[0::2], pos[1::2]):
        r = p1 - p2
        f = -1.0e4 * r
        stress += np.einsum('i,j', f, r) / np.prod(np.copy(system.box_l))
    return stress


def stress_nonbonded(particle_pairs):
    '''Analytical result for stress originating from non-bonded forces'''
    stress = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if (p1.type == 0 and p2.type == 0) or (p1.type == 1 and p2.type == 2):
            d = p1.pos - p2.pos
            r = np.linalg.norm(d)
            r_hat = d / r
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / r**13 - 1.0**6 / r**7)) * r_hat
            stress += np.einsum('i,j', f, d) / np.prod(np.copy(system.box_l))
    return stress


def stress_nonbonded_inter(particle_pairs):
    stress = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if p1.type == 1 and p2.type == 2 and p1.mol_id != p2.mol_id:
            r = p1.pos - p2.pos
            d = np.linalg.norm(r)
            r_hat = r / d
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / d**13 - 1.0**6 / d**7)) * r_hat
            stress += np.einsum('i,j', f, r) / np.prod(np.copy(system.box_l))
    return stress


def stress_nonbonded_intra(particle_pairs):
    stress = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if p1.type == 0 and p2.type == 0 and p1.mol_id == p2.mol_id:
            r = p1.pos - p2.pos
            d = np.linalg.norm(r)
            r_hat = r / d
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / d**13 - 1.0**6 / d**7)) * r_hat
            stress += np.einsum('i,j', f, r) / np.prod(np.copy(system.box_l))
    return stress


system = espressomd.System(box_l=[1.0, 1.0, 1.0])


@utx.skipIfMissingFeatures(['LENNARD_JONES'])
class Stress(ut.TestCase):

    def test(self):
        # system parameters
        system.box_l = 3 * [10.0]
        skin = 0.4
        time_step = 0.01
        system.time_step = time_step

        # thermostat and cell system
        system.thermostat.set_langevin(kT=0.0, gamma=1.0, seed=41)
        system.cell_system.skin = skin
        system.periodicity = [1, 1, 1]

        # particles and bond
        system.part.add(id=0, pos=[9.9, 9.75, 9.9], type=0, mol_id=0)
        system.part.add(id=1, pos=[9.9, 10.25, 9.9], type=0, mol_id=0)
        system.part.add(id=2, pos=[0.1, 9.7, 0.1], type=1, mol_id=1)
        system.part.add(id=3, pos=[0.1, 10.3, 0.1], type=2, mol_id=2)

        harmonic = HarmonicBond(k=1e4, r_0=0)
        system.bonded_inter.add(harmonic)
        system.part[0].add_bond((harmonic, 1))
        system.part[2].add_bond((harmonic, 3))

        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.non_bonded_inter[1, 2].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)

        system.integrator.run(steps=0)

        system.part[0].v = [10.0, 20.0, 30.0]
        system.part[1].v = [-15, -25, -35]
        system.part[2].v = [27.0, 23.0, 17.0]
        system.part[3].v = [13.0, 11.0, 19.0]

        pos = system.part[:].pos
        vel = system.part[:].v

        sim_stress_kinetic = np.copy(
            system.analysis.stress_tensor()['kinetic'])
        sim_stress_bonded = np.copy(system.analysis.stress_tensor()['bonded'])
        sim_stress_bonded_harmonic = np.copy(system.analysis.stress_tensor()[
            'bonded', len(system.bonded_inter) - 1])
        sim_stress_nonbonded = np.copy(
            system.analysis.stress_tensor()['non_bonded'])
        sim_stress_nonbonded_inter = np.copy(system.analysis.stress_tensor()[
            'non_bonded_inter'])
        sim_stress_nonbonded_inter12 = np.copy(system.analysis.stress_tensor()[
            'non_bonded_inter', 1, 2])
        sim_stress_nonbonded_intra = np.copy(system.analysis.stress_tensor()[
            'non_bonded_intra'])
        sim_stress_nonbonded_intra00 = np.copy(system.analysis.stress_tensor()[
            'non_bonded_intra', 0, 0])
        sim_stress_total = np.copy(system.analysis.stress_tensor()['total'])
        sim_pressure_kinetic = system.analysis.pressure()['kinetic']
        sim_pressure_bonded = system.analysis.pressure()['bonded']
        sim_pressure_bonded_harmonic = system.analysis.pressure()[
            'bonded', len(system.bonded_inter) - 1]
        sim_pressure_nonbonded = system.analysis.pressure()['non_bonded']
        sim_pressure_nonbonded_inter = system.analysis.pressure()[
            'non_bonded_inter']
        sim_pressure_nonbonded_inter12 = system.analysis.pressure()[
            'non_bonded_inter', 1, 2]
        sim_pressure_nonbonded_intra = system.analysis.pressure()[
            'non_bonded_intra']
        sim_pressure_nonbonded_intra00 = system.analysis.pressure()[
            'non_bonded_intra', 0, 0]
        sim_pressure_total = system.analysis.pressure()['total']

        anal_stress_kinetic = stress_kinetic(vel)
        anal_stress_bonded = stress_bonded(pos)
        anal_stress_nonbonded = stress_nonbonded(system.part.pairs())
        anal_stress_nonbonded_inter = stress_nonbonded_inter(
            system.part.pairs())
        anal_stress_nonbonded_intra = stress_nonbonded_intra(
            system.part.pairs())
        anal_stress_total = anal_stress_kinetic + \
            anal_stress_bonded + anal_stress_nonbonded
        anal_pressure_kinetic = np.einsum('ii', anal_stress_kinetic) / 3.0
        anal_pressure_bonded = np.einsum('ii', anal_stress_bonded) / 3.0
        anal_pressure_nonbonded = np.einsum('ii', anal_stress_nonbonded) / 3.0
        anal_pressure_nonbonded_inter = np.einsum(
            'ii', anal_stress_nonbonded_inter) / 3.0
        anal_pressure_nonbonded_intra = np.einsum(
            'ii', anal_stress_nonbonded_intra) / 3.0
        anal_pressure_total = anal_pressure_kinetic + \
            anal_pressure_bonded + anal_pressure_nonbonded

        system.part.clear()

        np.testing.assert_allclose(
            sim_stress_kinetic, anal_stress_kinetic, rtol=0, atol=tol,
            err_msg='kinetic stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_bonded, anal_stress_bonded, rtol=0, atol=tol,
            err_msg='bonded stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_bonded_harmonic, anal_stress_bonded, rtol=0, atol=tol,
            err_msg='bonded stress harmonic bond does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_nonbonded, anal_stress_nonbonded, rtol=0, atol=tol,
            err_msg='non-bonded stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_nonbonded_inter, anal_stress_nonbonded_inter, rtol=0, atol=tol,
            err_msg='non-bonded intermolecular stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_nonbonded_inter12, anal_stress_nonbonded_inter, rtol=0, atol=tol,
            err_msg='non-bonded intermolecular stress molecules 1 and 2 does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_nonbonded_intra, anal_stress_nonbonded_intra, rtol=0, atol=tol,
            err_msg='non-bonded intramolecular stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_nonbonded_intra00, anal_stress_nonbonded_intra, rtol=0, atol=tol,
            err_msg='non-bonded intramolecular stress molecule 0 does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_total, anal_stress_total, rtol=0, atol=tol,
            err_msg='total stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_total, sim_stress_kinetic + sim_stress_bonded + sim_stress_nonbonded, rtol=0, atol=tol,
            err_msg='total stress is not given as the sum of all major stress components')
        self.assertAlmostEqual(
            sim_pressure_kinetic, anal_pressure_kinetic, delta=tol,
            msg='kinetic pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_bonded, anal_pressure_bonded, delta=tol,
            msg='bonded pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_bonded_harmonic, anal_pressure_bonded, delta=tol,
            msg='bonded pressure harmonic bond does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_nonbonded, anal_pressure_nonbonded, delta=tol,
            msg='non-bonded pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_nonbonded_inter, anal_pressure_nonbonded_inter, delta=tol,
            msg='non-bonded intermolecular pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_nonbonded_inter12, anal_pressure_nonbonded_inter, delta=tol,
            msg='non-bonded intermolecular pressure molecule 1 and 2 does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_nonbonded_intra, anal_pressure_nonbonded_intra, delta=tol,
            msg='non-bonded intramolecular pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_nonbonded_intra00, anal_pressure_nonbonded_intra, delta=tol,
            msg='non-bonded intramolecular pressure molecule 0 does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_total, anal_pressure_total, delta=tol,
            msg='total pressure does not match analytical result')
        self.assertAlmostEqual(
            sim_pressure_total, sim_pressure_kinetic + sim_pressure_bonded + sim_pressure_nonbonded, delta=tol,
            msg='total pressure is not given as the sum of all major pressure components')

        # Compare stress tensor observable to stress tensor from analysis
        np.testing.assert_allclose(
            StressTensor().calculate(),
            system.analysis.stress_tensor()["total"],
            atol=1E-10)


@utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
class StressFENE(ut.TestCase):

    def get_anal_stress_fene(self, pos_1, pos_2, k, d_r_max, r_0):
        stress = np.zeros([3, 3])
        vec_r = pos_1 - pos_2
        f = -fene_force2(vec_r, k, d_r_max, r_0)
        stress += np.einsum('i,j', f, vec_r) / np.prod(np.copy(system.box_l))
        return stress

    def test_fene(self):
        # system parameters
        system.box_l = 3 * [10.0]
        skin = 0.4
        time_step = 0.01
        system.time_step = time_step

        # thermostat and cell system
        system.cell_system.skin = skin
        system.periodicity = [1, 1, 1]

        # particles and bond
        system.part.add(
            id=0, pos=[9.9, 9.75, 9.9], type=0, mol_id=0, fix=[1, 1, 1])
        system.part.add(
            id=1, pos=[9.9, 10.25, 9.9], type=0, mol_id=0, fix=[1, 1, 1])

        k = 1e4
        d_r_max = 1.5
        r_0 = 0.1

        fene = FeneBond(k=k, d_r_max=d_r_max, r_0=r_0)
        system.bonded_inter.add(fene)
        system.part[0].add_bond((fene, 1))
        system.integrator.run(steps=0)

        sim_stress_bonded = system.analysis.stress_tensor()['bonded']
        sim_stress_fene = system.analysis.stress_tensor()[
            'bonded', len(system.bonded_inter) - 1]

        total_bonded_stresses = np.zeros([3, 3])
        for i in range(len(system.bonded_inter)):
            total_bonded_stresses += system.analysis.stress_tensor()[
                'bonded', i]

        anal_stress_fene = self.get_anal_stress_fene(
            system.part[0].pos, system.part[1].pos, k, d_r_max, r_0)
        np.testing.assert_allclose(
            sim_stress_bonded, anal_stress_fene, atol=tol,
            err_msg='bonded stress does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_fene, anal_stress_fene, atol=tol,
            err_msg='bonded stress for fene does not match analytical result')
        np.testing.assert_allclose(
            sim_stress_bonded, total_bonded_stresses, atol=tol,
            err_msg='bonded stresses do not sum up to the total value')

        sim_pressure_fene = system.analysis.pressure()[
            'bonded', len(system.bonded_inter) - 1]
        anal_pressure_fene = np.einsum("ii", anal_stress_fene) / 3.0
        np.testing.assert_allclose(
            sim_pressure_fene, anal_pressure_fene, atol=tol,
            err_msg='bonded pressure for fene does not match analytical result')

        # Compare stress tensor observable to stress tensor from analysis
        np.testing.assert_allclose(
            StressTensor().calculate(),
            system.analysis.stress_tensor()["total"],
            rtol=0, atol=1E-10)

        system.part.clear()


if __name__ == "__main__":
    ut.main()
