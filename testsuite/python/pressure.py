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
import espressomd.interactions
import espressomd.observables

import tests_common

import numpy as np

# allowed deviation from analytical results
tol = 1.0e-13


def pressure_tensor_kinetic(vel):
    '''Analytical result for convective pressure tensor'''
    return np.einsum('ij,ik->jk', vel, vel) / system.volume()


def pressure_tensor_bonded(pos):
    '''Analytical result for pressure tensor originating from bonded forces'''
    tensor = np.zeros([3, 3])
    for p1, p2 in zip(pos[0::2], pos[1::2]):
        r = p1 - p2
        f = -1.0e4 * r
        tensor += np.einsum('i,j', f, r) / system.volume()
    return tensor


def pressure_tensor_nonbonded(particle_pairs):
    '''Analytical result for pressure tensor originating from non-bonded forces'''
    tensor = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if (p1.type == 0 and p2.type == 0) or (p1.type == 1 and p2.type == 2):
            d = p1.pos - p2.pos
            r = np.linalg.norm(d)
            r_hat = d / r
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / r**13 - 1.0**6 / r**7)) * r_hat
            tensor += np.einsum('i,j', f, d) / system.volume()
    return tensor


def pressure_tensor_nonbonded_inter(particle_pairs):
    tensor = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if p1.type == 1 and p2.type == 2 and p1.mol_id != p2.mol_id:
            r = p1.pos - p2.pos
            d = np.linalg.norm(r)
            r_hat = r / d
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / d**13 - 1.0**6 / d**7)) * r_hat
            tensor += np.einsum('i,j', f, r) / system.volume()
    return tensor


def pressure_tensor_nonbonded_intra(particle_pairs):
    tensor = np.zeros([3, 3])
    for p1, p2 in particle_pairs:
        if p1.type == 0 and p2.type == 0 and p1.mol_id == p2.mol_id:
            r = p1.pos - p2.pos
            d = np.linalg.norm(r)
            r_hat = r / d
            f = (24.0 * 1.0 * (2.0 * 1.0**12 / d**13 - 1.0**6 / d**7)) * r_hat
            tensor += np.einsum('i,j', f, r) / system.volume()
    return tensor


system = espressomd.System(box_l=[1.0, 1.0, 1.0])


@utx.skipIfMissingFeatures(['LENNARD_JONES'])
class PressureLJ(ut.TestCase):

    def tearDown(self):
        system.part.clear()

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
        p0 = system.part.add(pos=[9.9, 9.75, 9.9], type=0, mol_id=0)
        p1 = system.part.add(pos=[9.9, 10.25, 9.9], type=0, mol_id=0)
        p2 = system.part.add(pos=[0.1, 9.7, 0.1], type=1, mol_id=1)
        p3 = system.part.add(pos=[0.1, 10.3, 0.1], type=2, mol_id=2)

        harmonic = espressomd.interactions.HarmonicBond(k=1e4, r_0=0)
        system.bonded_inter.add(harmonic)
        p0.add_bond((harmonic, p1))
        p2.add_bond((harmonic, p3))

        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.non_bonded_inter[1, 2].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)

        system.integrator.run(steps=0)

        p0.v = [10.0, 20.0, 30.0]
        p1.v = [-15., -25., -35.]
        p2.v = [27.0, 23.0, 17.0]
        p3.v = [13.0, 11.0, 19.0]

        partcls = system.part.all()
        pos = partcls.pos
        vel = partcls.v

        sim_pressure_tensor = system.analysis.pressure_tensor()
        sim_pressure_tensor_kinetic = np.copy(sim_pressure_tensor['kinetic'])
        sim_pressure_tensor_bonded = np.copy(sim_pressure_tensor['bonded'])
        sim_pressure_tensor_bonded_harmonic = np.copy(
            sim_pressure_tensor['bonded', len(system.bonded_inter) - 1])
        sim_pressure_tensor_nonbonded = np.copy(
            sim_pressure_tensor['non_bonded'])
        sim_pressure_tensor_nonbonded_inter = np.copy(
            sim_pressure_tensor['non_bonded_inter'])
        sim_pressure_tensor_nonbonded_inter12 = np.copy(
            sim_pressure_tensor['non_bonded_inter', 1, 2])
        sim_pressure_tensor_nonbonded_intra = np.copy(
            sim_pressure_tensor['non_bonded_intra'])
        sim_pressure_tensor_nonbonded_intra00 = np.copy(
            sim_pressure_tensor['non_bonded_intra', 0, 0])
        sim_pressure_tensor_total = np.copy(sim_pressure_tensor['total'])

        sim_pressure = system.analysis.pressure()
        sim_pressure_kinetic = sim_pressure['kinetic']
        sim_pressure_bonded = sim_pressure['bonded']
        sim_pressure_bonded_harmonic = sim_pressure[
            'bonded', len(system.bonded_inter) - 1]
        sim_pressure_nonbonded = sim_pressure['non_bonded']
        sim_pressure_nonbonded_inter = sim_pressure['non_bonded_inter']
        sim_pressure_nonbonded_inter12 = sim_pressure['non_bonded_inter', 1, 2]
        sim_pressure_nonbonded_intra = sim_pressure['non_bonded_intra']
        sim_pressure_nonbonded_intra00 = sim_pressure['non_bonded_intra', 0, 0]
        sim_pressure_total = sim_pressure['total']

        anal_pressure_tensor_kinetic = pressure_tensor_kinetic(vel)
        anal_pressure_tensor_bonded = pressure_tensor_bonded(pos)
        anal_pressure_tensor_nonbonded = pressure_tensor_nonbonded(
            system.part.pairs())
        anal_pressure_tensor_nonbonded_inter = pressure_tensor_nonbonded_inter(
            system.part.pairs())
        anal_pressure_tensor_nonbonded_intra = pressure_tensor_nonbonded_intra(
            system.part.pairs())
        anal_pressure_tensor_total = anal_pressure_tensor_kinetic + \
            anal_pressure_tensor_bonded + anal_pressure_tensor_nonbonded
        anal_pressure_kinetic = np.einsum(
            'ii', anal_pressure_tensor_kinetic) / 3.0
        anal_pressure_bonded = np.einsum(
            'ii', anal_pressure_tensor_bonded) / 3.0
        anal_pressure_nonbonded = np.einsum(
            'ii', anal_pressure_tensor_nonbonded) / 3.0
        anal_pressure_nonbonded_inter = np.einsum(
            'ii', anal_pressure_tensor_nonbonded_inter) / 3.0
        anal_pressure_nonbonded_intra = np.einsum(
            'ii', anal_pressure_tensor_nonbonded_intra) / 3.0
        anal_pressure_total = anal_pressure_kinetic + \
            anal_pressure_bonded + anal_pressure_nonbonded

        np.testing.assert_allclose(
            sim_pressure_tensor_kinetic, anal_pressure_tensor_kinetic, rtol=0, atol=tol,
            err_msg='kinetic pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_bonded, anal_pressure_tensor_bonded, rtol=0, atol=tol,
            err_msg='bonded pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_bonded_harmonic, anal_pressure_tensor_bonded, rtol=0, atol=tol,
            err_msg='bonded pressure tensor harmonic bond does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_nonbonded, anal_pressure_tensor_nonbonded, rtol=0, atol=tol,
            err_msg='non-bonded pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_nonbonded_inter, anal_pressure_tensor_nonbonded_inter, rtol=0, atol=tol,
            err_msg='non-bonded intermolecular pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_nonbonded_inter12, anal_pressure_tensor_nonbonded_inter, rtol=0, atol=tol,
            err_msg='non-bonded intermolecular pressure tensor molecules 1 and 2 does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_nonbonded_intra, anal_pressure_tensor_nonbonded_intra, rtol=0, atol=tol,
            err_msg='non-bonded intramolecular pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_nonbonded_intra00, anal_pressure_tensor_nonbonded_intra, rtol=0, atol=tol,
            err_msg='non-bonded intramolecular pressure tensor molecule 0 does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_total, anal_pressure_tensor_total, rtol=0, atol=tol,
            err_msg='total pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_total, sim_pressure_tensor_kinetic + sim_pressure_tensor_bonded + sim_pressure_tensor_nonbonded, rtol=0, atol=tol,
            err_msg='total pressure tensor is not given as the sum of all major pressure components')
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

        # Compare pressure observables to pressure from analysis
        np.testing.assert_allclose(
            espressomd.observables.PressureTensor().calculate(),
            sim_pressure_tensor["total"],
            rtol=0, atol=1E-10)
        self.assertAlmostEqual(
            espressomd.observables.Pressure().calculate(),
            sim_pressure["total"],
            delta=tol)


@utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
class PressureFENE(ut.TestCase):

    def tearDown(self):
        system.part.clear()

    def get_anal_pressure_tensor_fene(self, pos_1, pos_2, k, d_r_max, r_0):
        tensor = np.zeros([3, 3])
        vec_r = pos_1 - pos_2
        f = -tests_common.fene_force2(vec_r, k, d_r_max, r_0)
        tensor += np.einsum('i,j', f, vec_r) / system.volume()
        return tensor

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
        p0 = system.part.add(
            pos=[9.9, 9.75, 9.9], type=0, mol_id=0, fix=[1, 1, 1])
        p1 = system.part.add(
            pos=[9.9, 10.25, 9.9], type=0, mol_id=0, fix=[1, 1, 1])

        k = 1e4
        d_r_max = 1.5
        r_0 = 0.1

        fene = espressomd.interactions.FeneBond(k=k, d_r_max=d_r_max, r_0=r_0)
        system.bonded_inter.add(fene)
        p0.add_bond((fene, p1))
        system.integrator.run(steps=0)

        sim_pressure_tensor = system.analysis.pressure_tensor()
        sim_pressure_tensor_bonded = sim_pressure_tensor['bonded']
        sim_pressure_tensor_fene = sim_pressure_tensor['bonded', len(
            system.bonded_inter) - 1]

        total_bonded_pressure_tensor = np.zeros([3, 3])
        for i in range(len(system.bonded_inter)):
            total_bonded_pressure_tensor += sim_pressure_tensor['bonded', i]

        anal_pressure_tensor_fene = self.get_anal_pressure_tensor_fene(
            p0.pos, p1.pos, k, d_r_max, r_0)
        np.testing.assert_allclose(
            sim_pressure_tensor_bonded, anal_pressure_tensor_fene, atol=tol,
            err_msg='bonded pressure tensor does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_fene, anal_pressure_tensor_fene, atol=tol,
            err_msg='bonded pressure tensor for fene does not match analytical result')
        np.testing.assert_allclose(
            sim_pressure_tensor_bonded, total_bonded_pressure_tensor, atol=tol,
            err_msg='bonded pressure tensor do not sum up to the total value')

        sim_pressure = system.analysis.pressure()
        sim_pressure_fene = sim_pressure['bonded', len(
            system.bonded_inter) - 1]
        anal_pressure_fene = np.einsum("ii", anal_pressure_tensor_fene) / 3.0
        np.testing.assert_allclose(
            sim_pressure_fene, anal_pressure_fene, atol=tol,
            err_msg='bonded pressure for fene does not match analytical result')

        # Compare pressure observables to pressure from analysis
        np.testing.assert_allclose(
            espressomd.observables.PressureTensor().calculate(),
            sim_pressure_tensor["total"],
            rtol=0, atol=1E-10)
        self.assertAlmostEqual(
            espressomd.observables.Pressure().calculate(),
            sim_pressure["total"],
            delta=tol)


if __name__ == "__main__":
    ut.main()
