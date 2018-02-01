#!/usr/bin/env python

from __future__ import print_function

import unittest as ut
import espressomd
from espressomd.interactions import HarmonicBond
from espressomd import analyze
import itertools

import numpy as np

# allowed deviation from analytical results
tol = 1.0e-13

# analytical result for convective stress
def stress_kinetic(vel, box_l):
  return np.einsum('ij,ik->jk', vel, vel) / box_l**3

# analytical result for stress originating from bond force
def stress_bonded(pos, box_l):
  stress = np.zeros([3,3])
  for p1, p2 in zip(pos[0::2], pos[1::2]):
    r = p1 - p2
    f = -1.0e4*r
    stress += np.einsum('i,j', f, r) / box_l**3
  return stress

# analytical result for stress originating from non-bonded force
def stress_nonbonded(particle_pairs, box_l):
  stress = np.zeros([3,3])
  for p1, p2 in particle_pairs:
    if (p1.type == 0 and p2.type == 0) or (p1.type == 1 and p2.type == 2):
      d = p1.pos - p2.pos
      r = np.sqrt(np.sum(d**2))
      r_hat = d/r
      f = (24.0 * 1.0 * (2.0*1.0**12/r**13 - 1.0**6/r**7)) * r_hat
      stress += np.einsum('i,j', f, d) / box_l**3
  return stress

def stress_nonbonded_inter(particle_pairs, box_l):
  stress = np.zeros([3,3])
  for p1, p2 in particle_pairs:
    if p1.type == 1 and p2.type == 2 and p1.mol_id != p2.mol_id:
      r = p1.pos - p2.pos
      d = np.sqrt(np.sum(r**2))
      r_hat = r/d
      f = (24.0 * 1.0 * (2.0*1.0**12/d**13 - 1.0**6/d**7)) * r_hat
      stress += np.einsum('i,j', f, r) / box_l**3
  return stress

def stress_nonbonded_intra(particle_pairs, box_l):
  stress = np.zeros([3,3])
  for p1, p2 in particle_pairs:
    if p1.type == 0 and p2.type == 0 and p1.mol_id == p2.mol_id:
      r = p1.pos - p2.pos
      d = np.sqrt(np.sum(r**2))
      r_hat = r/d
      f = (24.0 * 1.0 * (2.0*1.0**12/d**13 - 1.0**6/d**7)) * r_hat
      stress += np.einsum('i,j', f, r) / box_l**3
  return stress

system = espressomd.System()

@ut.skipIf(not espressomd.has_features(['LENNARD_JONES']),
'Features not available, skipping test!')
class stress_test(ut.TestCase):

  def test(self):
    # system parameters
    box_l        = 10.0
    system.box_l = [box_l, box_l, box_l]
    skin        = 0.4
    time_step   = 0.01
    system.time_step = time_step

    # thermostat and cell system
    system.thermostat.set_langevin(kT=0.0, gamma=1.0)
    system.cell_system.skin = skin
    system.periodicity = [1,1,1]

    # particles and bond
    system.part.add(id=0, pos=[9.9,9.75,9.9], type=0, mol_id=0)
    system.part.add(id=1, pos=[9.9,10.25,9.9], type=0, mol_id=0)
    system.part.add(id=2, pos=[0.1,9.7,0.1], type=1, mol_id=1)
    system.part.add(id=3, pos=[0.1,10.3,0.1], type=2, mol_id=2)

    harmonic=HarmonicBond(k=1e4, r_0=0)
    system.bonded_inter.add(harmonic)
    system.part[0].add_bond((harmonic, 1))
    system.part[2].add_bond((harmonic, 3))

    system.non_bonded_inter[0, 0].lennard_jones.set_params(
               epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
    system.non_bonded_inter[1, 2].lennard_jones.set_params(
                epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)

    # calculate forces and stress without Lees-Edwards offset
    system.integrator.run(steps=0)

    system.part[0].v = [10.0,20.0,30.0]
    system.part[1].v = [-15,-25,-35]
    system.part[2].v = [27.0,23.0,17.0]
    system.part[3].v = [13.0,11.0,19.0]

    pos = system.part[:].pos
    vel = system.part[:].v

    sim_stress_kinetic = system.analysis.stress_tensor()['kinetic']
    sim_stress_bonded = system.analysis.stress_tensor()['bonded']
    sim_stress_bonded_harmonic = system.analysis.stress_tensor()['bonded', len(system.bonded_inter)-1]
    sim_stress_nonbonded = system.analysis.stress_tensor()['non_bonded']
    sim_stress_nonbonded_inter = system.analysis.stress_tensor()['non_bonded_inter']
    sim_stress_nonbonded_inter12 = system.analysis.stress_tensor()['non_bonded_inter', 1, 2]
    sim_stress_nonbonded_intra = system.analysis.stress_tensor()['non_bonded_intra']
    sim_stress_nonbonded_intra00 = system.analysis.stress_tensor()['non_bonded_intra', 0, 0]
    sim_stress_total = system.analysis.stress_tensor()['total']
    sim_pressure_kinetic = system.analysis.pressure()['kinetic']
    sim_pressure_bonded = system.analysis.pressure()['bonded']
    sim_pressure_bonded_harmonic = system.analysis.pressure()['bonded', len(system.bonded_inter)-1]
    sim_pressure_nonbonded = system.analysis.pressure()['non_bonded']
    sim_pressure_nonbonded_inter = system.analysis.pressure()['non_bonded_inter']
    sim_pressure_nonbonded_inter12 = system.analysis.pressure()['non_bonded_inter', 1, 2]
    sim_pressure_nonbonded_intra = system.analysis.pressure()['non_bonded_intra']
    sim_pressure_nonbonded_intra00 = system.analysis.pressure()['non_bonded_intra', 0, 0]
    sim_pressure_total = system.analysis.pressure()['total']

    anal_stress_kinetic = stress_kinetic(vel, box_l)
    anal_stress_bonded = stress_bonded(pos, box_l)
    anal_stress_nonbonded = stress_nonbonded(system.part.pairs(), box_l)
    anal_stress_nonbonded_inter = stress_nonbonded_inter(system.part.pairs(), box_l)
    anal_stress_nonbonded_intra = stress_nonbonded_intra(system.part.pairs(), box_l)
    anal_stress_total = anal_stress_kinetic + anal_stress_bonded + anal_stress_nonbonded
    anal_pressure_kinetic = np.einsum('ii', anal_stress_kinetic)/3.0
    anal_pressure_bonded = np.einsum('ii', anal_stress_bonded)/3.0
    anal_pressure_nonbonded = np.einsum('ii', anal_stress_nonbonded)/3.0
    anal_pressure_nonbonded_inter = np.einsum('ii', anal_stress_nonbonded_inter)/3.0
    anal_pressure_nonbonded_intra = np.einsum('ii', anal_stress_nonbonded_intra)/3.0
    anal_pressure_total = anal_pressure_kinetic + anal_pressure_bonded + anal_pressure_nonbonded

    print('particle positions')
    print(pos)
    print('particle velocities')
    print(vel)

    print('\nsimulated kinetic pressure=')
    print(sim_pressure_kinetic)
    print('analytical kinetic pressure=')
    print(anal_pressure_kinetic)

    print('\nsimulated bonded pressure=')
    print(sim_pressure_bonded)
    print('analytical bonded pressure=')
    print(anal_pressure_bonded)

    print('\nsimulated bonded pressure harmonic bond =')
    print(sim_pressure_bonded_harmonic)
    print('analytical bonded pressure bond harmonic=')
    print(anal_pressure_bonded)

    print('\nsimulated non-bonded pressure=')
    print(sim_pressure_nonbonded)
    print('analytic non-bonded pressure=')
    print(anal_pressure_nonbonded)

    print('\nsimulated non-bonded intermolecular pressure=')
    print(sim_pressure_nonbonded_inter)
    print('analytic non-bonded intermolecular pressure=')
    print(anal_pressure_nonbonded_inter)

    print('\nsimulated non-bonded intermolecular pressure molecule 1 and 2=')
    print(sim_pressure_nonbonded_inter12)
    print('analytic non-bonded intermolecular pressure molecule 1 and 2=')
    print(anal_pressure_nonbonded_inter)

    print('\nsimulated non-bonded intramolecular pressure=')
    print(sim_pressure_nonbonded_intra)
    print('analytic non-bonded intramolecular pressure=')
    print(anal_pressure_nonbonded_intra)

    print('\nsimulated non-bonded intramolecular pressure molecule 0=')
    print(sim_pressure_nonbonded_intra00)
    print('analytic non-bonded intramolecular pressure molecule 0=')
    print(anal_pressure_nonbonded_intra)

    print('\nsimulated total pressure=')
    print(sim_pressure_total)
    print('analytic total pressure=')
    print(anal_pressure_total)

    print('\nsimulated kinetic stress=')
    print(sim_stress_kinetic)
    print('analytic kinetic stress=')
    print(anal_stress_kinetic)

    print('\nsimulated bonded stress=')
    print(sim_stress_bonded)
    print('analytic bonded stress=')
    print(anal_stress_bonded)

    print('\nsimulated bonded stress harmonic bond=')
    print(sim_stress_bonded_harmonic)
    print('analytic bonded stress harmonic bond=')
    print(anal_stress_bonded)

    print('\nsimulated non-bonded stress=')
    print(sim_stress_nonbonded)
    print('analytic non-bonded stress=')
    print(anal_stress_nonbonded)

    print('\nsimulated non-bonded intermolecular stress=')
    print(sim_stress_nonbonded_inter)
    print('analytic non-bonded intermolecular stress=')
    print(anal_stress_nonbonded_inter)

    print('\nsimulated non-bonded intermolecular stress between molecule 1 and 2=')
    print(sim_stress_nonbonded_inter12)
    print('analytic non-bonded intermolecular stress between molecule 1 and 2=')
    print(anal_stress_nonbonded_inter)

    print('\nsimulated total non-bonded intramolecular stress=')
    print(sim_stress_nonbonded_intra)
    print('analytic total non-bonded intramolecular stress=')
    print(anal_stress_nonbonded_intra)

    print('\nsimulated non-bonded intramolecular stress molecule 0=')
    print(sim_stress_nonbonded_intra00)
    print('analytic non-bonded intramolecular stress molecule 0=')
    print(anal_stress_nonbonded_intra)

    print('\nsimulated total stress=')
    print(sim_stress_total)
    print('analytic total stress=')
    print(anal_stress_total)

    system.part.clear()

    self.assertTrue(np.max(np.abs(sim_stress_kinetic-anal_stress_kinetic)) < tol, 'kinetic stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded-anal_stress_bonded)) < tol, 'bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded_harmonic-anal_stress_bonded)) < tol, 'bonded stress harmonic bond does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded-anal_stress_nonbonded)) < tol, 'non-bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_inter-anal_stress_nonbonded_inter)) < tol, 'non-bonded intermolecular stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_inter12-anal_stress_nonbonded_inter)) < tol, 'non-bonded intermolecular stress molecules 1 and 2 does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_intra-anal_stress_nonbonded_intra)) < tol, 'non-bonded intramolecular stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_intra00-anal_stress_nonbonded_intra)) < tol, 'non-bonded intramolecular stress molecule 0 does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_total-anal_stress_total)) < tol, 'total stress does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_kinetic-anal_pressure_kinetic) < tol, 'kinetic pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded-anal_pressure_bonded) < tol, 'bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded_harmonic-anal_pressure_bonded) < tol, 'bonded pressure harmonic bond does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded-anal_pressure_nonbonded) < tol, 'non-bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_inter-anal_pressure_nonbonded_inter) < tol, 'non-bonded intermolecular pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_inter12-anal_pressure_nonbonded_inter) < tol, 'non-bonded intermolecular pressure molecule 1 and 2 does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_intra-anal_pressure_nonbonded_intra) < tol, 'non-bonded intramolecular pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_intra00-anal_pressure_nonbonded_intra) < tol, 'non-bonded intramolecular pressure molecule 0 does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_total-anal_pressure_total) < tol, 'total pressure does not match analytical result')

@ut.skipIf(not espressomd.has_features(['LEES_EDWARDS', 'LENNARD_JONES']),
'Features not available, skipping test!')
class stress_lees_edwards_test(ut.TestCase):

  def test(self):
    # system parameters
    box_l        = 10.0
    system.box_l = [box_l, box_l, box_l]
    skin        = 0.4
    time_step   = 0.01
    system.time_step = time_step

    # thermostat and cell system
    system.thermostat.set_langevin(kT=0.0, gamma=1.0)
    system.cell_system.skin = skin
    system.periodicity = [1,1,1]

    # particles and bond
    system.part.add(id=0, pos=[9.9,9.75,9.9], type=0, mol_id=0)
    system.part.add(id=1, pos=[9.9,10.25,9.9], type=0, mol_id=0)
    system.part.add(id=2, pos=[0.1,9.7,0.1], type=1, mol_id=1)
    system.part.add(id=3, pos=[0.1,10.3,0.1], type=2, mol_id=2)

    harmonic=HarmonicBond(k=1e4, r_0=0)
    system.bonded_inter.add(harmonic)
    system.part[0].add_bond((harmonic, 1))
    system.part[2].add_bond((harmonic, 3))

    system.non_bonded_inter[0, 0].lennard_jones.set_params(
               epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
    system.non_bonded_inter[1, 2].lennard_jones.set_params(
                epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)

    # calculate forces and stress with Lees-Edwards offset
    le_offset = 0.3
    system.lees_edwards_offset = le_offset
    system.integrator.run(steps=0)

    system.part[0].v = [10.0,20.0,30.0]
    system.part[1].v = [-15,-25,-35]
    system.part[2].v = [27.0,23.0,17.0]
    system.part[3].v = [13.0,11.0,19.0]

    pos = system.part[:].pos
    vel = system.part[:].v

    sim_stress_kinetic = system.analysis.stress_tensor()['kinetic']
    sim_stress_bonded = system.analysis.stress_tensor()['bonded']
    sim_stress_bonded_harmonic = system.analysis.stress_tensor()['bonded', len(system.bonded_inter)-1]
    sim_stress_nonbonded = system.analysis.stress_tensor()['non_bonded']
    sim_stress_nonbonded_inter = system.analysis.stress_tensor()['non_bonded_inter']
    sim_stress_nonbonded_inter12 = system.analysis.stress_tensor()['non_bonded_inter', 1, 2]
    sim_stress_nonbonded_intra = system.analysis.stress_tensor()['non_bonded_intra']
    sim_stress_nonbonded_intra00 = system.analysis.stress_tensor()['non_bonded_intra', 0, 0]
    sim_stress_total = system.analysis.stress_tensor()['total']
    sim_pressure_kinetic = system.analysis.pressure()['kinetic']
    sim_pressure_bonded = system.analysis.pressure()['bonded']
    sim_pressure_bonded_harmonic = system.analysis.pressure()['bonded', len(system.bonded_inter)-1]
    sim_pressure_nonbonded = system.analysis.pressure()['non_bonded']
    sim_pressure_nonbonded_inter = system.analysis.pressure()['non_bonded_inter']
    sim_pressure_nonbonded_inter12 = system.analysis.pressure()['non_bonded_inter', 1, 2]
    sim_pressure_nonbonded_intra = system.analysis.pressure()['non_bonded_intra']
    sim_pressure_nonbonded_intra00 = system.analysis.pressure()['non_bonded_intra', 0, 0]
    sim_pressure_total = system.analysis.pressure()['total']

    anal_stress_kinetic = stress_kinetic(vel, box_l)
    anal_stress_bonded = stress_bonded(pos, box_l)
    anal_stress_nonbonded = stress_nonbonded(system.part.pairs(), box_l)
    anal_stress_nonbonded_inter = stress_nonbonded_inter(system.part.pairs(), box_l)
    anal_stress_nonbonded_intra = stress_nonbonded_intra(system.part.pairs(), box_l)
    anal_stress_total = anal_stress_kinetic + anal_stress_bonded + anal_stress_nonbonded
    anal_pressure_kinetic = np.einsum('ii', anal_stress_kinetic)/3.0
    anal_pressure_bonded = np.einsum('ii', anal_stress_bonded)/3.0
    anal_pressure_nonbonded = np.einsum('ii', anal_stress_nonbonded)/3.0
    anal_pressure_nonbonded_inter = np.einsum('ii', anal_stress_nonbonded_inter)/3.0
    anal_pressure_nonbonded_intra = np.einsum('ii', anal_stress_nonbonded_intra)/3.0
    anal_pressure_total = anal_pressure_kinetic + anal_pressure_bonded + anal_pressure_nonbonded

    print('le_offset = {}'.format(system.lees_edwards_offset))

    print('particle positions')
    print(pos)
    print('particle velocities')
    print(vel)

    print('\nsimulated kinetic pressure=')
    print(sim_pressure_kinetic)
    print('analytical kinetic pressure=')
    print(anal_pressure_kinetic)

    print('\nsimulated bonded pressure=')
    print(sim_pressure_bonded)
    print('analytical bonded pressure=')
    print(anal_pressure_bonded)

    print('\nsimulated bonded pressure harmonic bond=')
    print(sim_pressure_bonded_harmonic)
    print('analytical bonded pressure harmonic bond=')
    print(anal_pressure_bonded)

    print('\nsimulated non-bonded pressure=')
    print(sim_pressure_nonbonded)
    print('analytic non-bonded pressure=')
    print(anal_pressure_nonbonded)

    print('\nsimulated non-bonded intermolecular pressure=')
    print(sim_pressure_nonbonded_inter)
    print('analytic non-bonded intermolecular pressure=')
    print(anal_pressure_nonbonded_inter)

    print('\nsimulated non-bonded intermolecular pressure molecule 1 and 2=')
    print(sim_pressure_nonbonded_inter12)
    print('analytic non-bonded intermolecular pressure molecule 1 and 2=')
    print(anal_pressure_nonbonded_inter)

    print('\nsimulated non-bonded intramolecular pressure=')
    print(sim_pressure_nonbonded_intra)
    print('analytic non-bonded intramolecular pressure=')
    print(anal_pressure_nonbonded_intra)

    print('\nsimulated non-bonded intramolecular pressure molecule 0=')
    print(sim_pressure_nonbonded_intra00)
    print('analytic non-bonded intramolecular pressure molecule 0=')
    print(anal_pressure_nonbonded_intra)

    print('\nsimulated total pressure=')
    print(sim_pressure_total)
    print('analytic total pressure=')
    print(anal_pressure_total)

    print('\nsimulated kinetic stress=')
    print(sim_stress_kinetic)
    print('analytic kinetic stress=')
    print(anal_stress_kinetic)

    print('\nsimulated bonded stress=')
    print(sim_stress_bonded)
    print('analytic bonded stress=')
    print(anal_stress_bonded)

    print('\nsimulated bonded stress harmonic bond=')
    print(sim_stress_bonded_harmonic)
    print('analytic bonded stress harmonic bond=')
    print(anal_stress_bonded)

    print('\nsimulated non-bonded stress=')
    print(sim_stress_nonbonded)
    print('analytic non-bonded stress=')
    print(anal_stress_nonbonded)

    print('\nsimulated non-bonded intermolecular stress=')
    print(sim_stress_nonbonded_inter)
    print('analytic non-bonded intermolecular stress=')
    print(anal_stress_nonbonded_inter)

    print('\nsimulated non-bonded intermolecular stress between molecule 1 and 2=')
    print(sim_stress_nonbonded_inter12)
    print('analytic non-bonded intermolecular stress between molecule 1 and 2=')
    print(anal_stress_nonbonded_inter)

    print('\nsimulated total non-bonded intramolecular stress=')
    print(sim_stress_nonbonded_intra)
    print('analytic total non-bonded intramolecular stress=')
    print(anal_stress_nonbonded_intra)

    print('\nsimulated non-bonded intramolecular stress molecule 0=')
    print(sim_stress_nonbonded_intra00)
    print('analytic non-bonded intramolecular stress molecule 0=')
    print(anal_stress_nonbonded_intra)

    print('\nsimulated total stress=')
    print(sim_stress_total)
    print('analytic total stress=')
    print(anal_stress_total)

    system.part.clear()

    self.assertTrue(np.max(np.abs(sim_stress_kinetic-anal_stress_kinetic)) < tol, 'kinetic stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded-anal_stress_bonded)) < tol, 'bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded_harmonic-anal_stress_bonded)) < tol, 'bonded stress harmonic bond does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded-anal_stress_nonbonded)) < tol, 'non-bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_inter-anal_stress_nonbonded_inter)) < tol, 'non-bonded intermolecular stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_inter12-anal_stress_nonbonded_inter)) < tol, 'non-bonded intermolecular stress molecules 1 and 2 does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_intra-anal_stress_nonbonded_intra)) < tol, 'non-bonded intramolecular stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded_intra00-anal_stress_nonbonded_intra)) < tol, 'non-bonded intramolecular stress molecule 0 does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_total-anal_stress_total)) < tol, 'total stress does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_kinetic-anal_pressure_kinetic) < tol, 'kinetic pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded-anal_pressure_bonded) < tol, 'bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded_harmonic-anal_pressure_bonded) < tol, 'bonded pressure harmonic bond does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded-anal_pressure_nonbonded) < tol, 'non-bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_inter-anal_pressure_nonbonded_inter) < tol, 'non-bonded intermolecular pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_inter12-anal_pressure_nonbonded_inter) < tol, 'non-bonded intermolecular pressure molecule 1 and 2 does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_intra-anal_pressure_nonbonded_intra) < tol, 'non-bonded intramolecular pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded_intra00-anal_pressure_nonbonded_intra) < tol, 'non-bonded intramolecular pressure molecule 0 does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_total-anal_pressure_total) < tol, 'total pressure does not match analytical result')

if __name__ == "__main__":
  ut.main()
