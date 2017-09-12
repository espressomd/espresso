#!/usr/bin/env python

from __future__ import print_function

import unittest as ut
import espressomd
from espressomd.interactions import HarmonicBond
from espressomd import analyze

import numpy as np
import scipy as sp

# allowed deviation from analytical results
tol = 1e-14

# analytical result for convective stress
def stress_ideal(vel, box_l):
  return np.einsum('ij,ik->jk', vel, vel) / box_l**3

# analytical result for stress originating from bond force
def stress_bonded(pos, box_l):
  r = pos[1] - pos[0]
  f = -1.0*r
  return np.einsum('i,j', f, r) / box_l**3

# analytical result for stress originating from non-bonded force
def stress_nonbonded(pos, box_l):
  d = pos[1] - pos[0]
  r = np.sqrt(np.sum(d**2))
  r_hat = d/r
  f = (24.0 * 1.0 * (2.0*1.0**12/r**13 - 1.0**6/r**7)) * r_hat
  return np.einsum('i,j', f, d) / box_l**3

@ut.skipIf(not espressomd.has_features(['LEES_EDWARDS', 'LENNARD_JONES']),
'Features not available, skipping test!')
class lees_edwards_test(ut.TestCase):

  def test(self):
    # system parameters
    system = espressomd.System()
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
    system.part.add(id=0, pos=[9.5,9.75,9.5], type=0)
    system.part.add(id=1, pos=[9.5,10.25,9.5], type=0)

    harmonic=HarmonicBond(k=1, r_0=0)
    system.bonded_inter.add(harmonic)
    system.part[0].add_bond((harmonic, 1))

    system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)

    # calculate forces and stress without Lees-Edwards offset
    le_offset = 0.0
    system.lees_edwards_offset = le_offset
    system.integrator.run(steps=0)

    system.part[0].v = [1.0,2.0,3.0]
    system.part[1].v = [-0.1,-0.2,-0.3]

    pos = system.part[:].pos
    vel = system.part[:].v

    sim_stress_ideal = system.analysis.stress_tensor()['ideal']
    sim_stress_bonded = system.analysis.stress_tensor()['bonded']
    sim_stress_nonbonded = system.analysis.stress_tensor()['non_bonded']
    sim_stress_total = system.analysis.stress_tensor()['total']
    sim_pressure_ideal = system.analysis.pressure()['ideal']
    sim_pressure_bonded = system.analysis.pressure()['bonded']
    sim_pressure_nonbonded = system.analysis.pressure()['non_bonded']
    sim_pressure_total = system.analysis.pressure()['total']

    anal_stress_ideal = stress_ideal(vel, box_l)
    anal_stress_bonded = stress_bonded(pos, box_l)
    anal_stress_nonbonded = stress_nonbonded(pos, box_l)
    anal_stress_total = anal_stress_ideal + anal_stress_bonded + anal_stress_nonbonded
    anal_pressure_ideal = np.einsum('ii', anal_stress_ideal)/3.0
    anal_pressure_bonded = np.einsum('ii', anal_stress_bonded)/3.0
    anal_pressure_nonbonded = np.einsum('ii', anal_stress_nonbonded)/3.0
    anal_pressure_total = anal_pressure_ideal + anal_pressure_bonded + anal_pressure_nonbonded

    print('le_offset = {}'.format(system.lees_edwards_offset))

    print('\nparticle positions')
    print(pos)
    print('particle velocities')
    print(vel)

    print('\nsimulated ideal pressure=')
    print(sim_pressure_ideal)
    print('analytical ideal pressure=')
    print(anal_pressure_ideal)

    print('\nsimulated bonded pressure=')
    print(sim_pressure_bonded)
    print('analytical bonded pressure=')
    print(anal_pressure_bonded)

    print('\nsimulated non-bonded pressure=')
    print(sim_pressure_nonbonded)
    print('analytic non-bonded pressure=')
    print(anal_pressure_nonbonded)

    print('\nsimulated total pressure=')
    print(sim_pressure_total)
    print('analytic total pressure=')
    print(anal_pressure_total)

    print('\nsimulated ideal stress=')
    print(sim_stress_ideal)
    print('analytic ideal stress=')
    print(anal_stress_ideal)

    print('\nsimulated bonded stress=')
    print(sim_stress_bonded)
    print('analytic bonded stress=')
    print(anal_stress_bonded)

    print('\nsimulated non-bonded stress=')
    print(sim_stress_nonbonded)
    print('analytic non-bonded stress=')
    print(anal_stress_nonbonded)

    print('\nsimulated total stress=')
    print(sim_stress_total)
    print('analytic total stress=')
    print(anal_stress_total)

    self.assertTrue(np.max(np.abs(sim_stress_ideal-anal_stress_ideal)) < tol, 'ideal stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded-anal_stress_bonded)) < tol, 'bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded-anal_stress_nonbonded)) < tol, 'non-bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_total-anal_stress_total)) < tol, 'total stress does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_ideal-anal_pressure_ideal) < tol, 'ideal pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded-anal_pressure_bonded) < tol, 'bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded-anal_pressure_nonbonded) < tol, 'non-bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_total-anal_pressure_total) < tol, 'total pressure does not match analytical result')

    # calculate forces and stress with Lees-Edwards offset
    le_offset = 0.3
    system.lees_edwards_offset = le_offset
    system.integrator.run(steps=0)

    system.part[0].v = [1.0,2.0,3.0]
    system.part[1].v = [-0.1,-0.2,-0.3]

    pos = system.part[:].pos
    vel = system.part[:].v

    sim_stress_ideal = system.analysis.stress_tensor()['ideal']
    sim_stress_bonded = system.analysis.stress_tensor()['bonded']
    sim_stress_nonbonded = system.analysis.stress_tensor()['non_bonded']
    sim_stress_total = system.analysis.stress_tensor()['total']
    sim_pressure_ideal = system.analysis.pressure()['ideal']
    sim_pressure_bonded = system.analysis.pressure()['bonded']
    sim_pressure_nonbonded = system.analysis.pressure()['non_bonded']
    sim_pressure_total = system.analysis.pressure()['total']

    anal_stress_ideal = stress_ideal(vel, box_l)
    anal_stress_bonded = stress_bonded(pos, box_l)
    anal_stress_nonbonded = stress_nonbonded(pos, box_l)
    anal_stress_total = anal_stress_ideal + anal_stress_bonded + anal_stress_nonbonded
    anal_pressure_ideal = np.einsum('ii', anal_stress_ideal)/3.0
    anal_pressure_bonded = np.einsum('ii', anal_stress_bonded)/3.0
    anal_pressure_nonbonded = np.einsum('ii', anal_stress_nonbonded)/3.0
    anal_pressure_total = anal_pressure_ideal + anal_pressure_bonded + anal_pressure_nonbonded

    print('le_offset = {}'.format(system.lees_edwards_offset))

    print('\nparticle positions')
    print(pos)
    print('particle velocities')
    print(vel)

    print('\nsimulated ideal pressure=')
    print(sim_pressure_ideal)
    print('analytical ideal pressure=')
    print(anal_pressure_ideal)

    print('\nsimulated bonded pressure=')
    print(sim_pressure_bonded)
    print('analytical bonded pressure=')
    print(anal_pressure_bonded)

    print('\nsimulated non-bonded pressure=')
    print(sim_pressure_nonbonded)
    print('analytic non-bonded pressure=')
    print(anal_pressure_nonbonded)

    print('\nsimulated total pressure=')
    print(sim_pressure_total)
    print('analytic total pressure=')
    print(anal_pressure_total)

    print('\nsimulated ideal stress=')
    print(sim_stress_ideal)
    print('analytic ideal stress=')
    print(anal_stress_ideal)

    print('\nsimulated bonded stress=')
    print(sim_stress_bonded)
    print('analytic bonded stress=')
    print(anal_stress_bonded)

    print('\nsimulated non-bonded stress=')
    print(sim_stress_nonbonded)
    print('analytic non-bonded stress=')
    print(anal_stress_nonbonded)

    print('\nsimulated total stress=')
    print(sim_stress_total)
    print('analytic total stress=')
    print(anal_stress_total)

    self.assertTrue(np.max(np.abs(sim_stress_ideal-anal_stress_ideal)) < tol, 'ideal stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_bonded-anal_stress_bonded)) < tol, 'bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_nonbonded-anal_stress_nonbonded)) < tol, 'non-bonded stress does not match analytical result')
    self.assertTrue(np.max(np.abs(sim_stress_total-anal_stress_total)) < tol, 'total stress does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_ideal-anal_pressure_ideal) < tol, 'ideal pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_bonded-anal_pressure_bonded) < tol, 'bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_nonbonded-anal_pressure_nonbonded) < tol, 'non-bonded pressure does not match analytical result')
    self.assertTrue(np.abs(sim_pressure_total-anal_pressure_total) < tol, 'total pressure does not match analytical result')

if __name__ == "__main__":
  ut.main()
