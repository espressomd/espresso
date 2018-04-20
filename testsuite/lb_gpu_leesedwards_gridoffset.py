#!/usr/bin/env python

from __future__ import print_function

import espressomd as md
from espressomd import lb
import numpy as np
import unittest as ut
from tests_common import abspath

tol = 1.0e-16

# LB parameter for offset and particle coupling test
box_l = 9
eta = 1.0
rho = 1.0
nu = eta / rho
agrid = 1.0
fric = 1.0
lb1 = md.lb.LBFluidGPU(agrid=agrid, dens=rho, visc=eta, fric=fric, tau=time_step)

# LB parameter for test against Navier-Stokes equation
box_l = 9
eta = 0.27
rho = 0.97
nu = eta / rho
v = 0.0087
k_max = 100
agrid = 1.0
fric = 1.1
lb2 = md.lb.LBFluidGPU(agrid=agrid, dens=rho, visc=eta, fric=fric, tau=time_step)

system = md.System(box_l=[1.0, 1.0, 1.0])

# Integration parameters
time_step = 1.0
system.time_step = time_step
system.cell_system.skin = 0.4
total_time = 1000.0

@ut.skipIf(not md.has_features(['LEES_EDWARDS', 'LB_GPU', 'EXTERNAL_FORCES']),
  'Features not available, skipping test!')
class LeesEdwardsOneGridOffsetTest(ut.TestCase):

  def test(self):
    
    """The velocity components at all LB nodes are tested against a table stored in testsuite/data/"""   

    # Add a fixed particle with a velovity in the y-direction
    system.part.add(pos=[4.5, 7.5, 4.5], fix=[1, 1, 1], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.lees_edwards_offset = 1.0
    system.integrator.run(2)

    # Add LB fluid
    system.actors.add(lb1)

    system.integrator.run(2)

    node_list_pos_v = np.empty([1, 6])

    for i in range(int(box_l / agrid)):
      for j in range(int(box_l / agrid)):
        for k in range(int(box_l / agrid)):
          node_pos = [i, j, k]
          node_pos_v = np.concatenate((node_pos, lb[i, j, k].velocity))
          node_list_pos_v = np.vstack((node_list_pos_v, node_pos_v))

    node_list_pos_v = node_list_pos_v[1:, :]
    saved_data = np.loadtxt(abspath("data/lb_gpu_leesedwards_gridoffset.dat"))

    nonzero_column_calc, nonzero_row_calc = np.nonzero(node_list_pos_v[:, 3:])
    nonzero_position_calc = np.transpose(np.vstack((nonzero_column_calc, nonzero_row_calc)))

    nonzero_column_saved, nonzero_row_saved = np.nonzero(saved_data[:, 3:])
    nonzero_position_saved = np.transpose(np.vstack((nonzero_column_saved, nonzero_row_saved)))

    self.assertTrue(np.allclose(nonzero_position_calc, nonzero_position_saved, rtol=tol, atol=tol), 'nodes with non-zero populations do not match')
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do not match')

    system.part.clear()
    system.actors.remove(lb_1)

class LeesEdwardsHalfGridOffsetTest(ut.TestCase):

  def test(self):
    
    """The velocity components at all LB nodes are tested against a table stored in testsuite/data/"""   
    
    # Add a fixed particle with a velovity in the y-direction
    system.part.add(pos=[4.5, 7.5, 4.5], fix=[1, 1, 1], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.lees_edwards_offset = 0.5
    system.integrator.run(2)

    # Add LB fluid
    system.actors.add(lb1)

    system.integrator.run(2)

    node_list_pos_v = np.empty([1, 6])

    for i in range(int(box_l / agrid)):
      for j in range(int(box_l / agrid)):
        for k in range(int(box_l / agrid)):
          node_pos = [i, j, k]
          node_pos_v = np.concatenate((node_pos, lb[i, j, k].velocity))
          node_list_pos_v = np.vstack((node_list_pos_v, node_pos_v))

    node_list_pos_v = node_list_pos_v[1:, :]
    saved_data = np.loadtxt(abspath("data/lb_gpu_leesedwards_halfgridoffset.dat"))

    nonzero_column_calc, nonzero_row_calc = np.nonzero(node_list_pos_v[:, 3:])
    nonzero_position_calc = np.transpose(np.vstack((nonzero_column_calc, nonzero_row_calc)))

    nonzero_column_saved, nonzero_row_saved = np.nonzero(saved_data[:, 3:])
    nonzero_position_saved = np.transpose(np.vstack((nonzero_column_saved, nonzero_row_saved)))

    self.assertTrue(np.allclose(nonzero_position_calc, nonzero_position_saved, rtol=tol, atol=tol), 'nodes with non-zero populations do not match')
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do not match')

    system.part.clear()
    system.actors.remove(lb_1)

class LeesEdwardsOddGridOffsetTest(ut.TestCase):

  def test(self):
    
    """The velocity components at all LB nodes are tested against a table stored in testsuite/data/"""   
    
    # Add a fixed particle with a velovity in the y-direction
    system.part.add(pos=[4.5, 7.5, 4.5], fix=[1, 1, 1], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.lees_edwards_offset = 1.7
    system.integrator.run(2)

    # Add LB fluid
    system.actors.add(lb1)

    system.integrator.run(2)

    node_list_pos_v = np.empty([1, 6])

    for i in range(int(box_l / agrid)):
      for j in range(int(box_l / agrid)):
        for k in range(int(box_l / agrid)):
          node_pos = [i, j, k]
          node_pos_v = np.concatenate((node_pos, lb[i, j, k].velocity))
          node_list_pos_v = np.vstack((node_list_pos_v, node_pos_v))

    node_list_pos_v = node_list_pos_v[1:, :]
    saved_data = np.loadtxt(abspath("data/lb_gpu_leesedwards_oddgridoffset.dat"))

    nonzero_column_calc, nonzero_row_calc = np.nonzero(node_list_pos_v[:, 3:])
    nonzero_position_calc = np.transpose(np.vstack((nonzero_column_calc, nonzero_row_calc)))

    nonzero_column_saved, nonzero_row_saved = np.nonzero(saved_data[:, 3:])
    nonzero_position_saved = np.transpose(np.vstack((nonzero_column_saved, nonzero_row_saved)))

    self.assertTrue(np.allclose(nonzero_position_calc, nonzero_position_saved, rtol=tol, atol=tol), 'nodes with non-zero populations do not match')
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do not match')

    system.part.clear()
    system.actors.remove(lb_1)

class LeesEdwardsNoOffsetTest(ut.TestCase):

  def test(self):

    """The velocity components at all LB nodes are tested against a table stored in testsuite/data/""" 

    # Add a fixed particle with a velovity in the y-direction
    system.part.add(pos=[4.5, 7.5, 4.5], fix=[1, 1, 1], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.lees_edwards_offset = 0.0
    system.integrator.run(2)

    # Add LB fluid
    system.actors.add(lb1)

    system.integrator.run(2)

    node_list_pos_v = np.empty([1, 6])

    for i in range(int(box_l / agrid)):
      for j in range(int(box_l / agrid)):
        for k in range(int(box_l / agrid)):
          node_pos = [i, j, k]
          node_pos_v = np.concatenate((node_pos, lb[i, j, k].velocity))
          node_list_pos_v = np.vstack((node_list_pos_v, node_pos_v))

    node_list_pos_v = node_list_pos_v[1:, :]
    saved_data = np.loadtxt(abspath("data/lb_gpu_leesedwards_nooffset.dat"))

    nonzero_column_calc, nonzero_row_calc = np.nonzero(node_list_pos_v[:, 3:])
    nonzero_position_calc = np.transpose(np.vstack((nonzero_column_calc, nonzero_row_calc)))

    nonzero_column_saved, nonzero_row_saved = np.nonzero(saved_data[:, 3:])
    nonzero_position_saved = np.transpose(np.vstack((nonzero_column_saved, nonzero_row_saved)))

    self.assertTrue(np.allclose(nonzero_position_calc, nonzero_position_saved, rtol=tol, atol=tol), 'nodes with non-zero populations do not match')
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do not match')

    system.part.clear()
    system.actors.remove(lb_1)

class LeesEdwardsParticleCouplingTest(ut.TestCase):

  def test(self):

    """The velocity components at all LB nodes are tested against a table
    stored in testsuite/data/"""

    # Add LB fluid
    system.actors.add(lb1)

    system.lees_edwards_offset = 2.5
    counter = 0

    # Add a particle with a velocity in the y-direction
    system.part.add(pos=[4.5, 8.9, 4.5], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.integrator.run(10)

    node_list_pos_v = np.empty([1, 6])

    for i in range(int(box_l / agrid)):
      for j in range(int(box_l / agrid)):
        for k in range(int(box_l / agrid)):
          node_pos = [i, j, k]
          node_pos_v = np.concatenate((node_pos, lb[i, j, k].velocity))
          node_list_pos_v = np.vstack((node_list_pos_v, node_pos_v))

    node_list_pos_v = node_list_pos_v[1:, :]

    saved_data = np.loadtxt(abspath("data/lb_gpu_leesedwards_particlecoupling.dat"))

    nonzero_column_calc, nonzero_row_calc = np.nonzero(node_list_pos_v[:, 3:])
    nonzero_position_calc = np.transpose(np.vstack((nonzero_column_calc, nonzero_row_calc)))

    nonzero_column_saved, nonzero_row_saved = np.nonzero(saved_data[:, 3:])
    nonzero_position_saved = np.transpose(np.vstack((nonzero_column_saved, nonzero_row_saved)))

    self.assertTrue(np.allclose(nonzero_position_calc, nonzero_position_saved, rtol=tol, atol=tol), 'nodes with non-zero populations do not match')
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do notmatch')

    system.part.clear()
    system.actors.remove(lb_1)

class LBGPULeesEdwardsTest(ut.TestCase):
  def test(self):

    """In this test the velocity profile of a LB fluid under steady shear is compared to the velocity profile that is obtained by solving the Navier-Stokes equation with Fourier series."""

    # Analytical solution with Fourier series of Navier-Stokes equation
    def u(x, t, nu, v, h, k_max):
      u = x / h - 0.5
      for k in np.arange(1, k_max + 1):
        u += 1.0 / (np.pi * k) * np.exp(-4 * np.pi ** 2 * nu * k ** 2 / h ** 2 * t) * np.sin(2 * np.pi / h * k * x)
      return v * u

    # Add LB fluid
    system.actors.add(lb2)

    X = np.arange(0, box_l) + 0.5
    x_vel = np.empty(box_l)

    while system.time <= total_time:

      # Compute analytical solution
      U = u(X, system.time, nu, v / time_step, box_l, k_max)

      # Read the data from the lattice nodes
      for i in np.arange(0, box_l):
        w = lb[0, i, 0].velocity
        x_vel[i] = w[0]

      # Compare deviation
      quad_dev = 0.0
      dev = U - x_vel

      for i in np.arange(0, len(x_vel)):
        quad_dev += dev[i] ** 2

      self.assertTrue(quad_dev < tol)

      system.integrator.run(1)
      system.lees_edwards_offset += v * time_step

    system.part.clear()
    system.actors.remove(lb_1)

if __name__ == "__main__":
  ut.main()
