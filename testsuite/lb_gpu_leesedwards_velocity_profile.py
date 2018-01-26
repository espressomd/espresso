#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import unittest as ut
import espressomd as md
from espressomd import lb


@ut.skipIf(not md.has_features(['LB_GPU', 'LEES_EDWARDS']),
  'Features not available, skipping test!')
class LBGPULeesEdwardsTest(ut.TestCase):
  def test(self):

    """In this test the velocity profile of a LB fluid under steady shear is compared to the velocity profile that is obtained by solving the Navier-Stokes equation with Fourier series."""

    # Analytical solution with Fourier series of Navier-Stikes equation
    def u(x, t, nu, v, h, k_max):
      u = x / h - 0.5
      for k in np.arange(1, k_max + 1):
        u += 1.0 / (np.pi * k) * np.exp(-4 * np.pi ** 2 * nu * k ** 2 / h ** 2 * t) * np.sin(2 * np.pi / h * k * x)
      return v * u

    # Systemclass
    system = md.System()

    tol = 1.0e-5

    # LB parameter
    box_l = 9
    eta = 0.27
    rho = 0.97
    nu = eta / rho
    v = 0.0087
    k_max = 100
    system.box_l = [box_l, box_l, box_l]

    # Integration parameters
    time_step = 1.0
    total_time = 1000.0
    system.time_step = time_step
    system.cell_system.skin = 0.4

    # Add LB fluid
    lb = md.lb.LBFluidGPU(agrid=1.0, dens=rho, visc=eta, fric=1.1, tau=time_step)
    system.actors.add(lb)

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

if __name__ == "__main__":
  ut.main()
