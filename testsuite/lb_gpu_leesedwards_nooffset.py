#!/usr/bin/env python

from __future__ import print_function

import espressomd as md
from espressomd import lb
import numpy as np
import unittest as ut
from tests_common import abspath

# Systemclass
system = md.System(box_l=[9,9,9])

# allowed deviation
tol = 1.0e-16


@ut.skipIf(not md.has_features(['LEES_EDWARDS', 'LB_GPU', 'EXTERNAL_FORCES']),
  'Features not available, skipping test!')
class LeesEdwardsOffsetTest(ut.TestCase):

  def test(self):

    """The velocity components at all LB nodes are tested against a table stored in testsuite/data/""" 

    # LB parameter
    box_l = 9
    eta = 1.0
    rho = 1.0
    nu = eta / rho
    v = 1.0
    k_max = 100
    agrid = 1.0
    system.box_l = [box_l, box_l, box_l]

    # Integration parameters
    time_step = 1.0
    system.time_step = time_step
    system.cell_system.skin = 0.4

    # Add a fixed particle with a velovity in the y-direction
    system.part.add(pos=[4.5, 7.5, 4.5], fix=[1, 1, 1], id=0, type=0)
    system.part[0].v = [0, 0.1, 0]

    system.lees_edwards_offset = 0.0
    system.integrator.run(2)

    # Add LB fluid
    lb = md.lb.LBFluidGPU(agrid=agrid, dens=rho, visc=eta, fric=1.0, tau=time_step)
    system.actors.add(lb)

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
    self.assertTrue(np.allclose(node_list_pos_v, saved_data, rtol=tol, atol=tol), 'populations do notmatch')

if __name__ == "__main__":
  ut.main()
