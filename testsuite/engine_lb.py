from __future__ import print_function

import unittest as ut
import tests_common

import espressomd
import espressomd.lb

import os
import numpy as np
try:
    import vtk
except ImportError:
    print("Module \"vtk\" not available, skipping test!")
    exit()


@ut.skipIf(not espressomd.has_features(["ENGINE", "LB"]),
           "Features not available, skipping test!")
class SwimmerTest(ut.TestCase):
    def prepare(self, S):
        boxl = 12
        tstep = 0.01

        S.box_l = [boxl, boxl, boxl]
        S.cell_system.skin = 0.1
        S.time_step = tstep

        S.part.add(id=0, pos=[6.0, 3.0, 2.0],
                   swimming={"mode": "pusher", "v_swim": 0.10,
                             "dipole_length": 1.0, "rotational_friction": 2.0},
                   quat=[np.sqrt(.5), np.sqrt(.5), 0, 0])
        S.part.add(
            id=1,
            pos=[
                2.0,
                3.0,
                6.0],
            swimming={
                "mode": "pusher",
                "f_swim": 0.03,
                "dipole_length": 2.0,
                "rotational_friction": 20.0},
            quat=[
                np.sqrt(.5),
                0,
                np.sqrt(.5),
                0])
        S.part.add(
            id=2,
            pos=[
                3.0,
                2.0,
                6.0],
            swimming={
                "mode": "puller",
                "v_swim": 0.15,
                "dipole_length": 0.5,
                "rotational_friction": 15.0},
            quat=[
                np.sqrt(.5),
                0,
                0,
                np.sqrt(.5)])
        S.part.add(id=3, pos=[3.0, 6.0, 2.0],
                   swimming={"mode": "puller", "f_swim": 0.05,
                             "dipole_length": 1.5, "rotational_friction": 6.0},
                   quat=[0, 0, np.sqrt(.5), np.sqrt(.5)])
        S.part[:].rotation=1,1,1

    def run_and_check(self, S, lbm, vtk_name):
        S.integrator.run(self.sampsteps)

        if self.new_configuration:
            lbm.print_vtk_velocity(vtk_name)
            self.assertTrue(True)
        else:
            lbm.print_vtk_velocity("engine_test_tmp.vtk")
            different, difference = tests_common.calculate_vtk_max_pointwise_difference(
                vtk_name, "engine_test_tmp.vtk", tol=2.0e-7)
            os.remove("engine_test_tmp.vtk")
            print(
                "Maximum deviation to the reference point is: {}".format(difference))
            self.assertTrue(different)

    def test(self):
        self.new_configuration = False
        self.sampsteps = 2000

        S = espressomd.System(box_l=[1.0, 1.0, 1.0])
        S.seed = S.cell_system.get_state()['n_nodes'] * [1234]
        self.prepare(S)

        lbm = espressomd.lb.LBFluid(
            agrid=1.0, tau=S.time_step, fric=0.5, visc=1.0, dens=1.0)
        S.actors.add(lbm)

        if (S.cell_system.get_state()["n_nodes"] > 1):
            print("NOTE: Ignoring testcase for n_nodes > 1")
        else:
            self.run_and_check(S, lbm, tests_common.abspath("data/engine_lb.vtk"))


if __name__ == '__main__':
    ut.main()
