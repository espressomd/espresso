#
# Copyright (C) 2010-2022 The ESPResSo project
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
import espressomd.lb
import unittest as ut
import unittest_decorators as utx
import numpy as np
import tests_common


class SwimmerTest():
    system = espressomd.System(box_l=3 * [6])
    system.cell_system.skin = 1
    system.time_step = 1e-2
    LB_params = {'agrid': 1,
                 'density': 1.1,
                 'viscosity': 1.2,
                 'kT': 0,
                 'tau': system.time_step}
    gamma = 0.3
    n_nodes = system.cell_system.get_state()['n_nodes']

    def add_all_types_of_swimmers(
            self,
            fix=False,
            rotation=False,
            put_in_corners=True):
        """Places all combinations of pusher/puller and f_swim/v_swim
        in a box, either in the corners or around the center
        """
        system = self.system
        plus_x = np.sqrt([.5, 0, .5, 0])
        plus_y = np.sqrt([0, 0, .5, .5])
        plus_z = np.sqrt([.5, 0, 0, .5])
        minus_y = np.sqrt([.5, .5, 0, 0])

        pos0 = [2, 0.01, 3] if put_in_corners else [2, 3, 2.5]
        pos1 = [5.99, 2, 3] if put_in_corners else [3.1, 2.1, 2.2]
        pos2 = [2, 3, 5.99] if put_in_corners else [2.9, 2.5, 3]
        pos3 = [1.5, 5.99, 1] if put_in_corners else [2, 2, 2.5]

        # particle type is only relevant for engine_lb_hybrid test
        system.part.add(pos=pos0, quat=minus_y, fix=3 * [fix],
                        mass=0.9, rinertia=3 * [7], rotation=3 * [rotation],
                        swimming={"mode": "pusher", "f_swim": 0.10,
                                  "dipole_length": 0.5}, type=1)
        system.part.add(pos=pos1, quat=plus_x, fix=3 * [fix],
                        mass=1.9, rinertia=3 * [8], rotation=3 * [rotation],
                        swimming={"mode": "pusher", "v_swim": 0.02,
                                  "dipole_length": 0.6}, type=0)
        system.part.add(pos=pos2, quat=plus_z, fix=3 * [fix],
                        mass=2.9, rinertia=3 * [9], rotation=3 * [rotation],
                        swimming={"mode": "puller", "f_swim": 0.08,
                                  "dipole_length": 0.7}, type=1)
        system.part.add(pos=pos3, quat=plus_y, fix=3 * [fix],
                        mass=3.9, rinertia=3 * [10], rotation=3 * [rotation],
                        swimming={"mode": "puller", "v_swim": 0.05,
                                  "dipole_length": 0.8}, type=0)

    def setUp(self):
        self.set_cellsystem()
        self.lbf = self.lb_class(**self.LB_params, **self.lb_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=self.gamma)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test_momentum_conservation(self):
        """friction as well as 'active' forces apply to particles
        and to the fluid, so total momentum is conserved
        """
        if self.lbf.get_params().get('single_precision', False):
            self.skipTest('Momentum is not conserved on single precision')

        self.add_all_types_of_swimmers(rotation=False)

        # Comments by Christoph Lohrmann from #3514:
        # - why I used `reuse_forces=True` : If I don't use it, `force_calc()`
        #   is called the first time without LB-coupling. That means no friction
        #   for any swimmer and no additional force for the `v_swim` type swimmers.
        #   The active force for the `f_swim` swimmers gets added anyway because
        #   it is not derived from the coupling to LB. With `reuse_forces` at
        #   least both types are treated the same.
        # - Therefore, in the first halfstep, the active forces on the particles
        #   are missing. This creates the first half of the missing momentum.
        # - The LB fluid is always ahead by a half step (as checked by
        #   `test_ext_force_density()` in `lb.py`). It is also not affected by
        #   the `reuse_forces` in the first halfstep because `force_calc()`
        #   with coupling is called in the main integration loop before
        #   `lb_lbfluid_propagate()`
        # - => in total, the fluid momentum is ahead by a full time step
        self.system.integrator.run(100, reuse_forces=True)
        tot_mom = self.system.analysis.linear_momentum(include_particles=True,
                                                       include_lbfluid=True)
        # compensate offset between force calculation and LB-update
        for part in self.system.part:
            tot_mom += part.f * self.system.time_step

        np.testing.assert_allclose(tot_mom, 3 * [0.], atol=self.tol)

    def test_particle_forces(self):
        """run through all swimmers to check expected forces
        """
        self.add_all_types_of_swimmers(rotation=False)
        self.system.integrator.run(10)
        for swimmer in self.system.part:

            f_swim = swimmer.swimming['f_swim']
            v_swim = swimmer.swimming['v_swim']
            director = swimmer.director

            # due to dt/2 time-shift between force calculation and LB-update,
            # v_swimmer has to be calculated at the half step
            v_swimmer = swimmer.v + \
                0.5 * self.system.time_step * swimmer.f / swimmer.mass
            # for friction coupling, the old fluid at the new position is used
            v_fluid = self.lbf.get_interpolated_velocity(
                swimmer.pos + self.system.time_step * v_swimmer)
            force = -self.gamma * (v_swimmer - v_fluid) + \
                f_swim * director + self.gamma * v_swim * director

            self.system.integrator.run(1, reuse_forces=True)
            np.testing.assert_allclose(
                np.copy(swimmer.f), force, atol=self.tol)

    def test_fluid_force(self):
        """ forces on particles are already checked (self.test_particle_forces)
        total force on the fluid matches (self.test_momentum_conservation)
        only thing left to check is the location of the fluid force.
        """
        f_swim = 0.11
        dip_length = 2 * self.LB_params['agrid']

        sw0_pos = np.array([3.8, 1.1, 1.1])
        sw1_pos = np.array([1.8, 4.1, 4.1])
        sw0 = self.system.part.add(pos=sw0_pos, quat=np.sqrt([.5, 0, .5, 0]),
                                   mass=0.9, rotation=3 * [False],
                                   swimming={"mode": "pusher", "f_swim": f_swim,
                                             "dipole_length": dip_length})
        sw1 = self.system.part.add(pos=sw1_pos, quat=np.sqrt([.5, 0, .5, 0]),
                                   mass=0.7, rotation=3 * [False],
                                   swimming={"mode": "puller", "f_swim": f_swim,
                                             "dipole_length": dip_length})

        self.system.integrator.run(2)

        for sw in [sw0, sw1]:
            push_pull = -1 if sw.swimming['mode'] == 'pusher' else 1
            sw_source_pos = sw.pos + self.system.time_step * \
                sw.v + push_pull * dip_length * sw.director
            # fold into box
            sw_source_pos -= np.floor(sw_source_pos /
                                      self.system.box_l) * np.array(self.system.box_l)
            sw_source_nodes = tests_common.get_lb_nodes_around_pos(
                sw_source_pos, self.lbf)
            sw_source_forces = np.array(
                [n.last_applied_force for n in sw_source_nodes])
            np.testing.assert_allclose(
                np.sum(sw_source_forces, axis=0),
                -f_swim * np.array(sw.director), atol=self.tol)


@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestDomDecWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestDomDecWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestNSquareWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestNSquareWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestHybrid0CPUWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestHybrid0CPUWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestHybrid1CPUWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
@utx.skipIfMissingFeatures(
    ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA"])
class SwimmerTestHybrid1CPUWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


if __name__ == "__main__":
    ut.main()
