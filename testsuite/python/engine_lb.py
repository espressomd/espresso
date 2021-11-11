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
import numpy as np
import espressomd
import espressomd.lb


class SwimmerTest():
    system = espressomd.System(box_l=3 * [6])
    system.cell_system.skin = 1
    system.time_step = 1e-2
    LB_params = {'agrid': 1,
                 'dens': 1.1,
                 'visc': 1.2,
                 'kT': 0,
                 'tau': system.time_step}
    gamma = 0.3

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

        system.part.add(pos=pos0, quat=minus_y, fix=3 * [fix],
                        mass=0.9, rinertia=3 * [7], rotation=3 * [rotation],
                        swimming={"mode": "pusher", "f_swim": 0.10,
                                  "dipole_length": 0.5})
        system.part.add(pos=pos1, quat=plus_x, fix=3 * [fix],
                        mass=1.9, rinertia=3 * [8], rotation=3 * [rotation],
                        swimming={"mode": "pusher", "v_swim": 0.02,
                                  "dipole_length": 0.6})
        system.part.add(pos=pos2, quat=plus_z, fix=3 * [fix],
                        mass=2.9, rinertia=3 * [9], rotation=3 * [rotation],
                        swimming={"mode": "puller", "f_swim": 0.08,
                                  "dipole_length": 0.7})
        system.part.add(pos=pos3, quat=plus_y, fix=3 * [fix],
                        mass=3.9, rinertia=3 * [10], rotation=3 * [rotation],
                        swimming={"mode": "puller", "v_swim": 0.05,
                                  "dipole_length": 0.8})

    def setUp(self):
        self.lbf = self.lb_class(**self.LB_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=self.gamma)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test_conflicting_parameters(self):
        """v_swim and f_swim can't be set at the same time
        """
        swimmer = self.system.part.add(pos=[3] * 3)
        with self.assertRaises(Exception):
            swimmer.swimming = {"v_swim": 0.3, "f_swim": 0.6}

    def test_momentum_conservation(self):
        """friction as well as 'active' forces apply to particles
        and to the fluid, so total momentum is conserved
        """
        self.add_all_types_of_swimmers(rotation=False)
        self.system.integrator.run(20, reuse_forces=True)
        tot_mom = self.system.analysis.linear_momentum(include_particles=True,
                                                       include_lbfluid=True)
        # compensate half-step offset between force calculation and LB-update
        for part in self.system.part:
            tot_mom += part.f * self.system.time_step / 2.

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


@utx.skipIfMissingFeatures(["ENGINE", "ROTATIONAL_INERTIA", "MASS"])
class SwimmerTestCPU(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluid
    tol = 1e-10


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["ENGINE", "ROTATIONAL_INERTIA", "MASS"])
class SwimmerTestGPU(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidGPU
    tol = 1e-5


if __name__ == "__main__":
    ut.main()
