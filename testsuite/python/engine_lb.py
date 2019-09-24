# Copyright (C) 2010-2018 The ESPResSo project
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
from espressomd import System, lb, shapes, lbboundaries


class SwimmerTest():
    system = System(box_l=3 * [6])
    system.cell_system.skin = 1
    system.time_step = 1e-2
    LB_params = {'agrid': 1,
                 'dens': 1.1,
                 'visc': 1.2,
                 'kT': 0,
                 'tau': system.time_step}    
    gamma = 0.3
    lbf = None

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
                                  "dipole_length": 0.5, "rotational_friction": 0.3})
        system.part.add(pos=pos1, quat=plus_x, fix=3 * [fix],
                        mass=1.9, rinertia=3 * [8], rotation=3 * [rotation],
                        swimming={"mode": "pusher", "v_swim": 0.02,
                                  "dipole_length": 0.6, "rotational_friction": 0.4})
        system.part.add(pos=pos2, quat=plus_z, fix=3 * [fix],
                        mass=2.9, rinertia=3 * [9], rotation=3 * [rotation],
                        swimming={"mode": "puller", "f_swim": 0.08,
                                  "dipole_length": 0.7, "rotational_friction": 0.8})
        system.part.add(pos=pos3, quat=plus_y, fix=3 * [fix],
                        mass=3.9, rinertia=3 * [10], rotation=3 * [rotation],
                        swimming={"mode": "puller", "v_swim": 0.05,
                                  "dipole_length": 0.8, "rotational_friction": 0.3})

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.lbboundaries.clear()
        self.system.thermostat.turn_off()
        self.lbf = None

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
            np.testing.assert_allclose(swimmer.f, force, atol=self.tol)

    def check_fluid_force(self, swimmer):
        pass
        # forces on particles are checked
        # total force on the fluid matches (momentum conservation)
        # TODO: only thing left to check is the location of the fluid force.
        # There is no counter torque when using the rotational_friction-feature
        # so there is nothing to be tested


@utx.skipIfMissingFeatures(["ENGINE", "ROTATION", "MASS"])
class SwimmerTestCPU(SwimmerTest, ut.TestCase):    

    def setUp(self):
        self.tol = 1e-10
        self.lbf = lb.LBFluid(**self.LB_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=self.gamma)

    def test_rotfric_exception(self):
        """rotational_friction feature is disabled on CPU for more than one core
        """
        if self.system.cell_system.get_state()["n_nodes"] > 1:
            # swimming without rot_fric is fine
            self.system.part.add(pos=[0, 0, 0], rotation=3 * [True],
                                 swimming={"f_swim": 0.1, "mode": "pusher"})
            self.system.integrator.run(3)
            # with rot_fric it is not
            with self.assertRaises(Exception):
                self.system.part.add(pos=[0, 0, 0], rotation=3 * [True],
                                     swimming={"f_swim": 0.1, 
                                               "mode": "pusher", 
                                               "rotational_friction": 0.3})
                self.system.integrator.run(3)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["ENGINE", "ROTATION", "MASS"])
class SwimmerTestGPU(SwimmerTest, ut.TestCase):   

    def setUp(self):
        self.tol = 1e-5
        self.lbf = lb.LBFluidGPU(**self.LB_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=self.gamma)

    def test_particle_torques(self):
        """setup shear flow and check if resulting torques match 
        the formulae in the core
        """

        bottom = shapes.Wall(normal=[0, 0, 1], 
                             dist=self.LB_params['agrid'])
        top = shapes.Wall(normal=[0, 0, -1],
                          dist=-self.system.box_l[2] + self.LB_params['agrid'])
        self.system.lbboundaries.add(lbboundaries.LBBoundary(shape=bottom))
        self.system.lbboundaries.add(
            lbboundaries.LBBoundary(shape=top, velocity=[1e-3, 1e-3, 0]))
        self.system.integrator.run(100)

        # fix the particles so inaccuracies from position updates
        # before torque calculation don't matter
        self.add_all_types_of_swimmers(fix=True, rotation=True,
                                       put_in_corners=False)

        self.system.integrator.run(20)
        for swimmer in self.system.part:

            director = swimmer.director
            dip_len = swimmer.swimming["dipole_length"]
            mode_fac = 1. if swimmer.swimming["mode"] == "puller" else -1.
            source_pos = swimmer.pos + mode_fac * dip_len * director 

            v_center = self.lbf.get_interpolated_velocity(swimmer.pos)
            v_source = self.lbf.get_interpolated_velocity(source_pos)
            diff = v_center - v_source
            cross = np.cross(diff, director) 
            # half-step omega with isotropic rinertia
            omega_part = swimmer.omega_lab + 0.5 * self.system.time_step * \
                swimmer.torque_lab / swimmer.rinertia[0]
            omega_swim = cross / np.linalg.norm(cross) * \
                np.linalg.norm(diff) / dip_len
            torque = swimmer.swimming["rotational_friction"] * \
                (omega_swim - omega_part)

            self.system.integrator.run(1, reuse_forces=True)
            np.testing.assert_allclose(
                swimmer.torque_lab,
                torque,
                atol=self.tol)    


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(SwimmerTestGPU))
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(SwimmerTestCPU))
    result = ut.TextTestRunner(verbosity=4).run(suite)
