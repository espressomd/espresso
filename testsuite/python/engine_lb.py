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
import espressomd.virtual_sites
import espressomd.swimmer_helpers
import unittest as ut
import unittest_decorators as utx
import numpy as np
import tests_common


REQUIRED_FEATURES = ["ENGINE", "ROTATIONAL_INERTIA", "MASS", "WALBERLA",
                     "VIRTUAL_SITES_RELATIVE"]


class SwimmerTest:
    system = espressomd.System(box_l=3 * [8])
    system.cell_system.skin = 0.4
    system.time_step = 1e-2
    LB_params = {
        "agrid": 1,
        "density": 1.1,
        "kinematic_viscosity": 1.2,
        "kT": 0,
        "tau": system.time_step,
    }
    gamma = 0.3

    # largest dipole used in the test
    system.min_global_cut = 2.0001 * LB_params["agrid"]

    def add_all_types_of_swimmers(self, fix=False, rotation=False):
        """
        Places pushers and pullers in the corners of a box such that the
        dipole crosses the periodic boundary.
        """
        system = self.system
        dipole_partcl_type = 10

        # particle type is only relevant for engine_lb_hybrid test
        p0 = system.part.add(
            pos=[2, 0.01, 3],
            director=[0, -1, 0],
            fix=3 * [fix],
            mass=0.1,
            rinertia=3 * [7],
            rotation=3 * [rotation],
            swimming={"f_swim": 0.10},
            type=1,
        )
        p0_dip = espressomd.swimmer_helpers.add_dipole_particle(
            system, p0, 0.5, dipole_partcl_type)

        p1 = system.part.add(
            pos=[7.99, 2, 3],
            director=[1, 0, 0],
            fix=3 * [fix],
            mass=0.2,
            rinertia=3 * [8],
            rotation=3 * [rotation],
            swimming={"f_swim": 0.09},
            type=0,
        )
        p1_dip = espressomd.swimmer_helpers.add_dipole_particle(
            system, p1, 0.6, dipole_partcl_type)

        p2 = system.part.add(
            pos=[2, 3, 7.99],
            director=[0, 0, 1],
            fix=3 * [fix],
            mass=0.3,
            rinertia=3 * [9],
            rotation=3 * [rotation],
            swimming={"f_swim": 0.08},
            type=1,
        )
        p2_dip = espressomd.swimmer_helpers.add_dipole_particle(
            system, p2, 0.7, dipole_partcl_type, mode="puller")

        p3 = system.part.add(
            pos=[1.5, 7.99, 1],
            director=[0, 1, 0],
            fix=3 * [fix],
            mass=0.4,
            rinertia=3 * [10],
            rotation=3 * [rotation],
            swimming={"f_swim": 0.07},
            type=0,
        )
        p3_dip = espressomd.swimmer_helpers.add_dipole_particle(
            system, p3, 0.8, dipole_partcl_type, mode="puller")

        return [p0, p1, p2, p3], [p0_dip, p1_dip, p2_dip, p3_dip]

    def setUp(self):
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
            have_quaternion=True
        )
        self.set_cellsystem()
        self.lbf = self.lb_class(**self.LB_params, **self.lb_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=self.gamma)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test_dipole_particle_addition(self):
        swimmers, dipoles = self.add_all_types_of_swimmers(rotation=False)
        self.system.integrator.run(0)
        for sw, dip in zip(swimmers, dipoles):
            np.testing.assert_allclose(
                np.copy(dip.director), -np.copy(sw.director), atol=1e-15
            )
        with self.assertRaises(ValueError):
            espressomd.swimmer_helpers.add_dipole_particle(
                self.system, swimmers[0], 0.5, 1, mode="unknown_mode")

    def test_momentum_conservation(self):
        """
        friction as well as 'active' forces apply to particles
        and to the fluid, so total momentum is conserved
        """
        if self.lbf.get_params().get("single_precision", False):
            self.skipTest("Momentum is not conserved on single precision")

        swimmers, _ = self.add_all_types_of_swimmers(rotation=False)

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
        self.system.integrator.run(1000, reuse_forces=True)
        lb_mom = self.system.analysis.linear_momentum(
            include_particles=False, include_lbfluid=True
        )
        part_mom = np.sum([np.copy(s.v) * np.copy(s.mass)
                          for s in swimmers], axis=0)

        # compensate for one timestep offset between force calculation and LB -
        # update
        for part in swimmers:
            part_mom += np.copy(part.f) * self.system.time_step

        np.testing.assert_allclose(lb_mom, -part_mom, rtol=self.tol)

    def test_particle_forces(self):
        """run through all swimmers to check expected forces"""
        swimmers, _ = self.add_all_types_of_swimmers(rotation=False)
        self.system.integrator.run(100)
        for swimmer in swimmers:
            f_swim = swimmer.swimming["f_swim"]
            director = swimmer.director

            # due to dt/2 time-shift between force calculation and LB-update,
            # v_swimmer has to be calculated at the half step
            v_swimmer = (
                swimmer.v + 0.5 * self.system.time_step * swimmer.f / swimmer.mass
            )
            # for friction coupling, the old fluid at the new position is used
            v_fluid = self.lbf.get_interpolated_velocity(pos=swimmer.pos + self.system.time_step * v_swimmer
                                                         )
            force = -self.gamma * (v_swimmer - v_fluid) + f_swim * director

            self.system.integrator.run(1, reuse_forces=True)
            np.testing.assert_allclose(
                np.copy(swimmer.f), force, atol=self.tol)

    def test_fluid_force(self):
        """
        forces on particles are already checked (self.test_particle_forces)
        total force on the fluid matches (self.test_momentum_conservation)
        only thing left to check is the location of the fluid force.
        """
        f_swim = 0.11
        dip_length = 2 * self.LB_params["agrid"]
        modes = ["pusher", "puller"]

        sw0_pos = np.array([3.8, 1.1, 1.1])
        sw1_pos = np.array([1.8, 4.1, 4.1])

        sw_0 = self.system.part.add(
            pos=sw0_pos,
            director=[1, 0, 0],
            mass=0.9,
            rotation=3 * [False],
            swimming={"f_swim": f_swim},
        )
        espressomd.swimmer_helpers.add_dipole_particle(
            self.system, sw_0, dip_length, 10, mode=modes[0])
        sw_1 = self.system.part.add(
            pos=sw1_pos,
            director=[1, 0, 0],
            mass=0.7,
            rotation=3 * [False],
            swimming={"f_swim": f_swim},
        )
        espressomd.swimmer_helpers.add_dipole_particle(
            self.system, sw_1, dip_length, 10, mode=modes[1])

        self.system.integrator.run(2)

        for sw, mode in zip([sw_0, sw_1], modes):
            dip_sign = 1 if mode == "puller" else -1
            sw_source_pos = sw.pos + self.system.time_step * \
                sw.v + dip_sign * dip_length * sw.director
            # fold into box
            sw_source_pos -= np.floor(sw_source_pos / self.system.box_l) * np.array(
                self.system.box_l
            )
            sw_source_nodes = tests_common.get_lb_nodes_around_pos(
                sw_source_pos, self.lbf
            )
            sw_source_forces = np.array(
                [n.last_applied_force for n in sw_source_nodes])
            np.testing.assert_allclose(
                np.sum(sw_source_forces, axis=0),
                -f_swim * np.array(sw.director),
                atol=self.tol,
            )


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestDomDecWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestDomDecWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestNSquareWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestNSquareWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestHybrid0CPUWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestHybrid0CPUWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestHybrid1CPUWalberla(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


@utx.skipIfMissingFeatures(REQUIRED_FEATURES)
class SwimmerTestHybrid1CPUWalberlaSinglePrecision(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


if __name__ == "__main__":
    ut.main()
