#
# Copyright (C) 2026 The ESPResSo project
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

import itertools
import pathlib
import tempfile
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
import espressomd.lees_edwards
import espressomd.shapes

AGRID = 1.
BOX_L = 12. * AGRID
TIME_STEP = 0.01
LB_PARAMS = {"agrid": AGRID,
             "density": 0.9,
             "kinematic_viscosity": 1.3,
             "tau": TIME_STEP,
             "kT": 0.8,
             "seed": 41}
GAMMA = 2.
N_WARMUP = 20
N_CONTINUE = 10
# the offset is not a multiple of AGRID, such that the Lees-Edwards ghost
# layers are interpolated between two nodes, and it changes over time
LE_PROTOCOL_PARAMS = {"initial_pos_offset": 0.7 * AGRID,
                      "shear_velocity": 1.5,
                      "time_0": 0.}


class LBCheckpoint:
    """
    Check that an LB fluid restored from a checkpoint continues exactly
    like the uninterrupted simulation, for a thermalized fluid and for
    an unthermalized fluid with Lees-Edwards boundary conditions.
    """

    system = espressomd.System(box_l=3 * [BOX_L])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()
        self.system.lb = None
        self.system.lees_edwards.protocol = None
        self.system.time = 0.

    def make_fluid(self, with_walls, kT, wall_axis=0):
        lbf = espressomd.lb.LBFluid(**{**LB_PARAMS, "kT": kT},
                                    **self.lb_params)
        if with_walls:
            # boundary nodes on both sides of the periodic image, such that
            # the ghost layers contain boundary nodes too; one wall moves
            # tangentially, such that the boundary PDFs depend on the fluid
            # density
            normal = np.zeros(3)
            normal[wall_axis] = 1.
            velocity = np.roll([0., 0.02, -0.01], wall_axis)
            lbf.add_boundary_from_shape(
                espressomd.shapes.Wall(normal=normal, dist=1.5))
            lbf.add_boundary_from_shape(
                espressomd.shapes.Wall(normal=-normal, dist=-(BOX_L - 1.)),
                velocity=velocity)
        return lbf

    def get_fluid_state(self, lbf):
        return {"velocity": np.copy(lbf[:, :, :].velocity),
                "population": np.copy(lbf[:, :, :]._population),
                "last_applied_force":
                    np.copy(lbf[:, :, :].last_applied_force),
                "boundary_force": np.copy(lbf.boundary_force)}

    def get_particle_state(self):
        return {"pos": np.copy(self.system.part.all().pos),
                "v": np.copy(self.system.part.all().v),
                "f": np.copy(self.system.part.all().f)}

    def check_checkpoint(self, with_walls, with_particles, binary,
                         with_lees_edwards=False):
        system = self.system
        if with_lees_edwards:
            # Lees-Edwards LB is unthermalized and only supports a shear
            # along x with the shear plane normal along y, so the walls
            # are orthogonal to z
            kT = 0.
            wall_axis = 2
            protocol = espressomd.lees_edwards.LinearShear(
                **LE_PROTOCOL_PARAMS)
            system.lees_edwards.set_boundary_conditions(
                shear_direction="x", shear_plane_normal="y",
                protocol=protocol)
        else:
            kT = LB_PARAMS["kT"]
            wall_axis = 0
        lbf = self.make_fluid(with_walls, kT, wall_axis)
        system.lb = lbf
        if with_particles:
            # particles are far apart so that they never couple to the
            # same LB node, which makes the force summation order-independent
            if with_lees_edwards:
                # particles stay clear of the walls along z; the last two
                # are less than AGRID / 2 away from the shear plane and
                # couple to the Lees-Edwards-shifted ghost layer
                pos = np.array([[4.2, 1.1, 3.3], [4.8, 3.9, 5.6],
                                [3.3, 6.4, 7.7], [6.1, 0.2, 9.4],
                                [8.6, BOX_L - 0.3, 4.4]])
                vel = np.array([[0.1, -0.2, 0.3], [-0.3, 0.1, 0.2],
                                [0.2, 0.3, -0.1], [-0.1, 0.1, 0.2],
                                [0.2, -0.1, -0.3]])
            else:
                pos = np.array([[4.2, 1.1, 1.3], [4.8, 3.9, 4.6],
                                [3.3, 6.4, 6.7]])
                vel = np.array([[0.1, -0.2, 0.3], [-0.3, 0.1, 0.2],
                                [0.2, 0.3, -0.1]])
            system.part.add(pos=pos, v=vel)
            system.thermostat.set_lb(LB_fluid=lbf, seed=17, gamma=GAMMA)
        system.integrator.run(N_WARMUP)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path = pathlib.Path(tmp_dir) / "lb.cpt"
            lbf.save_checkpoint(str(path), binary)
            saved_fluid = self.get_fluid_state(lbf)
            saved_rng_state = lbf.rng_state
            saved_time = system.time
            if with_particles:
                saved_particles = self.get_particle_state()
                saved_philox_counter = system.thermostat.lb.philox_counter

            # uninterrupted reference run
            system.integrator.run(N_CONTINUE, reuse_forces=True)
            ref_fluid = self.get_fluid_state(lbf)
            if with_particles:
                ref_particles = self.get_particle_state()

            # restart from the checkpoint into a new fluid; the simulation
            # time is not part of the LB checkpoint and must be restored
            # first, because the Lees-Edwards offset of the ghost layers
            # is evaluated at the time the checkpoint is loaded
            system.time = saved_time
            system.thermostat.turn_off()
            system.lb = None
            lbf = self.make_fluid(with_walls, kT, wall_axis)
            system.lb = lbf
            lbf.load_checkpoint(str(path), binary)
        if with_lees_edwards:
            self.assertIsNone(saved_rng_state)
        else:
            lbf.rng_state = saved_rng_state
        if with_particles:
            system.thermostat.set_lb(LB_fluid=lbf, seed=17, gamma=GAMMA)
            system.thermostat.lb.call_method(
                "override_philox_counter", counter=saved_philox_counter)
            for key, value in saved_particles.items():
                setattr(system.part.all(), key, value)

        # the restored fluid matches the saved fluid, including the velocity
        loaded_fluid = self.get_fluid_state(lbf)
        for key, value in saved_fluid.items():
            np.testing.assert_allclose(loaded_fluid[key], value,
                                       rtol=self.rtol, atol=self.atol,
                                       err_msg=f"{key} differs after load")

        # the restarted simulation continues the uninterrupted run;
        # the restored particle forces must be reused, because a force
        # recalculation would skip the LB coupling for the first time step
        system.integrator.run(N_CONTINUE, reuse_forces=True)
        new_fluid = self.get_fluid_state(lbf)
        for key, value in ref_fluid.items():
            np.testing.assert_allclose(new_fluid[key], value,
                                       rtol=self.rtol, atol=self.atol,
                                       err_msg=f"{key} differs after run")
        if with_particles:
            new_particles = self.get_particle_state()
            for key, value in ref_particles.items():
                np.testing.assert_allclose(new_particles[key], value,
                                           rtol=self.rtol, atol=self.atol,
                                           err_msg=f"particle {key} differs")

    def test_checkpoint(self):
        for with_walls, with_particles, binary in itertools.product(
                (False, True), (False, True), (False, True)):
            with self.subTest(walls=with_walls, particles=with_particles,
                              binary=binary):
                try:
                    self.check_checkpoint(with_walls, with_particles, binary)
                finally:
                    self.tearDown()

    def test_checkpoint_lees_edwards(self):
        if self.lb_params["gpu"]:
            self.skipTest("Lees-Edwards LB doesn't support GPU")
        if self.system.cell_system.node_grid[0] != 1:
            self.skipTest("Lees-Edwards LB doesn't support domain "
                          "decomposition along the shear direction")
        for with_walls, with_particles in itertools.product(
                (False, True), (False, True)):
            with self.subTest(walls=with_walls, particles=with_particles):
                try:
                    self.check_checkpoint(with_walls, with_particles, True,
                                          with_lees_edwards=True)
                finally:
                    self.tearDown()


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCheckpointWalberlaDoublePrecisionCPU(LBCheckpoint,
                                             ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}
    rtol = 1e-12
    atol = 1e-12


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBCheckpointWalberlaSinglePrecisionGPU(LBCheckpoint,
                                             ut.TestCase):
    lb_params = {"single_precision": True, "gpu": True}
    rtol = 1e-6
    atol = 1e-6


if __name__ == "__main__":
    ut.main()
