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
import espressomd.shapes

AGRID = 1.
BOX_L = 8. * AGRID
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


class LBCheckpointThermalized:
    """
    Check that a thermalized LB fluid restored from a checkpoint
    continues exactly like the uninterrupted simulation.
    """

    system = espressomd.System(box_l=3 * [BOX_L])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()
        self.system.lb = None

    def make_fluid(self, with_walls):
        lbf = espressomd.lb.LBFluid(**LB_PARAMS, **self.lb_params)
        if with_walls:
            # boundary nodes on both sides of the periodic image, such that
            # the ghost layers contain boundary nodes too; one wall moves,
            # such that the boundary PDFs depend on the fluid density
            lbf.add_boundary_from_shape(
                espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5))
            lbf.add_boundary_from_shape(
                espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(BOX_L - 1.)),
                velocity=[0., 0.02, -0.01])
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

    def check_checkpoint(self, with_walls, with_particles, binary):
        system = self.system
        lbf = self.make_fluid(with_walls)
        system.lb = lbf
        if with_particles:
            # particles are far apart so that they never couple to the
            # same LB node, which makes the force summation order-independent
            pos = np.array([[4.2, 1.1, 1.3], [4.8, 3.9, 4.6], [3.3, 6.4, 6.7]])
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
            if with_particles:
                saved_particles = self.get_particle_state()
                saved_philox_counter = system.thermostat.lb.philox_counter

            # uninterrupted reference run
            system.integrator.run(N_CONTINUE, reuse_forces=True)
            ref_fluid = self.get_fluid_state(lbf)
            if with_particles:
                ref_particles = self.get_particle_state()

            # restart from the checkpoint into a new fluid
            system.thermostat.turn_off()
            system.lb = None
            lbf = self.make_fluid(with_walls)
            system.lb = lbf
            lbf.load_checkpoint(str(path), binary)
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


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCheckpointWalberlaDoublePrecisionCPU(LBCheckpointThermalized,
                                             ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}
    rtol = 1e-12
    atol = 1e-12


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCheckpointWalberlaSinglePrecisionCPU(LBCheckpointThermalized,
                                             ut.TestCase):
    lb_params = {"single_precision": True, "gpu": False}
    rtol = 1e-6
    atol = 1e-6


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBCheckpointWalberlaDoublePrecisionGPU(LBCheckpointThermalized,
                                             ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}
    rtol = 1e-12
    atol = 1e-12


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBCheckpointWalberlaSinglePrecisionGPU(LBCheckpointThermalized,
                                             ut.TestCase):
    lb_params = {"single_precision": True, "gpu": True}
    rtol = 1e-6
    atol = 1e-6


if __name__ == "__main__":
    ut.main()
