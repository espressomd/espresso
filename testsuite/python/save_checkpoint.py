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
import numpy as np
import os

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.interactions
import espressomd.virtual_sites
import espressomd.accumulators
import espressomd.observables
import espressomd.lb
import espressomd.electrokinetics

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

# use a box with 3 different dimensions
system = espressomd.System(box_l=[12.0, 14.0, 16.0])
system.cell_system.skin = 0.1
system.seed = system.cell_system.get_state()["n_nodes"] * [1234]
system.time_step = 0.01
system.min_global_cut = 2.0

# create checkpoint folder
idx = "mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace(".", "__")
checkpoint = espressomd.checkpointing.Checkpoint(
    checkpoint_id=idx, checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")

LB_implementation = None
if 'LB.CPU' in modes:
    LB_implementation = espressomd.lb.LBFluid
elif espressomd.gpu_available() and espressomd.has_features('CUDA') and 'LB.GPU' in modes:
    LB_implementation = espressomd.lb.LBFluidGPU
if LB_implementation:
    lbf = LB_implementation(agrid=0.5, visc=1.3, dens=1.5, tau=0.01)
    system.actors.add(lbf)

EK_implementation = None
if espressomd.gpu_available() and espressomd.has_features('ELECTROKINETICS') and 'EK.GPU' in modes:
    EK_implementation = espressomd.electrokinetics
    ek = EK_implementation.Electrokinetics(
        agrid=0.5,
          lb_density=26.15,
          viscosity=1.7,
          friction=0.0,
          T=1.1,
          prefactor=0.88,
          stencil="linkcentered")
    ek_species = EK_implementation.Species(
        density=0.4,
          D=0.02,
          valency=0.3,
          ext_force_density=[0.01, -0.08, 0.06])
    ek.add_species(ek_species)
    system.actors.add(ek)

system.part.add(pos=[1.0] * 3)
system.part.add(pos=[1.0, 1.0, 2.0])

if espressomd.has_features('EXCLUSIONS'):
    system.part.add(pos=[2.0] * 3, exclusions=[0, 1])

if espressomd.has_features('ELECTROSTATICS') and 'P3M.CPU' in modes:
    system.part[0].q = 1
    system.part[1].q = -1
    p3m = espressomd.electrostatics.P3M(
        prefactor=1.0,
        accuracy=0.1,
        mesh=10,
        cao=1,
        alpha=1.0,
        r_cut=1.0,
        tune=False)
    system.actors.add(p3m)

obs = espressomd.observables.ParticlePositions(ids=[0, 1])
acc = espressomd.accumulators.MeanVarianceCalculator(obs=obs)
acc.update()
system.part[0].pos = [1.0, 2.0, 3.0]
acc.update()

system.thermostat.set_langevin(kT=1.0, gamma=2.0, seed=42)

if espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']):
    system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
        have_velocity=True, have_quaternion=True)
    system.part[1].vs_auto_relate_to(0)

if espressomd.has_features(['LENNARD_JONES']) and 'LJ' in modes:
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.3, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[3, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.7, cutoff=2.0, shift=0.1)

harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0.0, k=1.0)
system.bonded_inter.add(harmonic_bond)
system.part[1].add_bond((harmonic_bond, 0))
checkpoint.register("system")
checkpoint.register("acc")
# calculate forces
system.integrator.run(0)
particle_force0 = np.copy(system.part[0].f)
particle_force1 = np.copy(system.part[1].f)
checkpoint.register("particle_force0")
checkpoint.register("particle_force1")
if espressomd.has_features("COLLISION_DETECTION"):
    system.collision_detection.set_params(
        mode="bind_centers", distance=0.11, bond_centers=harmonic_bond)

if LB_implementation:
    m = np.pi / 12
    nx = int(np.round(system.box_l[0] / lbf.get_params()["agrid"]))
    ny = int(np.round(system.box_l[1] / lbf.get_params()["agrid"]))
    nz = int(np.round(system.box_l[2] / lbf.get_params()["agrid"]))
    # Create a 3D grid with deterministic values to fill the LB fluid lattice
    grid_3D = np.fromfunction(
        lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
                              (nx, ny, nz), dtype=float)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                lbf[i, j, k].population = grid_3D[i, j, k] * np.arange(1, 20)
    cpt_mode = int("@TEST_BINARY@")
    # save LB checkpoint file
    lbf_cpt_path = checkpoint.checkpoint_dir + "/lb.cpt"
    lbf.save_checkpoint(lbf_cpt_path, cpt_mode)

if EK_implementation:
    m = np.pi / 12
    nx = int(np.round(system.box_l[0] / ek.get_params()["agrid"]))
    ny = int(np.round(system.box_l[1] / ek.get_params()["agrid"]))
    nz = int(np.round(system.box_l[2] / ek.get_params()["agrid"]))
    # Create a 3D grid with deterministic values to fill the LB fluid lattice
    grid_3D = np.fromfunction(
        lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
                              (nx, ny, nz), dtype=float)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ek_species[i, j, k].density = grid_3D[i, j, k]
    # save LB checkpoint file
    ek_cpt_path = checkpoint.checkpoint_dir + "/ek"
    ek.save_checkpoint(ek_cpt_path)

# save checkpoint file
checkpoint.save(0)


class TestLB(ut.TestCase):

    def test_checkpointing(self):
        self.assertTrue(os.path.isdir(checkpoint.checkpoint_dir),
                        "checkpoint directory not created")

        checkpoint_filepath = checkpoint.checkpoint_dir + "/0.checkpoint"
        self.assertTrue(os.path.isfile(checkpoint_filepath),
                        "checkpoint file not created")

        if LB_implementation:
            self.assertTrue(os.path.isfile(lbf_cpt_path),
                            "LB checkpoint file not created")

            with open(lbf_cpt_path, "rb") as f:
                lbf_cpt_str = f.read()
            # write an LB checkpoint with missing data
            with open(lbf_cpt_path[:-4] + "-corrupted.cpt", "wb") as f:
                f.write(lbf_cpt_str[:len(lbf_cpt_str) // 2])
            # write an LB checkpoint with different box dimensions
            with open(lbf_cpt_path[:-4] + "-wrong-boxdim.cpt", "wb") as f:
                if cpt_mode == 1:
                    # first dimension becomes 0
                    f.write(8 * b"\x00" + lbf_cpt_str[8:])
                else:
                    # first dimension becomes larger
                    f.write(b"1" + lbf_cpt_str)


if __name__ == '__main__':
    ut.main()
