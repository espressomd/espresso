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
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.interactions
import espressomd.virtual_sites
import espressomd.accumulators
import espressomd.observables
import espressomd.lb

checkpoint = espressomd.checkpointing.Checkpoint(
    checkpoint_id="mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace('.', '__'),
    checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

system = espressomd.System(box_l=[12.0, 12.0, 12.0])
system.cell_system.skin = 0.1
system.seed = system.cell_system.get_state()["n_nodes"] * [1234]
system.time_step = 0.01
system.min_global_cut = 2.0

LB_implementation = None
if espressomd.has_features('LB') and 'LB.CPU' in modes:
    LB_implementation = espressomd.lb.LBFluid
elif espressomd.has_features('LB_GPU') and 'LB.GPU' in modes:
    LB_implementation = espressomd.lb.LBFluidGPU
if LB_implementation:
    lbf = LB_implementation(agrid=0.5, visc=1.3, dens=1.5, tau=0.01, fric=2.0)
    system.actors.add(lbf)

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

system.thermostat.set_langevin(kT=1.0, gamma=2.0)

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
if LB_implementation:
    lbf[1, 1, 1].velocity = [0.1, 0.2, 0.3]
    lbf.save_checkpoint(
        "@CMAKE_CURRENT_BINARY_DIR@/lb_@TEST_COMBINATION@_@TEST_BINARY@.cpt",
        int("@TEST_BINARY@"))
if espressomd.has_features("COLLISION_DETECTION"):
        system.collision_detection.set_params(
            mode="bind_centers", distance=0.11, bond_centers=harmonic_bond)
checkpoint.save(0)
