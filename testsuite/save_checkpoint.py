import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.interactions
import espressomd.virtual_sites
import espressomd.accumulators
import espressomd.observables

checkpoint = espressomd.checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.cell_system.skin = 0.4
system.seed = system.cell_system.get_state()["n_nodes"] *[1234]
system.time_step = 0.01
system.min_global_cut = 2.0

system.part.add(pos=[1.0]*3)
system.part.add(pos=[1.0, 1.0, 2.0])
if espressomd.has_features('ELECTROSTATICS'):
    system.part[0].q = 1
    system.part[1].q = -1
    p3m = espressomd.electrostatics.P3M(prefactor=1.0, accuracy=0.1, mesh=10, cao=1, alpha=1.0, r_cut=1.0, tune=False)
    system.actors.add(p3m)
obs = espressomd.observables.ParticlePositions(ids=[0,1])
acc = espressomd.accumulators.MeanVarianceCalculator(obs=obs)
acc.update()
system.part[0].pos = [1.0, 2.0, 3.0]
acc.update()

system.thermostat.set_langevin(kT=1.0, gamma=2.0)

if espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']):
    system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_velocity = True,
                                            have_quaternion = True)
    system.part[1].vs_auto_relate_to(0)
if espressomd.has_features(['LENNARD_JONES']):
    system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=1.2, sigma=1.3, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[3, 0].lennard_jones.set_params(epsilon=1.2, sigma=1.7, cutoff=2.0, shift=0.1)

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
checkpoint.save(0)
