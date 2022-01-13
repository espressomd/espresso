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
import numpy as np
import os

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.magnetostatics
import espressomd.interactions
import espressomd.drude_helpers
import espressomd.virtual_sites
import espressomd.accumulators
import espressomd.observables
import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes
import espressomd.constraints

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

# use a box with 3 different dimensions, unless DipolarP3M is used
system = espressomd.System(box_l=[12.0, 14.0, 16.0])
if 'DP3M' in modes:
    system.box_l = 3 * [np.max(system.box_l)]
system.cell_system.skin = 0.1
system.time_step = 0.01
system.time = 1.5
system.force_cap = 1e8
system.min_global_cut = 2.0
system.max_oif_objects = 5

# create checkpoint folder
idx = "mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace(".", "__")
checkpoint = espressomd.checkpointing.Checkpoint(
    checkpoint_id=idx, checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")

# cleanup old checkpoint files
if checkpoint.has_checkpoints():
    for filepath in os.listdir(checkpoint.checkpoint_dir):
        if filepath.endswith((".checkpoint", ".cpt")):
            os.remove(os.path.join(checkpoint.checkpoint_dir, filepath))

n_nodes = system.cell_system.get_state()["n_nodes"]

LB_implementation = None
if 'LB.CPU' in modes:
    LB_implementation = espressomd.lb.LBFluid
elif 'LB.GPU' in modes and espressomd.gpu_available():
    LB_implementation = espressomd.lb.LBFluidGPU
if LB_implementation:
    lbf = LB_implementation(agrid=0.5, visc=1.3, dens=1.5, tau=0.01,
                            gamma_odd=0.2, gamma_even=0.3)
    system.actors.add(lbf)
    if 'THERM.LB' in modes:
        system.thermostat.set_lb(LB_fluid=lbf, seed=23, gamma=2.0)
    if n_nodes == 1 and (espressomd.has_features(
            "LB_BOUNDARIES") or espressomd.has_features("LB_BOUNDARIES_GPU")):
        system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
            shape=espressomd.shapes.Wall(normal=(1, 0, 0), dist=0.5), velocity=(1e-4, 1e-4, 0)))
        system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
            shape=espressomd.shapes.Wall(normal=(-1, 0, 0), dist=-(system.box_l[0] - 0.5)), velocity=(0, 0, 0)))

p1 = system.part.add(id=0, pos=[1.0] * 3)
p2 = system.part.add(id=1, pos=[1.0, 1.0, 2.0])

if espressomd.has_features('ELECTROSTATICS'):
    p1.q = 1
    p2.q = -1

if espressomd.has_features('DIPOLES'):
    p1.dip = (1.3, 2.1, -6)
    p2.dip = (7.3, 6.1, -4)

if espressomd.has_features('EXCLUSIONS'):
    system.part.add(id=2, pos=[2.0] * 3, exclusions=[0, 1])

# place particles at the interface between 2 MPI nodes
p3 = system.part.add(id=3, pos=system.box_l / 2.0 - 1.0, type=1)
p4 = system.part.add(id=4, pos=system.box_l / 2.0 + 1.0, type=1)

if espressomd.has_features('P3M') and 'P3M' in modes:
    p3m = espressomd.electrostatics.P3M(
        prefactor=1.0,
        accuracy=0.1,
        mesh=10,
        cao=1,
        alpha=1.0,
        r_cut=1.0,
        timings=15,
        tune=False)
    if 'P3M.CPU' in modes:
        system.actors.add(p3m)
    elif 'P3M.ELC' in modes:
        elc = espressomd.electrostatics.ELC(
            p3m_actor=p3m,
            gap_size=6.0,
            maxPWerror=0.1,
            delta_mid_top=0.9,
            delta_mid_bot=0.1)
        system.actors.add(elc)

# accumulators
obs = espressomd.observables.ParticlePositions(ids=[0, 1])
acc_mean_variance = espressomd.accumulators.MeanVarianceCalculator(obs=obs)
acc_time_series = espressomd.accumulators.TimeSeries(obs=obs)
acc_correlator = espressomd.accumulators.Correlator(
    obs1=obs, tau_lin=10, tau_max=2, delta_N=1,
    corr_operation="componentwise_product")
acc_mean_variance.update()
acc_time_series.update()
acc_correlator.update()
p1.pos = [1.0, 2.0, 3.0]
acc_mean_variance.update()
acc_time_series.update()
acc_correlator.update()

if n_nodes == 1:
    system.auto_update_accumulators.add(acc_mean_variance)
    system.auto_update_accumulators.add(acc_time_series)
    system.auto_update_accumulators.add(acc_correlator)

# constraints
if n_nodes == 1:
    system.constraints.add(shape=espressomd.shapes.Sphere(center=system.box_l / 2, radius=0.1),
                           particle_type=17)
    system.constraints.add(
        shape=espressomd.shapes.Wall(
            normal=[1. / np.sqrt(3)] * 3, dist=0.5))
    system.constraints.add(espressomd.constraints.Gravity(g=[1., 2., 3.]))
    system.constraints.add(
        espressomd.constraints.HomogeneousMagneticField(H=[1., 2., 3.]))
    system.constraints.add(
        espressomd.constraints.HomogeneousFlowField(u=[1., 2., 3.], gamma=2.3))
    pot_field_data = espressomd.constraints.ElectricPotential.field_from_fn(
        system.box_l, np.ones(3), lambda x: np.linalg.norm(10 * np.ones(3) - x))
    checkpoint.register("pot_field_data")
    system.constraints.add(espressomd.constraints.PotentialField(
        field=pot_field_data, grid_spacing=np.ones(3), default_scale=1.6,
        particle_scales={5: 6.0}))
    vec_field_data = espressomd.constraints.ForceField.field_from_fn(
        system.box_l, np.ones(3), lambda x: 10 * np.ones(3) - x)
    checkpoint.register("vec_field_data")
    system.constraints.add(espressomd.constraints.ForceField(
        field=vec_field_data, grid_spacing=np.ones(3), default_scale=1.4))
    if espressomd.has_features("ELECTROSTATICS"):
        system.constraints.add(espressomd.constraints.ElectricPlaneWave(
            E0=[1., -2., 3.], k=[-.1, .2, .3], omega=5., phi=1.4))

if 'LB.OFF' in modes:
    # set thermostat
    if 'THERM.LANGEVIN' in modes:
        system.thermostat.set_langevin(kT=1.0, gamma=2.0, seed=42)
    elif 'THERM.BD' in modes:
        system.thermostat.set_brownian(kT=1.0, gamma=2.0, seed=42)
    elif 'THERM.NPT' in modes and espressomd.has_features('NPT'):
        system.thermostat.set_npt(kT=1.0, gamma0=2.0, gammav=0.1, seed=42)
    elif 'THERM.DPD' in modes and espressomd.has_features('DPD'):
        system.thermostat.set_dpd(kT=1.0, seed=42)
    elif 'THERM.SDM' in modes and espressomd.has_features('STOKESIAN_DYNAMICS'):
        system.periodicity = [0, 0, 0]
        system.thermostat.set_stokesian(kT=1.0, seed=42)
    # set integrator
    if 'INT.NPT' in modes and espressomd.has_features('NPT'):
        system.integrator.set_isotropic_npt(ext_pressure=2.0, piston=0.01,
                                            direction=[True, False, False])
    elif 'INT.SD' in modes:
        system.integrator.set_steepest_descent(f_max=2.0, gamma=0.1,
                                               max_displacement=0.01)
    elif 'INT.NVT' in modes:
        system.integrator.set_nvt()
    elif 'INT.BD' in modes:
        system.integrator.set_brownian_dynamics()
    elif 'INT.SDM' in modes and espressomd.has_features('STOKESIAN_DYNAMICS'):
        system.periodicity = [0, 0, 0]
        system.integrator.set_stokesian_dynamics(
            approximation_method='ft', viscosity=0.5, radii={0: 1.5},
            pair_mobility=False, self_mobility=True)

if espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']):
    system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
        have_quaternion=True)
    p2.vs_auto_relate_to(p1)

# non-bonded interactions
if espressomd.has_features(['LENNARD_JONES']) and 'LJ' in modes:
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.3, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[3, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.7, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[1, 17].lennard_jones.set_params(
        epsilon=1.2e5, sigma=1.7, cutoff=2.0, shift=0.1)

# bonded interactions
harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0.0, k=1.0)
system.bonded_inter.add(harmonic_bond)
p2.add_bond((harmonic_bond, p1))
if 'THERM.LB' not in modes:
    thermalized_bond = espressomd.interactions.ThermalizedBond(
        temp_com=0.0, gamma_com=0.0, temp_distance=0.2, gamma_distance=0.5,
        r_cut=2, seed=51)
    system.bonded_inter.add(thermalized_bond)
    p2.add_bond((thermalized_bond, p1))
    if espressomd.has_features(['ELECTROSTATICS', 'MASS', 'ROTATION']):
        dh = espressomd.drude_helpers.DrudeHelpers()
        dh.add_drude_particle_to_core(system, harmonic_bond, thermalized_bond,
                                      p2, 10, 1., 4.6, 0.8, 2.)
        checkpoint.register("dh")
strong_harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0.0, k=5e5)
system.bonded_inter.add(strong_harmonic_bond)
p4.add_bond((strong_harmonic_bond, p3))
ibm_volcons_bond = espressomd.interactions.IBM_VolCons(softID=15, kappaV=0.01)
ibm_tribend_bond = espressomd.interactions.IBM_Tribend(
    ind1=p1.id, ind2=p2.id, ind3=p3.id, ind4=p4.id, kb=2., refShape="Initial")
ibm_triel_bond = espressomd.interactions.IBM_Triel(
    ind1=p1.id, ind2=p2.id, ind3=p3.id, k1=1.1, k2=1.2, maxDist=1.6,
    elasticLaw="NeoHookean")

checkpoint.register("system")
checkpoint.register("acc_mean_variance")
checkpoint.register("acc_time_series")
checkpoint.register("acc_correlator")
checkpoint.register("ibm_volcons_bond")
checkpoint.register("ibm_tribend_bond")
checkpoint.register("ibm_triel_bond")

# calculate forces
system.integrator.run(0)
particle_force0 = np.copy(p1.f)
particle_force1 = np.copy(p2.f)
checkpoint.register("particle_force0")
checkpoint.register("particle_force1")
if espressomd.has_features("COLLISION_DETECTION"):
    system.collision_detection.set_params(
        mode="bind_centers", distance=0.11, bond_centers=harmonic_bond)

if espressomd.has_features('DP3M') and 'DP3M' in modes:
    dp3m = espressomd.magnetostatics.DipolarP3M(
        prefactor=1.,
        epsilon=2.,
        mesh_off=[0.5, 0.5, 0.5],
        r_cut=2.4,
        cao=1,
        mesh=[8, 8, 8],
        alpha=12,
        accuracy=0.01,
        timings=15,
        tune=False)
    system.actors.add(dp3m)

if espressomd.has_features('SCAFACOS') and 'SCAFACOS' in modes \
        and 'p3m' in espressomd.scafacos.available_methods():
    system.actors.add(espressomd.electrostatics.Scafacos(
        prefactor=0.5,
        method_name="p3m",
        method_params={
            "p3m_r_cut": 1.0,
            "p3m_grid": 64,
            "p3m_cao": 7,
            "p3m_alpha": 2.084652}))

if espressomd.has_features('SCAFACOS_DIPOLES') and 'SCAFACOS' in modes \
        and 'p2nfft' in espressomd.scafacos.available_methods():
    system.actors.add(espressomd.magnetostatics.Scafacos(
        prefactor=1.2,
        method_name='p2nfft',
        method_params={
            "p2nfft_verbose_tuning": "0",
            "pnfft_N": "32,32,32",
            "pnfft_n": "32,32,32",
            "pnfft_window_name": "bspline",
            "pnfft_m": "4",
            "p2nfft_ignore_tolerance": "1",
            "pnfft_diff_ik": "0",
            "p2nfft_r_cut": "11",
            "p2nfft_alpha": "0.37"}))

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

# save checkpoint file
checkpoint.save(0)


class TestCheckpointLB(ut.TestCase):

    def test_checkpointing(self):
        '''
        Check for the presence of the checkpoint files.
        '''
        self.assertTrue(os.path.isdir(checkpoint.checkpoint_dir),
                        "checkpoint directory not created")

        checkpoint_filepath = checkpoint.checkpoint_dir + "/0.checkpoint"
        self.assertTrue(os.path.isfile(checkpoint_filepath),
                        "checkpoint file not created")

        if LB_implementation:
            self.assertTrue(os.path.isfile(lbf_cpt_path),
                            "LB checkpoint file not created")
            self.check_lb_checkpointing()

    def check_lb_checkpointing(self):
        '''
        Check the LB checkpointing exception mechanism. Write corrupted
        LB checkpoint files that will be tested in ``test_checkpoint.py``.
        '''

        # check exception mechanism
        with self.assertRaisesRegex(RuntimeError, 'could not open file'):
            dirname, filename = os.path.split(lbf_cpt_path)
            invalid_path = os.path.join(dirname, 'unknown_dir', filename)
            lbf.save_checkpoint(invalid_path, cpt_mode)
        system.actors.remove(lbf)
        with self.assertRaisesRegex(RuntimeError, 'one needs to have already initialized the LB fluid'):
            lbf.load_checkpoint(lbf_cpt_path, cpt_mode)

        # read the valid LB checkpoint file
        with open(lbf_cpt_path, "rb") as f:
            lbf_cpt_str = f.read()
        cpt_path = checkpoint.checkpoint_dir + "/lb{}.cpt"
        # write checkpoint file with missing data
        with open(cpt_path.format("-missing-data"), "wb") as f:
            f.write(lbf_cpt_str[:len(lbf_cpt_str) // 2])
        # write checkpoint file with extra data
        with open(cpt_path.format("-extra-data"), "wb") as f:
            f.write(lbf_cpt_str + lbf_cpt_str[-8:])
        if cpt_mode == 0:
            boxsize, data = lbf_cpt_str.split(b"\n", 1)
            # write checkpoint file with incorrectly formatted data
            with open(cpt_path.format("-wrong-format"), "wb") as f:
                f.write(boxsize + b"\ntext string\n" + data)
            # write checkpoint file with different box dimensions
            with open(cpt_path.format("-wrong-boxdim"), "wb") as f:
                f.write(b"2" + boxsize + b"\n" + data)


if __name__ == '__main__':
    ut.main()
