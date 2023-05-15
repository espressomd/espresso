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

import unittest as ut
import unittest_generator as utg
import numpy as np
import pathlib
import tempfile

import espressomd
import espressomd.checkpointing
import espressomd.code_info
import espressomd.electrostatics
import espressomd.magnetostatics
import espressomd.interactions
import espressomd.lees_edwards
import espressomd.drude_helpers
import espressomd.virtual_sites
import espressomd.accumulators
import espressomd.observables
import espressomd.io.writer
import espressomd.lb
import espressomd.electrokinetics
import espressomd.shapes
import espressomd.constraints
import espressomd.bond_breakage
import espressomd.reaction_methods

config = utg.TestGenerator()
modes = config.get_modes()

# use a box with 3 different dimensions, unless DipolarP3M is used
system = espressomd.System(box_l=[12.0, 8.0, 16.0])
if 'DP3M' in modes:
    system.box_l = 3 * [float(np.max(system.box_l))]
system.cell_system.skin = 0.1
system.time_step = 0.01
system.time = 1.5
system.force_cap = 1e8
system.min_global_cut = 2.0
system.max_oif_objects = 5

# create checkpoint folder
config.cleanup_old_checkpoint()
checkpoint = espressomd.checkpointing.Checkpoint(
    **config.get_checkpoint_params())
path_cpt_root = pathlib.Path(checkpoint.checkpoint_dir)

# cleanup old checkpoint files
for filepath in path_cpt_root.iterdir():
    filepath.unlink(missing_ok=True)

# Lees-Edwards boundary conditions
if 'INT.NPT' not in modes:
    protocol = espressomd.lees_edwards.LinearShear(
        initial_pos_offset=0.1, time_0=0.2, shear_velocity=1.2)
    system.lees_edwards.set_boundary_conditions(
        shear_direction="x", shear_plane_normal="y", protocol=protocol)

lbf_class = None
lb_lattice = None
if espressomd.has_features('WALBERLA') and 'LB.WALBERLA' in modes:
    lbf_class = espressomd.lb.LBFluidWalberla
    lb_lattice = espressomd.lb.LatticeWalberla(agrid=2.0, n_ghost_layers=1)
if lbf_class:
    lbf_cpt_mode = 0 if 'LB.ASCII' in modes else 1
    lbf = lbf_class(
        lattice=lb_lattice, kinematic_viscosity=1.3, density=1.5, tau=0.01)
    wall1 = espressomd.shapes.Wall(normal=(1, 0, 0), dist=1.0)
    wall2 = espressomd.shapes.Wall(normal=(-1, 0, 0),
                                   dist=-(system.box_l[0] - 1.0))
    lbf.add_boundary_from_shape(wall1, (1e-4, 1e-4, 0))
    lbf.add_boundary_from_shape(wall2, (0, 0, 0))

    ek_solver = espressomd.electrokinetics.EKNone(lattice=lb_lattice)
    ek_species = espressomd.electrokinetics.EKSpecies(
        lattice=lb_lattice, density=1.5, kT=2.0, diffusion=0.2, valency=0.1,
        advection=False, friction_coupling=False, ext_efield=[0.1, 0.2, 0.3],
        single_precision=False, tau=system.time_step)
    system.ekcontainer.solver = ek_solver
    system.ekcontainer.tau = ek_species.tau
    system.ekcontainer.add(ek_species)
    ek_species.add_boundary_from_shape(
        shape=wall1, value=1e-3 * np.array([1., 2., 3.]),
        boundary_type=espressomd.electrokinetics.FluxBoundary)
    ek_species.add_boundary_from_shape(
        shape=wall2, value=1e-3 * np.array([4., 5., 6.]),
        boundary_type=espressomd.electrokinetics.FluxBoundary)
    ek_species.add_boundary_from_shape(
        shape=wall1, value=1.,
        boundary_type=espressomd.electrokinetics.DensityBoundary)
    ek_species.add_boundary_from_shape(
        shape=wall2, value=2.,
        boundary_type=espressomd.electrokinetics.DensityBoundary)

p1 = system.part.add(id=0, pos=[1.0, 1.0, 1.0])
p2 = system.part.add(id=1, pos=[1.0, 1.0, 2.0])

if espressomd.has_features('ELECTROSTATICS'):
    p1.q = 1
    p2.q = -1

if espressomd.has_features('DIPOLES'):
    p1.dip = (1.3, 2.1, -6)
    p2.dip = (7.3, 6.1, -4)

if espressomd.has_features('EXCLUSIONS'):
    system.part.add(id=2, pos=[2.0, 2.0, 2.0], exclusions=[0, 1])

# place particles at the interface between 2 MPI nodes
p3 = system.part.add(id=3, pos=system.box_l / 2.0 - 1.0, type=1)
p4 = system.part.add(id=4, pos=system.box_l / 2.0 + 1.0, type=1)

system.comfixed.types = [0, 2]
p_slice = system.part.by_ids([4, 1])

if espressomd.has_features('P3M') and ('P3M' in modes or 'ELC' in modes):
    if espressomd.gpu_available() and 'P3M.GPU' in modes:
        ActorP3M = espressomd.electrostatics.P3MGPU
    else:
        ActorP3M = espressomd.electrostatics.P3M
    p3m = ActorP3M(
        prefactor=1.0,
        accuracy=0.1,
        mesh=10,
        cao=1,
        alpha=1.0,
        r_cut=1.0,
        check_complex_residuals=False,
        timings=15,
        tune=False)
    if 'ELC' in modes:
        elc = espressomd.electrostatics.ELC(
            actor=p3m,
            gap_size=6.0,
            maxPWerror=0.1,
            delta_mid_top=0.9,
            delta_mid_bot=0.1)
        system.actors.add(elc)
        elc.charge_neutrality_tolerance = 7e-12
    else:
        system.actors.add(p3m)
        p3m.charge_neutrality_tolerance = 5e-12

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

system.auto_update_accumulators.add(acc_mean_variance)
system.auto_update_accumulators.add(acc_time_series)
system.auto_update_accumulators.add(acc_correlator)

# constraints
system.constraints.add(shape=espressomd.shapes.Sphere(center=system.box_l / 2, radius=0.1),
                       particle_type=7)
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
union = espressomd.shapes.Union()
union.add([espressomd.shapes.Wall(normal=[1., 0., 0.], dist=0.5),
           espressomd.shapes.Wall(normal=[0., 1., 0.], dist=1.5)])
system.constraints.add(shape=union, particle_type=2)
if espressomd.has_features("ELECTROSTATICS"):
    system.constraints.add(espressomd.constraints.ElectricPlaneWave(
        E0=[1., -2., 3.], k=[-.1, .2, .3], omega=5., phi=1.4))

if 'LB' not in modes:
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
        system.periodicity = [False, False, False]
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
        system.periodicity = [False, False, False]
        system.integrator.set_stokesian_dynamics(
            approximation_method='ft', viscosity=0.5, radii={0: 1.5},
            pair_mobility=False, self_mobility=True)

if espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']):
    system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(
        have_quaternion=True)
    system.virtual_sites.have_quaternion = True
    system.virtual_sites.override_cutoff_check = True
    p2.vs_auto_relate_to(p1)

# non-bonded interactions
if espressomd.has_features(['LENNARD_JONES']) and 'LJ' in modes:
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.3, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[3, 0].lennard_jones.set_params(
        epsilon=1.2, sigma=1.7, cutoff=2.0, shift=0.1)
    system.non_bonded_inter[1, 7].lennard_jones.set_params(
        epsilon=1.2e5, sigma=1.7, cutoff=2.0, shift=0.1)
    handle_ia = espressomd.interactions.NonBondedInteractionHandle(
        _types=(0, 0))
    checkpoint.register("handle_ia")
if espressomd.has_features(['DPD']):
    dpd_params = {"weight_function": 1, "gamma": 2., "trans_r_cut": 2., "k": 2.,
                  "trans_weight_function": 0, "trans_gamma": 1., "r_cut": 2.}
    dpd_ia = espressomd.interactions.DPDInteraction(**dpd_params)
    checkpoint.register("dpd_ia")
    checkpoint.register("dpd_params")

# bonded interactions
harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0.0, k=1.0)
system.bonded_inter.add(harmonic_bond)
p2.add_bond((harmonic_bond, p1))
# create 3 thermalized bonds that will overwrite each other's seed
therm_params = dict(temp_com=0.1, temp_distance=0.2, gamma_com=0.3,
                    gamma_distance=0.5, r_cut=2.)
therm_bond1 = espressomd.interactions.ThermalizedBond(seed=1, **therm_params)
therm_bond2 = espressomd.interactions.ThermalizedBond(seed=2, **therm_params)
therm_bond3 = espressomd.interactions.ThermalizedBond(seed=3, **therm_params)
system.bonded_inter.add(therm_bond1)
p2.add_bond((therm_bond1, p1))
checkpoint.register("therm_bond2")
checkpoint.register("therm_params")
# create Drude particles
if espressomd.has_features(['ELECTROSTATICS', 'MASS', 'ROTATION']):
    dh = espressomd.drude_helpers.DrudeHelpers()
    dh.add_drude_particle_to_core(
        system=system, harmonic_bond=harmonic_bond,
        thermalized_bond=therm_bond1, p_core=p2, type_drude=10,
        alpha=1., mass_drude=0.6, coulomb_prefactor=0.8, thole_damping=2.)
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
break_spec = espressomd.bond_breakage.BreakageSpec(
    breakage_length=5., action_type="delete_bond")
system.bond_breakage[strong_harmonic_bond._bond_id] = break_spec

checkpoint.register("system")
checkpoint.register("acc_mean_variance")
checkpoint.register("acc_time_series")
checkpoint.register("acc_correlator")
checkpoint.register("ibm_volcons_bond")
checkpoint.register("ibm_tribend_bond")
checkpoint.register("ibm_triel_bond")
checkpoint.register("break_spec")
checkpoint.register("p_slice")

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
        r_cut=2.4,
        cao=1,
        mesh=[8, 8, 8],
        alpha=12,
        accuracy=0.01,
        timings=15,
        tune=False)
    system.actors.add(dp3m)

if espressomd.has_features('SCAFACOS') and 'SCAFACOS' in modes \
        and 'p3m' in espressomd.code_info.scafacos_methods():
    system.actors.add(espressomd.electrostatics.Scafacos(
        prefactor=0.5,
        method_name="p3m",
        method_params={
            "p3m_r_cut": 1.0,
            "p3m_grid": 64,
            "p3m_cao": 7,
            "p3m_alpha": 2.084652}))

if espressomd.has_features('SCAFACOS_DIPOLES') and 'SCAFACOS' in modes \
        and 'p2nfft' in espressomd.code_info.scafacos_methods():
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

if lbf_class:
    system.actors.add(lbf)
    if 'THERM.LB' in modes:
        system.thermostat.set_lb(LB_fluid=lbf, seed=23, gamma=2.0)
    # Create a 3D grid with deterministic values to fill the LB fluid lattice
    m = np.pi / 12
    grid_3D = np.fromfunction(
        lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
        lbf.shape, dtype=float)
    lbf[:, :, :].population = np.einsum(
        'abc,d->abcd', grid_3D, np.arange(1, 20))
    lbf[:, :, :].last_applied_force = np.einsum(
        'abc,d->abcd', grid_3D, np.arange(1, 4))
    # save LB checkpoint file
    lbf_cpt_path = path_cpt_root / "lb.cpt"
    lbf.save_checkpoint(str(lbf_cpt_path), lbf_cpt_mode)
    # save EK checkpoint file
    ek_species[:, :, :].density = grid_3D
    ek_cpt_path = path_cpt_root / "ek.cpt"
    ek_species.save_checkpoint(str(ek_cpt_path), lbf_cpt_mode)
    # setup VTK folder
    vtk_suffix = config.test_name
    vtk_root = pathlib.Path("vtk_out")
    # create LB VTK callbacks
    lb_vtk_auto_id = f"auto_lb_{vtk_suffix}"
    lb_vtk_manual_id = f"manual_lb_{vtk_suffix}"
    config.recursive_unlink(vtk_root / lb_vtk_auto_id)
    config.recursive_unlink(vtk_root / lb_vtk_manual_id)
    lb_vtk_auto = espressomd.lb.VTKOutput(
        identifier=lb_vtk_auto_id, delta_N=1,
        observables=('density', 'velocity_vector'), base_folder=str(vtk_root))
    lbf.add_vtk_writer(vtk=lb_vtk_auto)
    lb_vtk_auto.disable()
    lb_vtk_manual = espressomd.lb.VTKOutput(
        identifier=lb_vtk_manual_id, delta_N=0,
        observables=('density',), base_folder=str(vtk_root))
    lbf.add_vtk_writer(vtk=lb_vtk_manual)
    lb_vtk_manual.write()
    # create EK VTK callbacks
    ek_vtk_auto_id = f"auto_ek_{vtk_suffix}"
    ek_vtk_manual_id = f"manual_ek_{vtk_suffix}"
    config.recursive_unlink(vtk_root / ek_vtk_auto_id)
    config.recursive_unlink(vtk_root / ek_vtk_manual_id)
    ek_vtk_auto = espressomd.electrokinetics.VTKOutput(
        identifier=ek_vtk_auto_id,
        observables=('density',), delta_N=1, base_folder=str(vtk_root))
    ek_species.add_vtk_writer(vtk=ek_vtk_auto)
    ek_vtk_auto.disable()
    ek_vtk_manual = espressomd.electrokinetics.VTKOutput(
        identifier=ek_vtk_manual_id,
        observables=('density',), delta_N=0, base_folder=str(vtk_root))
    ek_species.add_vtk_writer(vtk=ek_vtk_manual)
    ek_vtk_manual.write()


# set various properties
p8 = system.part.add(id=8, pos=[2.0] * 3 + system.box_l)
p8.lees_edwards_offset = 0.2
p4.v = [-1., 2., -4.]
if espressomd.has_features('MASS'):
    p3.mass = 1.5
if espressomd.has_features('ROTATION'):
    p3.quat = [1., 2., 3., 4.]
    p4.director = [3., 2., 1.]
    p4.omega_body = [0.3, 0.5, 0.7]
    p3.rotation = [True, False, True]
if espressomd.has_features('EXTERNAL_FORCES'):
    p3.fix = [False, True, False]
    p3.ext_force = [-0.6, 0.1, 0.2]
if espressomd.has_features(['EXTERNAL_FORCES', 'ROTATION']):
    p3.ext_torque = [0.3, 0.5, 0.7]
if espressomd.has_features('ROTATIONAL_INERTIA'):
    p3.rinertia = [2., 3., 4.]
if espressomd.has_features('THERMOSTAT_PER_PARTICLE'):
    gamma = 2.
    if espressomd.has_features('PARTICLE_ANISOTROPY'):
        gamma = np.array([2., 3., 4.])
    p4.gamma = gamma
    if espressomd.has_features('ROTATION'):
        p3.gamma_rot = 2. * gamma
if espressomd.has_features('ENGINE'):
    p3.swimming = {"f_swim": 0.03}
if espressomd.has_features('ENGINE') and lbf_class:
    p4.swimming = {"v_swim": 0.02, "mode": "puller", "dipole_length": 1.}
if espressomd.has_features('LB_ELECTROHYDRODYNAMICS') and lbf_class:
    p8.mu_E = [-0.1, 0.2, -0.3]

# h5md output
if espressomd.has_features("H5MD"):
    h5_units = espressomd.io.writer.h5md.UnitSystem(
        time="ps", mass="u", length="m", charge="e")
    h5 = espressomd.io.writer.h5md.H5md(
        file_path=str(path_cpt_root / "test.h5"),
        unit_system=h5_units)
    h5.write()
    h5.flush()
    h5.close()
    checkpoint.register("h5")
    checkpoint.register("h5_units")

# save checkpoint file
checkpoint.save(0)


class TestCheckpoint(ut.TestCase):

    def test_checkpointing(self):
        '''
        Check for the presence of the checkpoint files.
        '''
        self.assertTrue(path_cpt_root.is_dir(),
                        "checkpoint directory not created")

        checkpoint_filepath = path_cpt_root / "0.checkpoint"
        self.assertTrue(checkpoint_filepath.is_file(),
                        "checkpoint file not created")

        if lbf_class:
            self.assertTrue(lbf_cpt_path.is_file(),
                            "LB checkpoint file not created")

        # only objects at global scope can be checkpointed
        with self.assertRaisesRegex(KeyError, "The given object 'local_obj' was not found in the current scope"):
            local_obj = "local"  # pylint: disable=unused-variable
            checkpoint.register("local_obj")

    @ut.skipIf(lbf_class is None, "Skipping test due to missing mode.")
    def test_lb_checkpointing_exceptions(self):
        '''
        Check the LB checkpointing exception mechanism. Write corrupted
        LB checkpoint files that will be tested in ``test_checkpoint.py``.
        '''

        # check exception mechanism
        lbf_cpt_root = lbf_cpt_path.parent
        with self.assertRaisesRegex(RuntimeError, "could not open file"):
            invalid_path = lbf_cpt_root / "unknown_dir" / "lb.cpt"
            lbf.save_checkpoint(str(invalid_path), lbf_cpt_mode)
        with self.assertRaisesRegex(RuntimeError, "unit test error"):
            lbf.save_checkpoint(str(lbf_cpt_root / "lb_err.cpt"), -1)
        with self.assertRaisesRegex(RuntimeError, "could not write to"):
            lbf.save_checkpoint(str(lbf_cpt_root / "lb_err.cpt"), -2)
        with self.assertRaisesRegex(ValueError, "Unknown mode -3"):
            lbf.save_checkpoint(str(lbf_cpt_root / "lb_err.cpt"), -3)
        with self.assertRaisesRegex(ValueError, "Unknown mode 2"):
            lbf.save_checkpoint(str(lbf_cpt_root / "lb_err.cpt"), 2)

        # deactivate LB actor
        system.actors.remove(lbf)

        # read the valid LB checkpoint file
        lbf_cpt_data = lbf_cpt_path.read_bytes()
        cpt_path = str(path_cpt_root / "lb") + "{}.cpt"
        # write checkpoint file with missing data
        with open(cpt_path.format("-missing-data"), "wb") as f:
            f.write(lbf_cpt_data[:len(lbf_cpt_data) // 2])
        # write checkpoint file with extra data
        with open(cpt_path.format("-extra-data"), "wb") as f:
            f.write(lbf_cpt_data + lbf_cpt_data[-8:])
        if lbf_cpt_mode == 0:
            boxsize, popsize, data = lbf_cpt_data.split(b"\n", 2)
            # write checkpoint file with incorrectly formatted data
            with open(cpt_path.format("-wrong-format"), "wb") as f:
                f.write(boxsize + b"\n" + popsize + b"\ntext string\n" + data)
            # write checkpoint file with different box dimensions
            with open(cpt_path.format("-wrong-boxdim"), "wb") as f:
                f.write(b"2" + boxsize + b"\n" + popsize + b"\n" + data)
            # write checkpoint file with different population size
            with open(cpt_path.format("-wrong-popsize"), "wb") as f:
                f.write(boxsize + b"\n" + b"2" + popsize + b"\n" + data)

    @ut.skipIf(lbf_class is None, "Skipping test due to missing mode.")
    def test_ek_checkpointing_exceptions(self):
        '''
        Check the EK checkpointing exception mechanism. Write corrupted
        EK checkpoint files that will be tested in ``test_checkpoint.py``.
        '''

        # check exception mechanism
        ek_cpt_root = ek_cpt_path.parent
        with self.assertRaisesRegex(RuntimeError, "could not open file"):
            invalid_path = ek_cpt_root / "unknown_dir" / "ek.cpt"
            ek_species.save_checkpoint(str(invalid_path), lbf_cpt_mode)
        with self.assertRaisesRegex(RuntimeError, "unit test error"):
            ek_species.save_checkpoint(str(ek_cpt_root / "ek_err.cpt"), -1)
        with self.assertRaisesRegex(RuntimeError, "could not write to"):
            ek_species.save_checkpoint(str(ek_cpt_root / "ek_err.cpt"), -2)
        with self.assertRaisesRegex(ValueError, "Unknown mode -3"):
            ek_species.save_checkpoint(str(ek_cpt_root / "ek_err.cpt"), -3)
        with self.assertRaisesRegex(ValueError, "Unknown mode 2"):
            ek_species.save_checkpoint(str(ek_cpt_root / "ek_err.cpt"), 2)

        # read the valid EK checkpoint file
        ek_cpt_data = ek_cpt_path.read_bytes()
        cpt_path = str(path_cpt_root / "ek") + "{}.cpt"
        # write checkpoint file with missing data
        with open(cpt_path.format("-missing-data"), "wb") as f:
            f.write(ek_cpt_data[:len(ek_cpt_data) // 2])
        # write checkpoint file with extra data
        with open(cpt_path.format("-extra-data"), "wb") as f:
            f.write(ek_cpt_data + ek_cpt_data[-8:])
        if lbf_cpt_mode == 0:
            boxsize, data = ek_cpt_data.split(b"\n", 1)
            # write checkpoint file with incorrectly formatted data
            with open(cpt_path.format("-wrong-format"), "wb") as f:
                f.write(boxsize + b"\n" + b"\ntext string\n" + data)
            # write checkpoint file with different box dimensions
            with open(cpt_path.format("-wrong-boxdim"), "wb") as f:
                f.write(b"2" + boxsize + b"\n" + data)

    def test_generator_recursive_unlink(self):
        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory)
            tree = root / "level1" / "level2"
            tree.mkdir(parents=True, exist_ok=False)
            for dirname in root.iterdir():
                filepath = dirname / "file"
                filepath.write_text("")
            config.recursive_unlink(root)
            for path in root.iterdir():
                self.assertTrue(path.is_dir(),
                                f"Path '{path}' should be a folder")

    def test_reaction_methods_sanity_check(self):
        with self.assertRaisesRegex(RuntimeError, "Reaction methods do not support checkpointing"):
            widom = espressomd.reaction_methods.WidomInsertion(kT=1, seed=1)
            widom._serialize()


if __name__ == '__main__':
    config.bind_test_class(TestCheckpoint)
    ut.main()
