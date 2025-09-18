/*
 * Copyright (C) 2014-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config/config.hpp"

#include "System.hpp"
#include "System.impl.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "PropagationMode.hpp"
#include "accumulators/AutoUpdateAccumulators.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cell_system/HybridDecomposition.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "electrostatics/icc.hpp"
#include "errorhandling.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "npt.hpp"
#include "particle_node.hpp"
#include "short_range_cabana.hpp"
#include "thermostat.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/all_compare.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

namespace System {

static std::shared_ptr<System> instance = System::create();

std::shared_ptr<System> System::create() {
  auto handle = std::make_shared<System>(Private());
  handle->initialize();
  return handle;
}

System::System(Private) {
  box_geo = std::make_shared<BoxGeometry>();
  local_geo = std::make_shared<LocalBox>();
  cell_structure = std::make_shared<CellStructure>(*box_geo);
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  cell_structure->set_kokkos_handle(::kokkos_handle);
#endif
  propagation = std::make_shared<Propagation>();
  bonded_ias = std::make_shared<BondedInteractionsMap>();
  thermostat = std::make_shared<Thermostat::Thermostat>();
  nonbonded_ias = std::make_shared<InteractionsNonBonded>();
  comfixed = std::make_shared<ComFixed>();
  galilei = std::make_shared<Galilei>();
  oif_global = std::make_shared<OifGlobal>();
  immersed_boundaries = std::make_shared<ImmersedBoundaries>();
#ifdef ESPRESSO_COLLISION_DETECTION
  collision_detection =
      std::make_shared<CollisionDetection::CollisionDetection>();
#endif
  bond_breakage = std::make_shared<BondBreakage::BondBreakage>();
  lees_edwards = std::make_shared<LeesEdwards::LeesEdwards>();
  auto_update_accumulators =
      std::make_shared<Accumulators::AutoUpdateAccumulators>();
  constraints = std::make_shared<Constraints::Constraints>();
#ifdef ESPRESSO_NPT
  nptiso = std::make_shared<NptIsoParameters>();
  npt_inst_pressure = std::make_shared<InstantaneousPressure>();
#endif
  reinit_thermo = true;
  time_step = -1.;
  sim_time = 0.;
  force_cap = 0.;
  min_global_cut = inactive_cutoff;
}

void System::initialize() {
  auto handle = shared_from_this();
  cell_structure->bind_system(handle);
  lees_edwards->bind_system(handle);
  immersed_boundaries->bind_system(handle);
  bonded_ias->bind_system(handle);
  thermostat->bind_system(handle);
  nonbonded_ias->bind_system(handle);
  oif_global->bind_system(handle);
  immersed_boundaries->bind_system(handle);
#ifdef ESPRESSO_COLLISION_DETECTION
  collision_detection->bind_system(handle);
#endif
  auto_update_accumulators->bind_system(handle);
  constraints->bind_system(handle);
#ifdef ESPRESSO_CUDA
  gpu.bind_system(handle);
  gpu.initialize();
#endif
  lb.bind_system(handle);
  ek.bind_system(handle);
}

void reset_system() { instance.reset(); }

void set_system(std::shared_ptr<System> new_instance) {
  instance = new_instance;
}

System &get_system() { return *instance; }

void System::set_time_step(double value) {
  if (value <= 0.)
    throw std::domain_error("time_step must be > 0.");
  if (lb.is_solver_set()) {
    lb.veto_time_step(value);
  }
  if (ek.is_solver_set()) {
    ek.veto_time_step(value);
  }
  time_step = value;
  on_timestep_change();
}

void System::check_kT(double value) const {
  if (lb.is_solver_set()) {
    lb.veto_kT(value);
  }
  if (ek.is_solver_set()) {
    ek.veto_kT(value);
  }
}

void System::set_force_cap(double value) {
  force_cap = value;
  propagation->recalc_forces = true;
}

void System::set_min_global_cut(double value) {
  min_global_cut = value;
  on_verlet_skin_change();
}

void System::set_cell_structure_topology(CellStructureType topology) {
  if (topology == CellStructureType::REGULAR) {
    if (cell_structure->decomposition_type() == CellStructureType::REGULAR) {
      // get fully connected info from exising regular decomposition
      auto &old_regular_decomposition =
          dynamic_cast<RegularDecomposition const &>(
              std::as_const(*cell_structure).decomposition());
      cell_structure->set_regular_decomposition(
          get_interaction_range(),
          old_regular_decomposition.fully_connected_boundary());
    } else { // prev. decomposition is not a regular decomposition
      cell_structure->set_regular_decomposition(get_interaction_range(), {});
    }
  } else if (topology == CellStructureType::NSQUARE) {
    cell_structure->set_atom_decomposition();
  } else {
    assert(topology == CellStructureType::HYBRID);
    /* Get current HybridDecomposition to extract n_square_types */
    auto &old_hybrid_decomposition = dynamic_cast<HybridDecomposition const &>(
        std::as_const(*cell_structure).decomposition());
    cell_structure->set_hybrid_decomposition(
        old_hybrid_decomposition.get_cutoff_regular(),
        old_hybrid_decomposition.get_n_square_types());
  }
}

void System::rebuild_cell_structure() {
  set_cell_structure_topology(cell_structure->decomposition_type());
}

void System::on_boxl_change(bool skip_method_adaption) {
  update_local_geo();
  rebuild_cell_structure();

  /* Now give methods a chance to react to the change in box length */
  if (not skip_method_adaption) {
    lb.on_boxl_change();
    ek.on_boxl_change();
#ifdef ESPRESSO_ELECTROSTATICS
    coulomb.on_boxl_change();
#endif
#ifdef ESPRESSO_DIPOLES
    dipoles.on_boxl_change();
#endif
  }
  constraints->on_boxl_change();
}

void System::veto_boxl_change(bool skip_particle_checks) const {
  if (not skip_particle_checks) {
    auto const n_part = boost::mpi::all_reduce(
        ::comm_cart, cell_structure->local_particles().size(), std::plus<>());
    if (n_part > 0ul) {
      throw std::runtime_error(
          "Cannot reset the box length when particles are present");
    }
  }
  constraints->veto_boxl_change();
  lb.veto_boxl_change();
  ek.veto_boxl_change();
}

void System::on_node_grid_change() {
  update_local_geo();
  lb.on_node_grid_change();
  ek.on_node_grid_change();
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_node_grid_change();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.on_node_grid_change();
#endif
  rebuild_cell_structure();
}

void System::on_periodicity_change() {
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_periodicity_change();
#endif

#ifdef ESPRESSO_DIPOLES
  dipoles.on_periodicity_change();
#endif

#ifdef ESPRESSO_STOKESIAN_DYNAMICS
  if (propagation->integ_switch == INTEG_METHOD_SD) {
    if (box_geo->periodic(0u) or box_geo->periodic(1u) or box_geo->periodic(2u))
      runtimeErrorMsg() << "Stokesian Dynamics requires periodicity "
                        << "(False, False, False)\n";
  }
#endif
  on_verlet_skin_change();
}

void System::on_cell_structure_change() {
  clear_particle_node();
  lb.on_cell_structure_change();
  ek.on_cell_structure_change();
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_cell_structure_change();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.on_cell_structure_change();
#endif
}

void System::on_thermostat_param_change() { reinit_thermo = true; }

void System::on_verlet_skin_change() {
  rebuild_cell_structure();
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_coulomb_change();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void System::on_temperature_change() {
  lb.on_temperature_change();
  ek.on_temperature_change();
}

void System::on_timestep_change() {
  lb.on_timestep_change();
  ek.on_timestep_change();
  on_thermostat_param_change();
}

void System::on_short_range_ia_change() {
  rebuild_cell_structure();
  propagation->recalc_forces = true;
}

void System::on_non_bonded_ia_change() {
  nonbonded_ias->recalc_maximal_cutoffs();
  rebuild_cell_structure();
  on_thermostat_param_change();
  propagation->recalc_forces = true;
}

void System::on_coulomb_change() {
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_coulomb_change();
#endif
  on_short_range_ia_change();
}

void System::on_dipoles_change() {
#ifdef ESPRESSO_DIPOLES
  dipoles.on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void System::on_constraint_change() { propagation->recalc_forces = true; }

void System::on_lb_boundary_conditions_change() {
  propagation->recalc_forces = true;
}

void System::on_particle_local_change() {
  cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
  propagation->recalc_forces = true;
}

void System::on_particle_change() {
  if (cell_structure->decomposition_type() == CellStructureType::HYBRID) {
    cell_structure->set_resort_particles(Cells::RESORT_GLOBAL);
  } else {
    cell_structure->set_resort_particles(Cells::RESORT_LOCAL);
  }
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_particle_change();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.on_particle_change();
#endif
  propagation->recalc_forces = true;

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  cell_structure->clear_local_properties();
#endif
}

void System::on_particle_charge_change() {
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_particle_change();
#endif
}

void System::update_dependent_particles() {
#ifdef ESPRESSO_VIRTUAL_SITES
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  vs_relative_update_particles(*cell_structure, *box_geo);
#endif
  cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
#endif

#ifdef ESPRESSO_ELECTROSTATICS
  if (has_icc_enabled()) {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    rebuild_aosoa();
#endif
    update_icc_particles();
  }
#endif

  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries->init_volume_conservation(*cell_structure);
}

void System::on_observable_calc() {
  /* Prepare particle structure: Communication step: number of ghosts and ghost
   * information */
  cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
  update_dependent_particles();

#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.on_observable_calc();
#endif

#ifdef ESPRESSO_DIPOLES
  dipoles.on_observable_calc();
#endif

  clear_particle_node();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  rebuild_aosoa();
#endif
}

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
void System::rebuild_aosoa() {
#ifdef ESPRESSO_COLLISION_DETECTION
  auto const collision_detection_cutoff = collision_detection->cutoff();
#else
  auto const collision_detection_cutoff = inactive_cutoff;
#endif

  VerletCriterion<> const verlet_criterion{*this,
                                           cell_structure->get_verlet_skin(),
                                           get_interaction_range(),
                                           coulomb.cutoff(),
                                           dipoles.cutoff(),
                                           collision_detection_cutoff};

  update_cabana_state(*cell_structure, verlet_criterion,
                      get_interaction_range(), propagation->integ_switch);
}
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

void System::on_lees_edwards_change() { lb.on_lees_edwards_change(); }

void System::update_local_geo() {
  *local_geo = LocalBox::make_regular_decomposition(
      box_geo->length(), ::communicator.calc_node_index(),
      ::communicator.node_grid);
}

double System::maximal_cutoff() const {
  auto max_cut = inactive_cutoff;
  max_cut = std::max(max_cut, get_min_global_cut());
  max_cut = std::max(max_cut, coulomb.cutoff());
  max_cut = std::max(max_cut, dipoles.cutoff());
  if (::communicator.size > 1) {
    // If there is just one node, the bonded cutoff can be omitted
    // because bond partners are always on the local node.
    max_cut = std::max(max_cut, bonded_ias->maximal_cutoff());
  }
  max_cut = std::max(max_cut, nonbonded_ias->maximal_cutoff());

#ifdef ESPRESSO_COLLISION_DETECTION
  max_cut = std::max(max_cut, collision_detection->cutoff());
#endif
  return max_cut;
}

bool System::long_range_interactions_sanity_checks() const {
  try {
#ifdef ESPRESSO_ELECTROSTATICS
    coulomb.sanity_checks();
#endif
#ifdef ESPRESSO_DIPOLES
    dipoles.sanity_checks();
#endif
  } catch (std::runtime_error const &err) {
    runtimeErrorMsg() << err.what();
    return true;
  }
  return false;
}

double System::get_interaction_range() const {
  auto const max_cut = maximal_cutoff();
  auto const verlet_skin = cell_structure->get_verlet_skin();
  /* Consider skin only if there are actually interactions */
  return (max_cut > 0.) ? max_cut + verlet_skin : inactive_cutoff;
}

void System::set_box_l(Utils::Vector3d const &box_l) {
  box_geo->set_length(box_l);
  on_boxl_change();
}

void System::on_integration_start() {
  // sanity checks
  integrator_sanity_checks();
  long_range_interactions_sanity_checks();
  lb.sanity_checks();
  ek.sanity_checks();

#ifdef ESPRESSO_NPT
  if (propagation->integ_switch == INTEG_METHOD_NPT_ISO_AND ||
      propagation->integ_switch == INTEG_METHOD_NPT_ISO_MTK) {
    npt_ensemble_init(propagation->recalc_forces);
  }
#endif

  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermostat->recalc_prefactors(time_step);
    reinit_thermo = false;
    propagation->recalc_forces = true;
  }

  invalidate_fetch_cache();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  cell_structure->clear_local_properties();
#endif

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  if (!Utils::Mpi::all_compare(::comm_cart, cell_structure->use_verlet_list)) {
    runtimeErrorMsg() << "Nodes disagree about use of verlet lists.";
  }
#ifdef ESPRESSO_ELECTROSTATICS
  {
    auto const &actor = coulomb.impl->solver;
    if (not Utils::Mpi::all_compare(::comm_cart, static_cast<bool>(actor)) or
        (actor and not Utils::Mpi::all_compare(::comm_cart, (*actor).index())))
      runtimeErrorMsg() << "Nodes disagree about Coulomb long-range method";
  }
#endif
#ifdef ESPRESSO_DIPOLES
  {
    auto const &actor = dipoles.impl->solver;
    if (not Utils::Mpi::all_compare(::comm_cart, static_cast<bool>(actor)) or
        (actor and not Utils::Mpi::all_compare(::comm_cart, (*actor).index())))
      runtimeErrorMsg() << "Nodes disagree about dipolar long-range method";
  }
#endif
#endif /* ESPRESSO_ADDITIONAL_CHECKS */

  on_observable_calc();
}

/**
 * @brief Returns the ghost flags required for running pair
 *        kernels for the global state, e.g. the force calculation.
 * @return Required data parts;
 */
unsigned System::get_global_ghost_flags() const {
  /* Position and Properties are always requested. */
  unsigned data_parts = Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES;

  if (lb.is_solver_set())
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (thermostat->thermo_switch & THERMO_DPD)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (thermostat->thermo_switch & THERMO_BOND) {
    data_parts |= Cells::DATA_PART_MOMENTUM;
    data_parts |= Cells::DATA_PART_BONDS;
  }

#ifdef ESPRESSO_COLLISION_DETECTION
  if (not collision_detection->is_off()) {
    data_parts |= Cells::DATA_PART_BONDS;
  }
#endif

  return data_parts;
}

#ifdef ESPRESSO_NPT
bool System::System::has_npt_enabled() const {
  return (propagation->integ_switch == INTEG_METHOD_NPT_ISO_AND) or
         (propagation->integ_switch == INTEG_METHOD_NPT_ISO_MTK);
}
#endif

Utils::Vector3d *System::System::get_npt_virial() const {
#ifdef ESPRESSO_NPT
  if (has_npt_enabled()) {
    return &npt_inst_pressure->p_vir;
  }
#endif
  return nullptr;
}

#ifdef ESPRESSO_COLLISION_DETECTION
bool System::System::has_collision_detection_enabled() const {
  return not collision_detection->is_off();
}
#endif

} // namespace System
