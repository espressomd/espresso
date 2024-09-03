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
#include "npt.hpp"
#include "particle_node.hpp"
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
  propagation = std::make_shared<Propagation>();
  bonded_ias = std::make_shared<BondedInteractionsMap>();
  thermostat = std::make_shared<Thermostat::Thermostat>();
  nonbonded_ias = std::make_shared<InteractionsNonBonded>();
  comfixed = std::make_shared<ComFixed>();
  galilei = std::make_shared<Galilei>();
  oif_global = std::make_shared<OifGlobal>();
  immersed_boundaries = std::make_shared<ImmersedBoundaries>();
#ifdef COLLISION_DETECTION
  collision_detection =
      std::make_shared<CollisionDetection::CollisionDetection>();
#endif
  bond_breakage = std::make_shared<BondBreakage::BondBreakage>();
  lees_edwards = std::make_shared<LeesEdwards::LeesEdwards>();
  auto_update_accumulators =
      std::make_shared<Accumulators::AutoUpdateAccumulators>();
  constraints = std::make_shared<Constraints::Constraints>();
  reinit_thermo = true;
  time_step = -1.;
  sim_time = 0.;
  force_cap = 0.;
  min_global_cut = INACTIVE_CUTOFF;
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
#ifdef COLLISION_DETECTION
  collision_detection->bind_system(handle);
#endif
  auto_update_accumulators->bind_system(handle);
  constraints->bind_system(handle);
#ifdef CUDA
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
#ifdef ELECTROSTATICS
    coulomb.on_boxl_change();
#endif
#ifdef DIPOLES
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
#ifdef ELECTROSTATICS
  coulomb.on_node_grid_change();
#endif
#ifdef DIPOLES
  dipoles.on_node_grid_change();
#endif
  rebuild_cell_structure();
}

void System::on_periodicity_change() {
#ifdef ELECTROSTATICS
  coulomb.on_periodicity_change();
#endif

#ifdef DIPOLES
  dipoles.on_periodicity_change();
#endif

#ifdef STOKESIAN_DYNAMICS
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
#ifdef ELECTROSTATICS
  coulomb.on_cell_structure_change();
#endif
#ifdef DIPOLES
  dipoles.on_cell_structure_change();
#endif
}

void System::on_thermostat_param_change() { reinit_thermo = true; }

void System::on_verlet_skin_change() {
  rebuild_cell_structure();
#ifdef ELECTROSTATICS
  coulomb.on_coulomb_change();
#endif
#ifdef DIPOLES
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
  propagation->recalc_forces = true;
}

void System::on_coulomb_change() {
#ifdef ELECTROSTATICS
  coulomb.on_coulomb_change();
#endif
  on_short_range_ia_change();
}

void System::on_dipoles_change() {
#ifdef DIPOLES
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
#ifdef ELECTROSTATICS
  coulomb.on_particle_change();
#endif
#ifdef DIPOLES
  dipoles.on_particle_change();
#endif
  propagation->recalc_forces = true;

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void System::on_particle_charge_change() {
#ifdef ELECTROSTATICS
  coulomb.on_particle_change();
#endif
}

void System::update_dependent_particles() {
#ifdef VIRTUAL_SITES
#ifdef VIRTUAL_SITES_RELATIVE
  vs_relative_update_particles(*cell_structure, *box_geo);
#endif
  cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
#endif

#ifdef ELECTROSTATICS
  update_icc_particles();
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

#ifdef ELECTROSTATICS
  coulomb.on_observable_calc();
#endif

#ifdef DIPOLES
  dipoles.on_observable_calc();
#endif

  clear_particle_node();
}

void System::update_local_geo() {
  *local_geo = LocalBox::make_regular_decomposition(
      box_geo->length(), ::communicator.calc_node_index(),
      ::communicator.node_grid);
}

double System::maximal_cutoff() const {
  auto max_cut = INACTIVE_CUTOFF;
  max_cut = std::max(max_cut, get_min_global_cut());
#ifdef ELECTROSTATICS
  max_cut = std::max(max_cut, coulomb.cutoff());
#endif
#ifdef DIPOLES
  max_cut = std::max(max_cut, dipoles.cutoff());
#endif
  if (::communicator.size > 1) {
    // If there is just one node, the bonded cutoff can be omitted
    // because bond partners are always on the local node.
    max_cut = std::max(max_cut, bonded_ias->maximal_cutoff());
  }
  max_cut = std::max(max_cut, nonbonded_ias->maximal_cutoff());

#ifdef COLLISION_DETECTION
  max_cut = std::max(max_cut, collision_detection->cutoff());
#endif
  return max_cut;
}

bool System::long_range_interactions_sanity_checks() const {
  try {
#ifdef ELECTROSTATICS
    coulomb.sanity_checks();
#endif
#ifdef DIPOLES
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
  return (max_cut > 0.) ? max_cut + verlet_skin : INACTIVE_CUTOFF;
}

void System::set_box_l(Utils::Vector3d const &box_l) {
  box_geo->set_length(box_l);
  on_boxl_change();
}

void System::on_integration_start() {
  // sanity checks
  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks();
#endif
  long_range_interactions_sanity_checks();
  lb.sanity_checks();
  ek.sanity_checks();

  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermostat->recalc_prefactors(time_step);
    reinit_thermo = false;
    propagation->recalc_forces = true;
  }

#ifdef NPT
  if (propagation->integ_switch == INTEG_METHOD_NPT_ISO) {
    npt_ensemble_init(box_geo->length(), propagation->recalc_forces);
  }
#endif

  invalidate_fetch_cache();

#ifdef ADDITIONAL_CHECKS
  if (!Utils::Mpi::all_compare(::comm_cart, cell_structure->use_verlet_list)) {
    runtimeErrorMsg() << "Nodes disagree about use of verlet lists.";
  }
#ifdef ELECTROSTATICS
  {
    auto const &actor = coulomb.impl->solver;
    if (not Utils::Mpi::all_compare(::comm_cart, static_cast<bool>(actor)) or
        (actor and not Utils::Mpi::all_compare(::comm_cart, (*actor).index())))
      runtimeErrorMsg() << "Nodes disagree about Coulomb long-range method";
  }
#endif
#ifdef DIPOLES
  {
    auto const &actor = dipoles.impl->solver;
    if (not Utils::Mpi::all_compare(::comm_cart, static_cast<bool>(actor)) or
        (actor and not Utils::Mpi::all_compare(::comm_cart, (*actor).index())))
      runtimeErrorMsg() << "Nodes disagree about dipolar long-range method";
  }
#endif
#endif /* ADDITIONAL_CHECKS */

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

#ifdef COLLISION_DETECTION
  if (not collision_detection->is_off()) {
    data_parts |= Cells::DATA_PART_BONDS;
  }
#endif

  return data_parts;
}

} // namespace System
