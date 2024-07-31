/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "reaction_methods/ReactionAlgorithm.hpp"

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/serialization/serialization.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace boost::serialization {
template <typename Archive, typename... Types>
void serialize(Archive &ar, std::tuple<Types...> &tuple, const unsigned int) {
  std::apply([&](auto &...item) { ((ar & item), ...); }, tuple);
}
} // namespace boost::serialization

namespace ReactionMethods {

/**
 * Adds a reaction to the reaction system
 */
void ReactionAlgorithm::add_reaction(
    std::shared_ptr<SingleReaction> const &new_reaction) {

  // make ESPResSo count the particle numbers which take part in the reactions
  for (int reactant_type : new_reaction->reactant_types)
    init_type_map(reactant_type);
  for (int product_type : new_reaction->product_types)
    init_type_map(product_type);

  init_type_map(non_interacting_type);

  reactions.push_back(new_reaction);
}

/**
 * @details This method tries to keep the cell system overhead to a minimum.
 * Event callbacks are only called once after all particles are updated,
 * except for particle deletion (the cell structure is still reinitialized
 * after each deletion).
 */
void ReactionAlgorithm::restore_old_system_state() {
  auto &system = System::get_system();
  auto const &old_state = get_old_system_state();
  auto const &box_geo = *system.box_geo;
  // restore the properties of changed and hidden particles
  for (auto const &state : {old_state.changed, old_state.hidden}) {
    for (auto const &[p_id, p_type] : state) {
      on_particle_type_change(p_id, type_tracking::any_type, p_type);
      if (auto p = get_local_particle(p_id)) {
        p->type() = p_type;
#ifdef ELECTROSTATICS
        p->q() = charges_of_types.at(p_type);
#endif
      }
    }
  }
  // delete created particles
  for (auto const p_id : old_state.created) {
    delete_particle(p_id);
  }
  // restore original positions and velocities
  for (auto const &[p_id, pos, vel] : old_state.moved) {
    if (auto p = get_local_particle(p_id)) {
      p->v() = vel;
      auto folded_pos = pos;
      auto image_box = Utils::Vector3i{};
      box_geo.fold_position(folded_pos, image_box);
      p->pos() = folded_pos;
      p->image_box() = image_box;
    }
  }
  if (not old_state.moved.empty()) {
    auto const &system = System::get_system();
    auto &cell_structure = *system.cell_structure;
    cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  }
  system.on_particle_change();
  clear_old_system_state();
}

/**
 * Automatically sets the volume which is used by the reaction ensemble to the
 * volume of a cuboid box.
 */
void ReactionAlgorithm::update_volume() {
  auto const &box_geo = *System::get_system().box_geo;
  volume = box_geo.volume();
}

/**
 * Checks whether all particles exist for the provided reaction.
 */
bool ReactionAlgorithm::all_reactant_particles_exist(
    SingleReaction const &current_reaction) const {
  bool enough_particles = true;
  for (int i = 0; i < current_reaction.reactant_types.size(); i++) {
    int current_number =
        number_of_particles_with_type(current_reaction.reactant_types[i]);
    if (current_number < current_reaction.reactant_coefficients[i]) {
      enough_particles = false;
      break;
    }
  }
  return enough_particles;
}

void ReactionAlgorithm::make_reaction_attempt(SingleReaction const &reaction,
                                              ParticleChanges &bookkeeping) {
  // create or hide particles of types with corresponding types in reaction
  auto const n_product_types = reaction.product_types.size();
  auto const n_reactant_types = reaction.reactant_types.size();
  auto const get_random_p_id_of_type = [this](int type) {
    auto const random_index = i_random(number_of_particles_with_type(type));
    return get_random_p_id(type, random_index);
  };
  auto only_local_changes = true;
  for (int i = 0; i < std::min(n_product_types, n_reactant_types); i++) {
    auto const n_product_coef = reaction.product_coefficients[i];
    auto const n_reactant_coef = reaction.reactant_coefficients[i];
    // change std::min(reactant_coefficients(i),product_coefficients(i)) many
    // particles of reactant_types(i) to product_types(i)
    auto const old_type = reaction.reactant_types[i];
    auto const new_type = reaction.product_types[i];
#ifdef ELECTROSTATICS
    if (charges_of_types.at(new_type) != charges_of_types.at(old_type)) {
      only_local_changes = false;
    }
#endif
    for (int j = 0; j < std::min(n_product_coef, n_reactant_coef); j++) {
      auto const p_id = get_random_p_id_of_type(old_type);
      on_particle_type_change(p_id, old_type, new_type);
      if (auto p = get_local_particle(p_id)) {
        p->type() = new_type;
#ifdef ELECTROSTATICS
        p->q() = charges_of_types.at(new_type);
#endif
      }
      bookkeeping.changed.emplace_back(p_id, old_type);
    }
    // create product_coefficients(i)-reactant_coefficients(i) many product
    // particles iff product_coefficients(i)-reactant_coefficients(i)>0,
    // iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number
    // of reactant particles
    auto const delta_n = n_product_coef - n_reactant_coef;
    if (delta_n > 0) {
      auto const type = reaction.product_types[i];
      for (int j = 0; j < delta_n; j++) {
        auto const p_id = create_particle(type);
        check_exclusion_range(p_id, type);
        bookkeeping.created.emplace_back(p_id);
      }
      only_local_changes = false;
    } else if (delta_n < 0) {
      auto const type = reaction.reactant_types[i];
      for (int j = 0; j < -delta_n; j++) {
        auto const p_id = get_random_p_id_of_type(type);
        bookkeeping.hidden.emplace_back(p_id, type);
        check_exclusion_range(p_id, type);
        hide_particle(p_id, type);
      }
      only_local_changes = false;
    }
  }
  // create or hide particles of types with noncorresponding replacement types
  for (auto i = std::min(n_product_types, n_reactant_types);
       i < std::max(n_product_types, n_reactant_types); i++) {
    if (n_product_types < n_reactant_types) {
      // hide superfluous reactant_types particles
      auto const type = reaction.reactant_types[i];
      for (int j = 0; j < reaction.reactant_coefficients[i]; j++) {
        auto const p_id = get_random_p_id_of_type(type);
        bookkeeping.hidden.emplace_back(p_id, type);
        check_exclusion_range(p_id, type);
        hide_particle(p_id, type);
      }
    } else {
      // create additional product_types particles
      auto const type = reaction.product_types[i];
      for (int j = 0; j < reaction.product_coefficients[i]; j++) {
        auto const p_id = create_particle(type);
        check_exclusion_range(p_id, type);
        bookkeeping.created.emplace_back(p_id);
      }
    }
  }
  // determine which fine-grained event to trigger
  if (n_product_types == n_reactant_types and only_local_changes) {
    System::get_system().on_particle_local_change();
  } else {
    System::get_system().on_particle_change();
  }
}

std::unordered_map<int, int>
ReactionAlgorithm::get_particle_numbers(SingleReaction const &reaction) const {
  std::unordered_map<int, int> particle_numbers;
  // reactants
  for (int type : reaction.reactant_types) {
    particle_numbers[type] = ::number_of_particles_with_type(type);
  }
  // products
  for (int type : reaction.product_types) {
    particle_numbers[type] = ::number_of_particles_with_type(type);
  }
  return particle_numbers;
}

std::optional<double>
ReactionAlgorithm::create_new_trial_state(int reaction_id) {
  auto &reaction = *reactions[reaction_id];
  reaction.tried_moves++;
  particle_inside_exclusion_range_touched = false;
  if (!all_reactant_particles_exist(reaction)) {
    // make sure that no incomplete reaction is performed -> only need to
    // consider rollback of complete reactions
    return {};
  }
  auto &bookkeeping = make_new_system_state();
  bookkeeping.reaction_id = reaction_id;
  bookkeeping.old_particle_numbers = get_particle_numbers(reaction);
  make_reaction_attempt(reaction, bookkeeping);
  auto E_pot_new = std::numeric_limits<double>::max();
  if (not particle_inside_exclusion_range_touched) {
    E_pot_new = calculate_potential_energy();
  }
  return {E_pot_new};
}

double ReactionAlgorithm::make_reaction_mc_move_attempt(int reaction_id,
                                                        double bf,
                                                        double E_pot_old,
                                                        double E_pot_new) {
  auto const exponential = std::exp(-(E_pot_new - E_pot_old) / kT);
  auto &reaction = *reactions[reaction_id];
  reaction.accumulator_potential_energy_difference_exponential(
      std::vector<double>{exponential});
  if (get_random_uniform_number() >= bf) {
    // reject trial move: restore previous state, energy is unchanged
    restore_old_system_state();
    return E_pot_old;
  }
  // accept trial move: delete hidden particles and return new system energy
  for (auto const &[p_id, p_type] : get_old_system_state().hidden) {
    delete_particle(p_id);
  }
  reaction.accepted_moves++;
  clear_old_system_state();
  return E_pot_new;
}

/**
 * Hides a particle from short ranged interactions and from the electrostatic
 * interaction. Additional hiding from interactions would need to be implemented
 * here.
 *
 * Removing the particle charge and changing its type to a non existing one
 * deactivates all interactions with other particles, as if the particle was
 * inexistent (currently only type-based interactions are switched off, as well
 * as the electrostatic interaction).
 * This function does not break bonds for simple reactions, as long as there
 * are no reactions like 2A -->B where one of the reacting A particles occurs
 * in the polymer (think of bond breakages if the monomer in the polymer gets
 * deleted in the reaction). This constraint is not of fundamental reason, but
 * there would be a need for a rule for such "collision" reactions (a reaction
 * like the one above).
 */
void ReactionAlgorithm::hide_particle(int p_id, int p_type) const {
  on_particle_type_change(p_id, p_type, non_interacting_type);
  if (auto p = get_local_particle(p_id)) {
    p->type() = non_interacting_type;
#ifdef ELECTROSTATICS
    p->q() = 0.;
#endif
  }
}

/**
 * Check if the inserted particle is too close to neighboring particles.
 */
void ReactionAlgorithm::check_exclusion_range(int p_id, int p_type) {

  /* Check the exclusion radius of the inserted particle */
  if (exclusion_radius_per_type.count(p_type) != 0) {
    if (exclusion_radius_per_type[p_type] == 0.) {
      return;
    }
  }

  auto p1_ptr = get_real_particle(p_id);

  std::vector<int> particle_ids;
  if (neighbor_search_order_n) {
    auto all_ids = get_particle_ids_parallel();
    /* remove the inserted particle id */
    std::erase(all_ids, p_id);
    particle_ids = all_ids;
  } else {
    auto &system = System::get_system();
    system.on_observable_calc();
    auto const local_ids =
        get_short_range_neighbors(system, p_id, m_max_exclusion_range);
    assert(p1_ptr == nullptr or !!local_ids);
    if (local_ids) {
      particle_ids = std::move(*local_ids);
    }
  }

  if (p1_ptr != nullptr) {
    auto &p1 = *p1_ptr;
    auto const &system = System::get_system();
    auto const &box_geo = *system.box_geo;
    auto &cell_structure = *system.cell_structure;

    /* Check if the inserted particle within the exclusion radius of any other
     * particle */
    for (auto const p2_id : particle_ids) {
      if (auto const p2_ptr = cell_structure.get_local_particle(p2_id)) {
        auto const &p2 = *p2_ptr;
        double excluded_distance;
        if (exclusion_radius_per_type.count(p_type) == 0 ||
            exclusion_radius_per_type.count(p2.type()) == 0) {
          excluded_distance = exclusion_range;
        } else if (exclusion_radius_per_type[p2.type()] == 0.) {
          continue;
        } else {
          excluded_distance = exclusion_radius_per_type[p_type] +
                              exclusion_radius_per_type[p2.type()];
        }

        auto const d_min = box_geo.get_mi_vector(p2.pos(), p1.pos()).norm();

        if (d_min < excluded_distance) {
          particle_inside_exclusion_range_touched = true;
          break;
        }
      }
    }
    if (m_comm.rank() != 0) {
      m_comm.send(0, 1, particle_inside_exclusion_range_touched);
    }
  } else if (m_comm.rank() == 0) {
    m_comm.recv(boost::mpi::any_source, 1,
                particle_inside_exclusion_range_touched);
  }
  boost::mpi::broadcast(m_comm, particle_inside_exclusion_range_touched, 0);
}

/**
 * Deletes the particle with the given p_id and stores the id if the deletion
 * created a hole in the particle id range. This method is intended to only
 * delete unbonded particles since bonds are coupled to ids. This is used to
 * avoid the id range becoming excessively huge.
 */
void ReactionAlgorithm::delete_particle(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
  auto const old_max_seen_id = ::get_maximal_particle_id();
  if (p_id == old_max_seen_id) {
    // last particle, just delete
    remove_particle(p_id);
    // remove all saved empty p_ids which are greater than the max_seen_particle
    // this is needed in order to avoid the creation of holes
    for (auto p_id_iter = m_empty_p_ids_smaller_than_max_seen_particle.begin();
         p_id_iter != m_empty_p_ids_smaller_than_max_seen_particle.end();) {
      if ((*p_id_iter) >= old_max_seen_id)
        p_id_iter = m_empty_p_ids_smaller_than_max_seen_particle.erase(
            p_id_iter); // update iterator after container was modified
      else
        ++p_id_iter;
    }
  } else if (p_id <= old_max_seen_id) {
    remove_particle(p_id);
    m_empty_p_ids_smaller_than_max_seen_particle.push_back(p_id);
  } else {
    throw std::runtime_error(
        "Particle id is greater than the max seen particle id");
  }
}

void ReactionAlgorithm::set_cyl_constraint(double center_x, double center_y,
                                           double radius) {
  auto const &box_geo = *System::get_system().box_geo;
  if (center_x < 0. or center_x > box_geo.length()[0])
    throw std::domain_error("center_x is outside the box");
  if (center_y < 0. or center_y > box_geo.length()[1])
    throw std::domain_error("center_y is outside the box");
  if (radius < 0.)
    throw std::domain_error("radius is invalid");
  m_cyl_x = center_x;
  m_cyl_y = center_y;
  m_cyl_radius = radius;
  m_reaction_constraint = ReactionConstraint::CYL_Z;
}

void ReactionAlgorithm::set_slab_constraint(double slab_start_z,
                                            double slab_end_z) {
  auto const &box_geo = *System::get_system().box_geo;
  if (slab_start_z < 0. or slab_start_z > box_geo.length()[2])
    throw std::domain_error("slab_start_z is outside the box");
  if (slab_end_z < 0. or slab_end_z > box_geo.length()[2])
    throw std::domain_error("slab_end_z is outside the box");
  if (slab_end_z < slab_start_z)
    throw std::domain_error("slab_end_z must be >= slab_start_z");
  m_slab_start_z = slab_start_z;
  m_slab_end_z = slab_end_z;
  m_reaction_constraint = ReactionConstraint::SLAB_Z;
}

/**
 * Writes a random position inside the central box into the provided array.
 */
Utils::Vector3d ReactionAlgorithm::get_random_position_in_box() {
  auto const &box_geo = *System::get_system().box_geo;
  Utils::Vector3d out_pos{};

  if (m_reaction_constraint == ReactionConstraint::CYL_Z) {
    // see http://mathworld.wolfram.com/DiskPointPicking.html
    // for uniform disk point picking in cylinder
    auto const random_radius =
        m_cyl_radius * std::sqrt(m_uniform_real_distribution(m_generator));
    auto const random_phi =
        2. * std::numbers::pi * m_uniform_real_distribution(m_generator);
    out_pos[0] = m_cyl_x + random_radius * cos(random_phi);
    out_pos[1] = m_cyl_y + random_radius * sin(random_phi);
    out_pos[2] = box_geo.length()[2] * m_uniform_real_distribution(m_generator);
  } else if (m_reaction_constraint == ReactionConstraint::SLAB_Z) {
    out_pos[0] = box_geo.length()[0] * m_uniform_real_distribution(m_generator);
    out_pos[1] = box_geo.length()[1] * m_uniform_real_distribution(m_generator);
    out_pos[2] = m_slab_start_z + (m_slab_end_z - m_slab_start_z) *
                                      m_uniform_real_distribution(m_generator);
  } else {
    assert(m_reaction_constraint == ReactionConstraint::NONE);
    out_pos[0] = box_geo.length()[0] * m_uniform_real_distribution(m_generator);
    out_pos[1] = box_geo.length()[1] * m_uniform_real_distribution(m_generator);
    out_pos[2] = box_geo.length()[2] * m_uniform_real_distribution(m_generator);
  }
  return out_pos;
}

/**
 * Creates a particle at the end of the observed particle id range.
 */
int ReactionAlgorithm::create_particle(int p_type) {
  int p_id;
  if (!m_empty_p_ids_smaller_than_max_seen_particle.empty()) {
    auto p_id_iter = std::min_element(
        std::begin(m_empty_p_ids_smaller_than_max_seen_particle),
        std::end(m_empty_p_ids_smaller_than_max_seen_particle));
    p_id = *p_id_iter;
    m_empty_p_ids_smaller_than_max_seen_particle.erase(p_id_iter);
  } else {
    p_id = ::get_maximal_particle_id() + 1;
  }

  // create random velocity vector according to Maxwell-Boltzmann distribution
  auto pos = get_random_position_in_box();
  auto vel = get_random_velocity_vector();

  ::make_new_particle(p_id, pos);
  if (auto p = get_local_particle(p_id)) {
    p->v() = std::sqrt(kT / p->mass()) * vel;
    p->type() = p_type;
#ifdef ELECTROSTATICS
    p->q() = charges_of_types.at(p_type);
#endif
  }
  on_particle_type_change(p_id, ::type_tracking::new_part, p_type);
  return p_id;
}

void ReactionAlgorithm::displacement_mc_move(int type, int n_particles) {
  auto &bookkeeping = make_new_system_state();
  // draw particle ids at random without replacement
  int p_id = -1;
  std::vector<int> drawn_pids{p_id};
  for (int i = 0; i < n_particles; i++) {
    // draw a new particle id
    while (Utils::contains(drawn_pids, p_id)) {
      auto const random_index = i_random(number_of_particles_with_type(type));
      p_id = get_random_p_id(type, random_index);
    }
    drawn_pids.emplace_back(p_id);
    // write new position and new velocity
    typename decltype(ParticleChanges::moved)::value_type old_state;
    auto const new_pos = get_random_position_in_box();
    auto vel = get_random_velocity_vector();
    if (auto p = get_real_particle(p_id)) {
      old_state = {p_id, p->pos(), p->v()};
      p->v() = std::sqrt(kT / p->mass()) * vel;
      if (m_comm.rank() != 0) {
        m_comm.send(0, 42, old_state);
      }
    } else if (m_comm.rank() == 0) {
      m_comm.recv(boost::mpi::any_source, 42, old_state);
    }
    boost::mpi::broadcast(m_comm, old_state, 0);
    bookkeeping.moved.emplace_back(old_state);
    ::set_particle_pos(p_id, new_pos);

    check_exclusion_range(p_id, type);
    if (particle_inside_exclusion_range_touched) {
      break;
    }
  }
}

bool ReactionAlgorithm::make_displacement_mc_move_attempt(int type,
                                                          int n_particles) {

  if (type < 0) {
    throw std::domain_error("Parameter 'type_mc' must be >= 0");
  }
  if (n_particles < 0) {
    throw std::domain_error(
        "Parameter 'particle_number_to_be_changed' must be >= 0");
  }

  if (n_particles == 0) {
    // reject
    return false;
  }

  m_tried_configurational_MC_moves += 1;
  particle_inside_exclusion_range_touched = false;

  auto const n_particles_of_type = ::number_of_particles_with_type(type);
  if (n_particles > n_particles_of_type) {
    // reject
    return false;
  }

  auto const E_pot_old = calculate_potential_energy();
  displacement_mc_move(type, n_particles);
  auto const E_pot_new = (particle_inside_exclusion_range_touched)
                             ? std::numeric_limits<double>::max()
                             : calculate_potential_energy();

  // Metropolis algorithm since proposal density is symmetric
  auto const bf = std::min(1., std::exp(-(E_pot_new - E_pot_old) / kT));

  // // correct for enhanced proposal of small radii by using the
  // // Metropolis-Hastings algorithm for asymmetric proposal densities
  // double old_radius =
  //     std::sqrt(std::pow(particle_positions[0][0]-cyl_x,2) +
  //               std::pow(particle_positions[0][1]-cyl_y,2));
  // double new_radius =
  //     std::sqrt(std::pow(new_pos[0]-cyl_x,2)+std::pow(new_pos[1]-cyl_y,2));
  // auto const bf = std::min(1.0,
  //     exp(-beta*(E_pot_new-E_pot_old))*new_radius/old_radius);

  // Metropolis-Hastings algorithm for asymmetric proposal density
  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    m_accepted_configurational_MC_moves += 1;
    clear_old_system_state();
    return true;
  }
  // reject: restore original particle properties
  restore_old_system_state();
  return false;
}

/**
 * Cleans the list of empty pids and searches for empty pid in the system
 */
void ReactionAlgorithm::setup_bookkeeping_of_empty_pids() {
  // Clean-up the list of empty pids
  m_empty_p_ids_smaller_than_max_seen_particle.clear();

  auto particle_ids = get_particle_ids_parallel();
  std::sort(particle_ids.begin(), particle_ids.end());
  auto pid1 = -1;
  for (auto pid2 : particle_ids) {
    for (int pid = pid1 + 1; pid < pid2; ++pid) {
      m_empty_p_ids_smaller_than_max_seen_particle.push_back(pid);
    }
    pid1 = pid2;
  }
}

double ReactionAlgorithm::calculate_potential_energy() const {
  auto &system = System::get_system();
  auto const obs = system.calculate_energy();
  auto pot = obs->accumulate(-obs->kinetic[0]);
  boost::mpi::broadcast(m_comm, pot, 0);
  return pot;
}

Particle *ReactionAlgorithm::get_real_particle(int p_id) const {
  assert(p_id >= 0);
  auto const &system = System::get_system();
  auto ptr = system.cell_structure->get_local_particle(p_id);
  if (ptr != nullptr and ptr->is_ghost()) {
    ptr = nullptr;
  }
  assert(boost::mpi::all_reduce(m_comm, static_cast<int>(ptr != nullptr),
                                std::plus<>()) == 1);
  return ptr;
}

Particle *ReactionAlgorithm::get_local_particle(int p_id) const {
  assert(p_id >= 0);
  auto const &system = System::get_system();
  auto ptr = system.cell_structure->get_local_particle(p_id);
  assert(boost::mpi::all_reduce(
             m_comm, static_cast<int>(ptr != nullptr and not ptr->is_ghost()),
             std::plus<>()) == 1);
  return ptr;
}

} // namespace ReactionMethods
