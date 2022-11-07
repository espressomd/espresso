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

#include "analysis/statistics.hpp"
#include "cells.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace ReactionMethods {

/**
 * Performs a randomly selected reaction in the reaction ensemble
 */
void ReactionAlgorithm::do_reaction(int reaction_steps) {
  // calculate potential energy; only consider potential energy since we
  // assume that the kinetic part drops out in the process of calculating
  // ensemble averages (kinetic part may be separated and crossed out)
  auto current_E_pot = mpi_calculate_potential_energy();
  for (int i = 0; i < reaction_steps; i++) {
    int reaction_id = i_random(static_cast<int>(reactions.size()));
    generic_oneway_reaction(*reactions[reaction_id], current_E_pot);
  }
}

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
 * Checks whether all necessary variables for the reaction ensemble have been
 * set.
 */
void ReactionAlgorithm::check_reaction_method() const {
  if (reactions.empty()) {
    throw std::runtime_error("Reaction system not initialized");
  }

#ifdef ELECTROSTATICS
  // check for the existence of default charges for all types that take part in
  // the reactions

  for (const auto &current_reaction : reactions) {
    // check for reactants
    for (int reactant_type : current_reaction->reactant_types) {
      auto it = charges_of_types.find(reactant_type);
      if (it == charges_of_types.end()) {
        std::string message = std::string("Forgot to assign charge to type ") +
                              std::to_string(reactant_type);
        throw std::runtime_error(message);
      }
    }
    // check for products
    for (int product_type : current_reaction->product_types) {
      auto it = charges_of_types.find(product_type);
      if (it == charges_of_types.end()) {
        std::string message = std::string("Forgot to assign charge to type ") +
                              std::to_string(product_type);
        throw std::runtime_error(message);
      }
    }
  }
#endif
}

/**
 * Automatically sets the volume which is used by the reaction ensemble to the
 * volume of a cuboid box.
 */
void ReactionAlgorithm::update_volume() { volume = box_geo.volume(); }

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

/**
 * Stores the particle property of a random particle of the provided type into
 * the provided vector
 */
void ReactionAlgorithm::append_particle_property_of_random_particle(
    int type, std::vector<StoredParticleProperty> &list_of_particles) {
  int random_index_in_type_map = i_random(number_of_particles_with_type(type));
  int p_id = get_random_p_id(type, random_index_in_type_map);
  StoredParticleProperty properties = {p_id, type, charges_of_types[type]};
  list_of_particles.push_back(properties);
}

/**
 *Performs a trial reaction move
 */
std::tuple<std::vector<StoredParticleProperty>, std::vector<int>,
           std::vector<StoredParticleProperty>>
ReactionAlgorithm::make_reaction_attempt(
    SingleReaction const &current_reaction) {
  // create or hide particles of types with corresponding types in reaction
  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;

  for (int i = 0; i < std::min(current_reaction.product_types.size(),
                               current_reaction.reactant_types.size());
       i++) {
    // change std::min(reactant_coefficients(i),product_coefficients(i)) many
    // particles of reactant_types(i) to product_types(i)
    for (int j = 0; j < std::min(current_reaction.product_coefficients[i],
                                 current_reaction.reactant_coefficients[i]);
         j++) {
      append_particle_property_of_random_particle(
          current_reaction.reactant_types[i], changed_particles_properties);
      replace_particle(changed_particles_properties.back().p_id,
                       current_reaction.product_types[i]);
    }
    // create product_coefficients(i)-reactant_coefficients(i) many product
    // particles iff product_coefficients(i)-reactant_coefficients(i)>0,
    // iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number
    // of reactant particles
    if (current_reaction.product_coefficients[i] -
            current_reaction.reactant_coefficients[i] >
        0) {
      for (int j = 0; j < current_reaction.product_coefficients[i] -
                              current_reaction.reactant_coefficients[i];
           j++) {
        int p_id = create_particle(current_reaction.product_types[i]);
        check_exclusion_range(p_id);
        p_ids_created_particles.push_back(p_id);
      }
    } else if (current_reaction.reactant_coefficients[i] -
                   current_reaction.product_coefficients[i] >
               0) {
      for (int j = 0; j < current_reaction.reactant_coefficients[i] -
                              current_reaction.product_coefficients[i];
           j++) {
        append_particle_property_of_random_particle(
            current_reaction.reactant_types[i], hidden_particles_properties);
        check_exclusion_range(hidden_particles_properties.back().p_id);
        hide_particle(hidden_particles_properties.back().p_id);
      }
    }
  }
  // create or hide particles of types with noncorresponding replacement types
  for (auto i = std::min(current_reaction.product_types.size(),
                         current_reaction.reactant_types.size());
       i < std::max(current_reaction.product_types.size(),
                    current_reaction.reactant_types.size());
       i++) {
    if (current_reaction.product_types.size() <
        current_reaction.reactant_types.size()) {
      // hide superfluous reactant_types particles
      for (int j = 0; j < current_reaction.reactant_coefficients[i]; j++) {
        append_particle_property_of_random_particle(
            current_reaction.reactant_types[i], hidden_particles_properties);
        check_exclusion_range(hidden_particles_properties.back().p_id);
        hide_particle(hidden_particles_properties.back().p_id);
      }
    } else {
      // create additional product_types particles
      for (int j = 0; j < current_reaction.product_coefficients[i]; j++) {
        int p_id = create_particle(current_reaction.product_types[i]);
        check_exclusion_range(p_id);
        p_ids_created_particles.push_back(p_id);
      }
    }
  }

  return {changed_particles_properties, p_ids_created_particles,
          hidden_particles_properties};
}

/**
 * Restores the previously stored particle properties. This function is invoked
 * when a reaction attempt is rejected.
 */
void ReactionAlgorithm::restore_properties(
    std::vector<StoredParticleProperty> const &property_list) {
  // this function restores all properties of all particles provided in the
  // property list, the format of the property list is (p_id,charge,type)
  // repeated for each particle that occurs in that list
  for (auto &i : property_list) {
    int type = i.type;
#ifdef ELECTROSTATICS
    // set charge
    double charge = i.charge;
    set_particle_q(i.p_id, charge);
#endif
    // set type
    set_particle_type(i.p_id, type);
  }
}

std::map<int, int> ReactionAlgorithm::save_old_particle_numbers(
    SingleReaction const &current_reaction) const {
  std::map<int, int> old_particle_numbers;
  // reactants
  for (int type : current_reaction.reactant_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }

  // products
  for (int type : current_reaction.product_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }
  return old_particle_numbers;
}

void ReactionAlgorithm::generic_oneway_reaction(
    SingleReaction &current_reaction, double &E_pot_old) {

  current_reaction.tried_moves += 1;
  particle_inside_exclusion_range_touched = false;
  if (!all_reactant_particles_exist(current_reaction)) {
    // makes sure, no incomplete reaction is performed -> only need to consider
    // rollback of complete reactions
    return;
  }

  // find reacting molecules in reactants and save their properties for later
  // recreation if step is not accepted
  // do reaction
  // save old particle_numbers
  std::map<int, int> old_particle_numbers =
      save_old_particle_numbers(current_reaction);

  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;

  std::tie(changed_particles_properties, p_ids_created_particles,
           hidden_particles_properties) =
      make_reaction_attempt(current_reaction);

  auto const E_pot_new = (particle_inside_exclusion_range_touched)
                             ? std::numeric_limits<double>::max()
                             : mpi_calculate_potential_energy();

  auto const bf = calculate_acceptance_probability(
      current_reaction, E_pot_old, E_pot_new, old_particle_numbers);

  std::vector<double> exponential = {exp(-1.0 / kT * (E_pot_new - E_pot_old))};
  current_reaction.accumulator_potential_energy_difference_exponential(
      exponential);

  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    // delete hidden reactant_particles (remark: don't delete changed particles)
    // extract ids of to be deleted particles
    auto len_hidden_particles_properties =
        static_cast<int>(hidden_particles_properties.size());
    std::vector<int> to_be_deleted_hidden_ids(len_hidden_particles_properties);
    std::vector<int> to_be_deleted_hidden_types(
        len_hidden_particles_properties);
    for (int i = 0; i < len_hidden_particles_properties; i++) {
      auto p_id = hidden_particles_properties[i].p_id;
      to_be_deleted_hidden_ids[i] = p_id;
      to_be_deleted_hidden_types[i] = hidden_particles_properties[i].type;
      // change back type otherwise the bookkeeping algorithm is not working
      set_particle_type(p_id, hidden_particles_properties[i].type);
    }

    for (int i = 0; i < len_hidden_particles_properties; i++) {
      delete_particle(to_be_deleted_hidden_ids[i]); // delete particle
    }
    current_reaction.accepted_moves += 1;
    E_pot_old = E_pot_new; // Update the system energy

  } else {
    // reject
    // reverse reaction
    // 1) delete created product particles
    for (int p_ids_created_particle : p_ids_created_particles) {
      delete_particle(p_ids_created_particle);
    }
    // 2) restore previously hidden reactant particles
    restore_properties(hidden_particles_properties);
    // 3) restore previously changed reactant particles
    restore_properties(changed_particles_properties);
  }
}

/**
 * Replaces a particle with the given particle id to be of a certain type. This
 * especially means that the particle type and the particle charge are changed.
 */
void ReactionAlgorithm::replace_particle(int p_id, int desired_type) const {
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, charges_of_types.at(desired_type));
#endif
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
void ReactionAlgorithm::hide_particle(int p_id) const {
  set_particle_type(p_id, non_interacting_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, 0.0);
#endif
}

/**
 * Check if the inserted particle is too close to neighboring particles.
 */
void ReactionAlgorithm::check_exclusion_range(int inserted_particle_id) {

  auto const &inserted_particle = get_particle_data(inserted_particle_id);

  /* Check the exclusion radius of the inserted particle */
  if (exclusion_radius_per_type.count(inserted_particle.type()) != 0) {
    if (exclusion_radius_per_type[inserted_particle.type()] == 0.) {
      return;
    }
  }

  std::vector<int> particle_ids;
  if (neighbor_search_order_n) {
    particle_ids = get_particle_ids();
    /* remove the inserted particle id */
    particle_ids.erase(std::remove(particle_ids.begin(), particle_ids.end(),
                                   inserted_particle_id),
                       particle_ids.end());
  } else {
    particle_ids = mpi_get_short_range_neighbors(inserted_particle.id(),
                                                 m_max_exclusion_range);
  }

  /* Check if the inserted particle within the exclusion radius of any other
   * particle */
  for (auto const &particle_id : particle_ids) {
    auto const &p = get_particle_data(particle_id);
    double excluded_distance;
    if (exclusion_radius_per_type.count(inserted_particle.type()) == 0 ||
        exclusion_radius_per_type.count(p.type()) == 0) {
      excluded_distance = exclusion_range;
    } else if (exclusion_radius_per_type[p.type()] == 0.) {
      continue;
    } else {
      excluded_distance = exclusion_radius_per_type[inserted_particle.type()] +
                          exclusion_radius_per_type[p.type()];
    }

    auto const d_min =
        box_geo.get_mi_vector(p.pos(), inserted_particle.pos()).norm();

    if (d_min < excluded_distance) {
      particle_inside_exclusion_range_touched = true;
      break;
    }
  }
}

/**
 * Deletes the particle with the given p_id and stores the id if the deletion
 * created a hole in the particle id range. This method is intended to only
 * delete unbonded particles since bonds are coupled to ids. This is used to
 * avoid the id range becoming excessively huge.
 */
void ReactionAlgorithm::delete_particle(int p_id) {
  auto const old_max_seen_id = get_maximal_particle_id();
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
  Utils::Vector3d out_pos{};

  if (m_reaction_constraint == ReactionConstraint::CYL_Z) {
    // see http://mathworld.wolfram.com/DiskPointPicking.html
    // for uniform disk point picking in cylinder
    auto const random_radius =
        m_cyl_radius * std::sqrt(m_uniform_real_distribution(m_generator));
    auto const random_phi =
        2. * Utils::pi() * m_uniform_real_distribution(m_generator);
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
int ReactionAlgorithm::create_particle(int desired_type) {
  int p_id;
  if (!m_empty_p_ids_smaller_than_max_seen_particle.empty()) {
    auto p_id_iter = std::min_element(
        std::begin(m_empty_p_ids_smaller_than_max_seen_particle),
        std::end(m_empty_p_ids_smaller_than_max_seen_particle));
    p_id = *p_id_iter;
    m_empty_p_ids_smaller_than_max_seen_particle.erase(p_id_iter);
  } else {
    p_id = get_maximal_particle_id() + 1;
  }

  // we use mass=1 for all particles, think about adapting this
  move_particle(p_id, get_random_position_in_box(), std::sqrt(kT));
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, charges_of_types[desired_type]);
#endif
  return p_id;
}

void ReactionAlgorithm::move_particle(int p_id, Utils::Vector3d const &new_pos,
                                      double velocity_prefactor) {
  place_particle(p_id, new_pos);
  // create random velocity vector according to Maxwell-Boltzmann distribution
  Utils::Vector3d vel;
  vel[0] = velocity_prefactor * m_normal_distribution(m_generator);
  vel[1] = velocity_prefactor * m_normal_distribution(m_generator);
  vel[2] = velocity_prefactor * m_normal_distribution(m_generator);
  set_particle_v(p_id, vel);
}

/**
 * Performs particle displacement MC moves in the canonical ensemble
 * @param mc_steps Number of trial MC steps
 * @param particle_types_to_move  List of particle types from which particles
 * are selected. If empty, any particle can be chosen by default.
 *
 */

void ReactionAlgorithm::do_particle_displacement_MC_move(
    int mc_steps, std::vector<int> particle_types_to_move) {
  std::vector<signed int> ids_to_move;

  if (particle_types_to_move.empty()) {
    ids_to_move = get_particle_ids();
  } else {
    for (signed int p_id : get_particle_ids()) {
      if (std::find(
              particle_types_to_move.begin(), particle_types_to_move.end(),
              get_particle_data(p_id).type()) != particle_types_to_move.end()) {
        ids_to_move.push_back(p_id);
      }
    }
  }

  // If there are no particles available to move, return
  if (ids_to_move.empty()) {
    return;
  }

  for (int step = 0; step < mc_steps; step++) {

    N_trial_particle_displacement_MC_moves += 1;
    particle_inside_exclusion_range_touched = false;
    int p_id = ids_to_move[i_random(static_cast<int>(ids_to_move.size()))];
    auto const particle = get_particle_data(p_id);
    auto const old_position = particle.pos();
    auto const E_pot_old = mpi_calculate_potential_energy();
    auto const new_position = get_random_position_in_box();

    place_particle(p_id, new_position);
    check_exclusion_range(p_id);

    auto const E_pot_new = (particle_inside_exclusion_range_touched)
                               ? std::numeric_limits<double>::max()
                               : mpi_calculate_potential_energy();

    // Metropolis algorithm
    auto const bf = std::min(1.0, exp(-1.0 * (E_pot_new - E_pot_old) / kT));

    if (m_uniform_real_distribution(m_generator) < bf) {
      // accept
      N_accepted_particle_displacement_MC_moves += 1;
    } else {
      // reject: restore original particle position
      place_particle(p_id, old_position);
    }
  }
}

} // namespace ReactionMethods
