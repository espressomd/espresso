/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "config.hpp"

#include "reaction_methods/ReactionAlgorithm.hpp"

#include "energy.hpp"
#include "grid.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "statistics.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ReactionMethods {

/**
 * Performs a randomly selected reaction in the reaction ensemble
 */
int ReactionAlgorithm::do_reaction(int reaction_steps) {
  for (int i = 0; i < reaction_steps; i++) {
    int reaction_id = i_random(static_cast<int>(reactions.size()));
    generic_oneway_reaction(reaction_id);
  }
  return 0;
}

/**
 * Adds a reaction to the reaction system
 */
void ReactionAlgorithm::add_reaction(
    double gamma, const std::vector<int> &reactant_types,
    const std::vector<int> &reactant_coefficients,
    const std::vector<int> &product_types,
    const std::vector<int> &product_coefficients) {
  SingleReaction new_reaction(gamma, reactant_types, reactant_coefficients,
                              product_types, product_coefficients);

  // make ESPResSo count the particle numbers which take part in the reactions
  for (int reactant_type : new_reaction.reactant_types)
    init_type_map(reactant_type);
  for (int product_type : new_reaction.product_types)
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

  if (temperature < 0) {
    throw std::runtime_error("Temperatures cannot be negative. Please provide "
                             "a temperature (in k_B T) to the simulation. "
                             "Normally it should be 1.0. This will be used "
                             "directly to calculate beta:=1/(k_B T) which "
                             "occurs in the exp(-beta*E)\n");
  }
#ifdef ELECTROSTATICS
  // check for the existence of default charges for all types that take part in
  // the reactions

  for (const auto &current_reaction : reactions) {
    // check for reactants
    for (int reactant_type : current_reaction.reactant_types) {
      auto it = charges_of_types.find(reactant_type);
      if (it == charges_of_types.end()) {
        std::string message = std::string("Forgot to assign charge to type ") +
                              std::to_string(reactant_type);
        throw std::runtime_error(message);
      }
    }
    // check for products
    for (int product_type : current_reaction.product_types) {
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
 * volume of a cuboid box, if not already initialized with another value.
 */
void ReactionAlgorithm::set_cuboid_reaction_ensemble_volume() {
  if (volume < 0)
    volume = box_geo.volume();
}

/**
 * Checks whether all particles exist for the provided reaction.
 */
bool ReactionAlgorithm::all_reactant_particles_exist(int reaction_id) const {
  bool enough_particles = true;
  for (int i = 0; i < reactions[reaction_id].reactant_types.size(); i++) {
    int current_number =
        number_of_particles_with_type(reactions[reaction_id].reactant_types[i]);
    if (current_number < reactions[reaction_id].reactant_coefficients[i]) {
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
  StoredParticleProperty property_of_part = {p_id, charges_of_types[type],
                                             type};
  list_of_particles.push_back(property_of_part);
}

/**
 *Performs a trial reaction move
 */
void ReactionAlgorithm::make_reaction_attempt(
    SingleReaction const &current_reaction,
    std::vector<StoredParticleProperty> &changed_particles_properties,
    std::vector<int> &p_ids_created_particles,
    std::vector<StoredParticleProperty> &hidden_particles_properties) {
  // create or hide particles of types with corresponding types in reaction
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
        hide_particle(hidden_particles_properties.back().p_id,
                      current_reaction.reactant_types[i]);
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
        hide_particle(hidden_particles_properties.back().p_id,
                      current_reaction.reactant_types[i]);
      }
    } else {
      // create additional product_types particles
      for (int j = 0; j < current_reaction.product_coefficients[i]; j++) {
        int p_id = create_particle(current_reaction.product_types[i]);
        p_ids_created_particles.push_back(p_id);
      }
    }
  }
}

/**
 * Restores the previously stored particle properties. This function is invoked
 * when a reaction attempt is rejected.
 */
void ReactionAlgorithm::restore_properties(
    std::vector<StoredParticleProperty> const &property_list,
    const int number_of_saved_properties) {
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

std::map<int, int>
ReactionAlgorithm::save_old_particle_numbers(int reaction_id) {
  std::map<int, int> old_particle_numbers;
  // reactants
  for (int type : reactions[reaction_id].reactant_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }

  // products
  for (int type : reactions[reaction_id].product_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }
  return old_particle_numbers;
}

/**
 * Generic one way reaction
 * A+B+...+G +... --> K+...X + Z +...
 * You need to use 2A --> B instead of A+A --> B since in the last case you
 * assume distinctness of the particles, however both ways to describe the
 * reaction are equivalent in the thermodynamic limit (large particle numbers).
 * Furthermore, it is crucial for the function in which order you provide the
 * reactant and product types since particles will be replaced correspondingly!
 * If there are less reactants than products, new product particles are created
 * randomly in the box. Matching particles simply change the types. If there
 * are more reactants than products, old reactant particles are deleted.
 */
void ReactionAlgorithm::generic_oneway_reaction(int reaction_id) {

  SingleReaction &current_reaction = reactions[reaction_id];
  current_reaction.tried_moves += 1;
  particle_inside_exclusion_radius_touched = false;
  int old_state_index = -1; // for Wang-Landau algorithm
  on_reaction_entry(old_state_index);
  if (!all_reactant_particles_exist(reaction_id)) {
    // makes sure, no incomplete reaction is performed -> only need to consider
    // rollback of complete reactions
    on_reaction_rejection_directly_after_entry(old_state_index);
    return;
  }

  // calculate potential energy
  const double E_pot_old =
      calculate_current_potential_energy_of_system(); // only consider potential
                                                      // energy since we assume
                                                      // that the kinetic part
                                                      // drops out in the
                                                      // process of calculating
                                                      // ensemble averages
                                                      // (kinetic part may be
                                                      // separated and crossed
                                                      // out)

  // find reacting molecules in reactants and save their properties for later
  // recreation if step is not accepted
  // do reaction
  // save old particle_numbers
  std::map<int, int> old_particle_numbers =
      save_old_particle_numbers(reaction_id);

  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;
  const int number_of_saved_properties =
      3; // save p_id, charge and type of the reactant particle, only thing we
         // need to hide the particle and recover it
  make_reaction_attempt(current_reaction, changed_particles_properties,
                        p_ids_created_particles, hidden_particles_properties);

  double E_pot_new;
  if (particle_inside_exclusion_radius_touched)
    E_pot_new = std::numeric_limits<double>::max();
  else
    E_pot_new = calculate_current_potential_energy_of_system();

  int new_state_index = -1; // save new_state_index for Wang-Landau algorithm
  int accepted_state = -1;  // for Wang-Landau algorithm
  on_attempted_reaction(new_state_index);

  constexpr bool only_make_configuration_changing_move = false;
  auto const bf = calculate_acceptance_probability(
      current_reaction, E_pot_old, E_pot_new, old_particle_numbers,
      old_state_index, new_state_index, only_make_configuration_changing_move);

  std::vector<double> exponential = {
      exp(-1.0 / temperature * (E_pot_new - E_pot_old))};
  current_reaction.accumulator_exponentials(exponential);

  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    accepted_state = new_state_index;

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
  } else {
    // reject
    accepted_state = old_state_index;
    // reverse reaction
    // 1) delete created product particles
    for (int p_ids_created_particle : p_ids_created_particles) {
      delete_particle(p_ids_created_particle);
    }
    // 2) restore previously hidden reactant particles
    restore_properties(hidden_particles_properties, number_of_saved_properties);
    // 3) restore previously changed reactant particles
    restore_properties(changed_particles_properties,
                       number_of_saved_properties);
  }
  on_end_reaction(accepted_state);
}

/**
 * Replaces a particle with the given particle id to be of a certain type. This
 * especially means that the particle type and the particle charge are changed.
 */
void ReactionAlgorithm::replace_particle(int p_id, int desired_type) {
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, charges_of_types[desired_type]);
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
void ReactionAlgorithm::hide_particle(int p_id, int previous_type) {

  auto const part = get_particle_data(p_id);
  auto const d_min = distto(partCfg(), part.r.p, p_id);
  if (d_min < exclusion_radius)
    particle_inside_exclusion_radius_touched = true;

#ifdef ELECTROSTATICS
  // set charge
  set_particle_q(p_id, 0.0);
#endif
  // set type
  set_particle_type(p_id, non_interacting_type);
}

/**
 * Deletes the particle with the given p_id and stores the id if the deletion
 * created a hole in the particle id range. This method is intended to only
 * delete unbonded particles since bonds are coupled to ids. This is used to
 * avoid the id range becoming excessively huge.
 */
int ReactionAlgorithm::delete_particle(int p_id) {
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
  return 0;
}

/**
 * Writes a random position inside the central box into the provided array.
 */
Utils::Vector3d ReactionAlgorithm::get_random_position_in_box() {
  Utils::Vector3d out_pos{};

  if (box_is_cylindric_around_z_axis) {
    // see http://mathworld.wolfram.com/DiskPointPicking.html
    double random_radius =
        cyl_radius *
        std::sqrt(m_uniform_real_distribution(
            m_generator)); // for uniform disk point picking in cylinder
    double phi = 2.0 * Utils::pi() * m_uniform_real_distribution(m_generator);
    out_pos[0] = cyl_x + random_radius * cos(phi);
    out_pos[1] = cyl_y + random_radius * sin(phi);
    out_pos[2] = box_geo.length()[2] * m_uniform_real_distribution(m_generator);
  } else if (box_has_wall_constraints) {
    out_pos[0] = box_geo.length()[0] * m_uniform_real_distribution(m_generator);
    out_pos[1] = box_geo.length()[1] * m_uniform_real_distribution(m_generator);
    out_pos[2] = slab_start_z + (slab_end_z - slab_start_z) *
                                    m_uniform_real_distribution(m_generator);
  } else {
    // cubic case
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

  // create random velocity vector according to Maxwell-Boltzmann distribution
  Utils::Vector3d vel;
  // we use mass=1 for all particles, think about adapting this
  vel[0] = std::sqrt(temperature) * m_normal_distribution(m_generator);
  vel[1] = std::sqrt(temperature) * m_normal_distribution(m_generator);
  vel[2] = std::sqrt(temperature) * m_normal_distribution(m_generator);
#ifdef ELECTROSTATICS
  double charge = charges_of_types[desired_type];
#endif

  auto const pos_vec = get_random_position_in_box();
  place_particle(p_id, pos_vec);
  // set type
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  // set charge
  set_particle_q(p_id, charge);
#endif
  // set velocities
  set_particle_v(p_id, vel);
  double d_min = distto(partCfg(), pos_vec, p_id);
  if (d_min < exclusion_radius) {
    // setting of a minimal distance is allowed to avoid overlapping
    // configurations if there is a repulsive potential. States with
    // very high energies have a probability of almost zero and
    // therefore do not contribute to ensemble averages.
    particle_inside_exclusion_radius_touched = true;
  }
  return p_id;
}

/**
 * Performs a global MC move for a particle of the provided type.
 */
bool ReactionAlgorithm::do_global_mc_move_for_particles_of_type(
    int type, int particle_number_of_type_to_be_changed, bool use_wang_landau) {
  m_tried_configurational_MC_moves += 1;
  particle_inside_exclusion_radius_touched = false;

  int old_state_index = -1;
  if (use_wang_landau) {
    on_reaction_entry(old_state_index);
  }

  int particle_number_of_type = number_of_particles_with_type(type);
  if (particle_number_of_type == 0 or
      particle_number_of_type_to_be_changed == 0) {
    // reject
    if (use_wang_landau) {
      on_mc_rejection_directly_after_entry(old_state_index);
    }
    return false;
  }

  const double E_pot_old = calculate_current_potential_energy_of_system();

  std::vector<Utils::Vector3d> particle_positions(
      particle_number_of_type_to_be_changed);
  std::vector<int> p_id_s_changed_particles;

  // heuristic to find a particle that hasn't been touched already
  int random_index_in_type_map = i_random(number_of_particles_with_type(type));
  int p_id = get_random_p_id(type, random_index_in_type_map);
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++) {
    while (Utils::contains(p_id_s_changed_particles, p_id)) {
      random_index_in_type_map = i_random(number_of_particles_with_type(type));
      p_id = get_random_p_id(
          type, random_index_in_type_map); // check whether you already touched
                                           // this p_id, then reassign
    }

    auto part = get_particle_data(p_id);

    particle_positions[i] = part.r.p;
    p_id_s_changed_particles.push_back(p_id);
  }

  // propose new positions
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++) {
    p_id = p_id_s_changed_particles[i];
    // change particle position
    auto const new_pos = get_random_position_in_box();
    auto const &p = get_particle_data(p_id);
    auto const prefactor = std::sqrt(temperature / p.p.mass);
    Utils::Vector3d vel;
    vel[0] = prefactor * m_normal_distribution(m_generator);
    vel[1] = prefactor * m_normal_distribution(m_generator);
    vel[2] = prefactor * m_normal_distribution(m_generator);
    set_particle_v(p_id, vel);
    place_particle(p_id, new_pos);
    auto const d_min = distto(partCfg(), new_pos, p_id);
    if (d_min < exclusion_radius)
      particle_inside_exclusion_radius_touched = true;
  }

  double E_pot_new;
  if (particle_inside_exclusion_radius_touched)
    E_pot_new = std::numeric_limits<double>::max();
  else
    E_pot_new = calculate_current_potential_energy_of_system();

  double beta = 1.0 / temperature;

  int new_state_index = -1;
  double bf = 1.0;
  std::map<int, int> dummy_old_particle_numbers;
  SingleReaction temp_unimportant_arbitrary_reaction;

  if (use_wang_landau) {
    new_state_index = on_mc_use_WL_get_new_state();
    bf = calculate_acceptance_probability(
        temp_unimportant_arbitrary_reaction, E_pot_old, E_pot_new,
        dummy_old_particle_numbers, old_state_index, new_state_index, true);
  } else {
    // Metropolis algorithm since proposal density is symmetric
    bf = std::min(1.0, bf * exp(-beta * (E_pot_new - E_pot_old)));
  }

  // // correct for enhanced proposal of small radii by using the
  // // Metropolis-Hastings algorithm for asymmetric proposal densities
  // double old_radius =
  //     std::sqrt(std::pow(particle_positions[0][0]-cyl_x,2) +
  //               std::pow(particle_positions[0][1]-cyl_y,2));
  // double new_radius =
  //     std::sqrt(std::pow(new_pos[0]-cyl_x,2)+std::pow(new_pos[1]-cyl_y,2));
  // bf = std::min(1.0,
  //     bf*exp(-beta*(E_pot_new-E_pot_old))*new_radius/old_radius);

  // Metropolis-Hastings algorithm for asymmetric proposal density
  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    m_accepted_configurational_MC_moves += 1;
    if (use_wang_landau) {
      on_mc_accept(new_state_index);
    }
    return true;
  }
  // reject
  // modify wang_landau histogram and potential
  if (use_wang_landau) {
    on_mc_reject(old_state_index);
  }
  // create particles again at the positions they were
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++)
    place_particle(p_id_s_changed_particles[i], particle_positions[i]);
  return false;
}

} // namespace ReactionMethods
