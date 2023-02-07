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
#ifndef REACTION_METHODS_REACTION_ALGORITHM_HPP
#define REACTION_METHODS_REACTION_ALGORITHM_HPP

#include "config/config.hpp"

#include "SingleReaction.hpp"

#include "Particle.hpp"
#include "random.hpp"

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace boost::mpi {
class communicator;
} // namespace boost::mpi

namespace ReactionMethods {

/** Base class for reaction ensemble methods */
class ReactionAlgorithm {
private:
  boost::mpi::communicator const &m_comm;

public:
  ReactionAlgorithm(
      boost::mpi::communicator const &comm, int seed, double kT,
      double exclusion_range,
      std::unordered_map<int, double> const &exclusion_radius_per_type)
      : m_comm{comm}, kT{kT}, exclusion_range{exclusion_range},
        m_generator(Random::mt19937(std::seed_seq({seed, seed, seed}))),
        m_normal_distribution(0.0, 1.0), m_uniform_real_distribution(0.0, 1.0) {
    if (kT < 0.) {
      throw std::domain_error("Invalid value for 'kT'");
    }
    if (exclusion_range < 0.) {
      throw std::domain_error("Invalid value for 'exclusion_range'");
    }
    set_exclusion_radius_per_type(exclusion_radius_per_type);
    update_volume();
  }

  virtual ~ReactionAlgorithm() = default;

  std::vector<std::shared_ptr<SingleReaction>> reactions;
  std::unordered_map<int, double> charges_of_types;
  double kT;
  /**
   * Hard sphere radius. If particles are closer than this value,
   * it is assumed that their interaction energy gets approximately
   * infinite, therefore these configurations do not contribute
   * to the partition function and ensemble averages.
   */
  double exclusion_range;
  std::unordered_map<int, double> exclusion_radius_per_type;
  double volume;
  int non_interacting_type = 100;

  int m_accepted_configurational_MC_moves = 0;
  int m_tried_configurational_MC_moves = 0;
  double get_acceptance_rate_configurational_moves() const {
    return static_cast<double>(m_accepted_configurational_MC_moves) /
           static_cast<double>(m_tried_configurational_MC_moves);
  }

  auto get_kT() const { return kT; }
  auto get_exclusion_range() const { return exclusion_range; }
  auto get_volume() const { return volume; }
  void set_volume(double new_volume) {
    if (new_volume <= 0.) {
      throw std::domain_error("Invalid value for 'volume'");
    }
    volume = new_volume;
  }
  void update_volume();
  void
  set_exclusion_radius_per_type(std::unordered_map<int, double> const &map) {
    auto max_exclusion_range = exclusion_range;
    for (auto const &item : map) {
      auto const type = item.first;
      auto const exclusion_radius = item.second;
      if (exclusion_radius < 0.) {
        throw std::domain_error("Invalid excluded_radius value for type " +
                                std::to_string(type) + ": radius " +
                                std::to_string(exclusion_radius));
      }
      max_exclusion_range =
          std::max(max_exclusion_range, 2. * exclusion_radius);
    }
    exclusion_radius_per_type = map;
    m_max_exclusion_range = max_exclusion_range;
  }

  void remove_constraint() { m_reaction_constraint = ReactionConstraint::NONE; }
  void set_cyl_constraint(double center_x, double center_y, double radius);
  void set_slab_constraint(double slab_start_z, double slab_end_z);
  Utils::Vector2d get_slab_constraint_parameters() const {
    if (m_reaction_constraint != ReactionConstraint::SLAB_Z) {
      throw std::runtime_error("no slab constraint is currently active");
    }
    return {m_slab_start_z, m_slab_end_z};
  }

  void setup_bookkeeping_of_empty_pids();
  void delete_particle(int p_id);
  void add_reaction(std::shared_ptr<SingleReaction> const &new_reaction);
  void delete_reaction(int reaction_id) {
    reactions.erase(reactions.begin() + reaction_id);
  }

  bool particle_inside_exclusion_range_touched = false;
  bool neighbor_search_order_n = true;

protected:
  std::vector<int> m_empty_p_ids_smaller_than_max_seen_particle;

public:
  struct ParticleChanges {
    std::vector<int> created{};
    std::vector<std::tuple<int, int>> changed{};
    std::vector<std::tuple<int, int>> hidden{};
    std::vector<std::tuple<int, Utils::Vector3d, Utils::Vector3d>> moved{};
    std::unordered_map<int, int> old_particle_numbers{};
    int reaction_id{-1};
  };

  bool is_reaction_under_way() const { return m_system_changes != nullptr; }
  auto const &get_old_system_state() const {
    if (not is_reaction_under_way()) {
      throw std::runtime_error("No chemical reaction is currently under way");
    }
    return *m_system_changes;
  }

protected:
  /** @brief Restore last valid system state. */
  void restore_old_system_state();
  /** @brief Clear last valid system state. */
  void clear_old_system_state() {
    assert(is_reaction_under_way());
    m_system_changes = nullptr;
  }
  /** @brief Open new handle for system state tracking. */
  auto &make_new_system_state() {
    assert(not is_reaction_under_way());
    m_system_changes = std::make_shared<ParticleChanges>();
    return *m_system_changes;
  }
  std::shared_ptr<ParticleChanges> m_system_changes;

public:
  /**
   * Attempt a reaction move and calculate the new potential energy.
   * Particles are selected without replacement.
   * @returns Potential energy of the system after the move.
   */
  std::optional<double> generic_oneway_reaction_part_1(int reaction_id);
  /**
   * Accept or reject moves made by @ref generic_oneway_reaction_part_1 based
   * on a probability acceptance @c bf.
   * @returns Potential energy of the system after the move was accepted or
   * rejected.
   */
  double generic_oneway_reaction_part_2(int reaction_id, double bf,
                                        double E_pot_old, double E_pot_new);
  /**
   * Attempt a global MC move for particles of a given type.
   * Particles are selected without replacement.
   * @param type          Type of particles to move.
   * @param n_particles   Number of particles of that type to move.
   * @returns true if all moves were accepted.
   */
  bool make_displacement_mc_move_attempt(int type, int n_particles);
  /**
   * Perform a global MC move for particles of a given type.
   * Particles are selected without replacement.
   * @param type          Type of particles to move.
   * @param n_particles   Number of particles of that type to move.
   */
  void displacement_mc_move(int type, int n_particles);

  /** @brief Compute the system potential energy. */
  double calculate_potential_energy() const;

protected:
  /** @brief Carry out a chemical reaction and save the old system state. */
  void make_reaction_attempt(::ReactionMethods::SingleReaction const &reaction,
                             ParticleChanges &bookkeeping);

public:
  /**
   * @brief draws a random integer from the uniform distribution in the range
   * [0,maxint-1]
   *
   * @param maxint range.
   */
  int i_random(int maxint) {
    std::uniform_int_distribution<int> uniform_int_dist(0, maxint - 1);
    return uniform_int_dist(m_generator);
  }

protected:
  bool
  all_reactant_particles_exist(SingleReaction const &current_reaction) const;

private:
  std::mt19937 m_generator;
  std::normal_distribution<double> m_normal_distribution;
  std::uniform_real_distribution<double> m_uniform_real_distribution;

  std::unordered_map<int, int>
  get_particle_numbers(SingleReaction const &reaction) const;

  int create_particle(int p_type);
  void hide_particle(int p_id, int p_type) const;
  void check_exclusion_range(int p_id, int p_type);
  auto get_random_uniform_number() {
    return m_uniform_real_distribution(m_generator);
  }
  auto get_random_velocity_vector() {
    return Utils::Vector3d{m_normal_distribution(m_generator),
                           m_normal_distribution(m_generator),
                           m_normal_distribution(m_generator)};
  }

  enum class ReactionConstraint { NONE, CYL_Z, SLAB_Z };
  ReactionConstraint m_reaction_constraint = ReactionConstraint::NONE;
  double m_cyl_radius = -10.0;
  double m_cyl_x = -10.0;
  double m_cyl_y = -10.0;
  double m_slab_start_z = -10.0;
  double m_slab_end_z = -10.0;
  double m_max_exclusion_range = 0.;

  Particle *get_real_particle(int p_id) const;
  Particle *get_local_particle(int p_id) const;

protected:
  Utils::Vector3d get_random_position_in_box();
};

} // namespace ReactionMethods
#endif
