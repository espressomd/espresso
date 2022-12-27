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

#include "random.hpp"

#include <utils/Vector.hpp>

#include <functional>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ReactionMethods {

/** @brief Bookkeeping for particle changes. */
struct ParticleChangeRecorder {
  explicit ParticleChangeRecorder(std::function<void(int)> &&deleter)
      : m_deleter{deleter} {}

  struct StoredParticleProperty {
    int p_id;
    int type;
    double charge;
  };

  /** @brief Record particle creation. */
  void save_created_particle(int p_id) { m_created.push_back(p_id); }

  /** @brief Record particle state before it is hidden. */
  void save_hidden_particle(StoredParticleProperty &&item) {
    m_hidden.push_back(item);
  }

  /** @brief Record particle state before a property change. */
  void save_changed_particle(StoredParticleProperty &&item) {
    m_changed.push_back(item);
  }

  /** @brief Delete hidden particles from the system. */
  void delete_hidden_particles() const;

  /** @brief Restore original system state. */
  void restore_original_state() const;

private:
  std::function<void(int)> m_deleter;
  std::vector<int> m_created;
  std::vector<StoredParticleProperty> m_hidden;
  std::vector<StoredParticleProperty> m_changed;
};

/** Base class for reaction ensemble methods */
class ReactionAlgorithm {

public:
  ReactionAlgorithm(
      int seed, double kT, double exclusion_range,
      std::unordered_map<int, double> const &exclusion_radius_per_type)
      : kT{kT}, exclusion_range{exclusion_range},
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
  std::map<int, double> charges_of_types;
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

  void do_reaction(int reaction_steps);
  void check_reaction_method() const;
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

  bool displacement_move_for_particles_of_type(int type, int n_part);

  bool particle_inside_exclusion_range_touched = false;
  bool neighbor_search_order_n = true;

protected:
  std::vector<int> m_empty_p_ids_smaller_than_max_seen_particle;
  /**
   * @brief Carry out a generic one-way chemical reaction.
   *
   * Generic one way reaction of the type
   * <tt>A+B+...+G +... --> K+...X + Z +...</tt>
   * You need to use <tt>2A --> B</tt> instead of <tt>A+A --> B</tt> since
   * in the latter you assume distinctness of the particles, however both
   * ways to describe the reaction are equivalent in the thermodynamic limit
   * (large particle numbers). Furthermore, it is crucial for the function
   * in which order you provide the reactant and product types since particles
   * will be replaced correspondingly! If there are less reactants than
   * products, new product particles are created randomly in the box.
   * Matching particles simply change the types. If there are more reactants
   * than products, old reactant particles are deleted.
   *
   * @param[in,out] current_reaction  The reaction to attempt.
   * @param[in,out] E_pot_old         The current potential energy.
   */
  void generic_oneway_reaction(SingleReaction &current_reaction,
                               double &E_pot_old);

  ParticleChangeRecorder make_reaction_attempt(SingleReaction const &reaction);
  std::vector<std::tuple<int, Utils::Vector3d, Utils::Vector3d>>
  generate_new_particle_positions(int type, int n_particles);
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
  bool
  all_reactant_particles_exist(SingleReaction const &current_reaction) const;

  virtual double
  calculate_acceptance_probability(SingleReaction const &, double, double,
                                   std::map<int, int> const &) const {
    return -10.;
  }

private:
  std::mt19937 m_generator;
  std::normal_distribution<double> m_normal_distribution;
  std::uniform_real_distribution<double> m_uniform_real_distribution;

  std::map<int, int>
  save_old_particle_numbers(SingleReaction const &current_reaction) const;

  void replace_particle(int p_id, int desired_type) const;
  int create_particle(int desired_type);
  void hide_particle(int p_id) const;
  void check_exclusion_range(int inserted_particle_id);
  void move_particle(int p_id, Utils::Vector3d const &new_pos,
                     double velocity_prefactor);

  enum class ReactionConstraint { NONE, CYL_Z, SLAB_Z };
  ReactionConstraint m_reaction_constraint = ReactionConstraint::NONE;
  double m_cyl_radius = -10.0;
  double m_cyl_x = -10.0;
  double m_cyl_y = -10.0;
  double m_slab_start_z = -10.0;
  double m_slab_end_z = -10.0;
  double m_max_exclusion_range = 0.;

protected:
  Utils::Vector3d get_random_position_in_box();
};

} // namespace ReactionMethods
#endif
