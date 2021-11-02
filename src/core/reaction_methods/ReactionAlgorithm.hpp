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
#ifndef REACTION_METHODS_REACTION_ALGORITHM_HPP
#define REACTION_METHODS_REACTION_ALGORITHM_HPP

#include "config.hpp"

#include "SingleReaction.hpp"

#include "random.hpp"

#include <utils/Vector.hpp>

#include <map>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

namespace ReactionMethods {

struct StoredParticleProperty {
  int p_id;
  double charge;
  int type;
};

/** Base class for reaction ensemble methods */
class ReactionAlgorithm {

public:
  ReactionAlgorithm(int seed)
      : m_generator(Random::mt19937(std::seed_seq({seed, seed, seed}))),
        m_normal_distribution(0.0, 1.0), m_uniform_real_distribution(0.0, 1.0) {
  }

  virtual ~ReactionAlgorithm() = default;

  std::vector<SingleReaction> reactions;
  std::map<int, double> charges_of_types;
  double kT = -10.0;
  /**
   * Hard sphere radius. If particles are closer than this value,
   * it is assumed that their interaction energy gets approximately
   * infinite, therefore these configurations do not contribute
   * to the partition function and ensemble averages.
   */
  double exclusion_radius = 0.0;
  double volume = -10.0;
  int non_interacting_type = 100;

  int m_accepted_configurational_MC_moves = 0;
  int m_tried_configurational_MC_moves = 0;
  double get_acceptance_rate_configurational_moves() const {
    return static_cast<double>(m_accepted_configurational_MC_moves) /
           static_cast<double>(m_tried_configurational_MC_moves);
  }

  void set_cuboid_reaction_ensemble_volume();
  virtual int do_reaction(int reaction_steps);
  void check_reaction_method() const;
  void remove_constraint() { m_reaction_constraint = ReactionConstraint::NONE; }
  void set_cyl_constraint(double center_x, double center_y, double radius);
  void set_slab_constraint(double slab_start_z, double slab_end_z);
  Utils::Vector2d get_slab_constraint_parameters() const {
    return {m_slab_start_z, m_slab_end_z};
  }

  int delete_particle(int p_id);
  void add_reaction(double gamma, const std::vector<int> &reactant_types,
                    const std::vector<int> &reactant_coefficients,
                    const std::vector<int> &product_types,
                    const std::vector<int> &product_coefficients);
  void delete_reaction(int reaction_id) {
    reactions.erase(reactions.begin() + reaction_id);
  }

  bool do_global_mc_move_for_particles_of_type(int type,
                                               int particle_number_of_type);

  bool particle_inside_exclusion_radius_touched = false;

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

  std::tuple<std::vector<StoredParticleProperty>, std::vector<int>,
             std::vector<StoredParticleProperty>>
  make_reaction_attempt(SingleReaction const &current_reaction);
  std::vector<std::pair<int, Utils::Vector3d>>
  generate_new_particle_positions(int type, int n_particles);
  void
  restore_properties(std::vector<StoredParticleProperty> const &property_list,
                     int number_of_saved_properties);

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

protected:
  virtual double calculate_acceptance_probability(
      SingleReaction const &current_reaction, double E_pot_old,
      double E_pot_new, std::map<int, int> const &old_particle_numbers) const {
    return -10;
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
  void check_exclusion_radius(int p_id);
  void move_particle(int p_id, Utils::Vector3d const &new_pos,
                     double velocity_prefactor);

  void append_particle_property_of_random_particle(
      int type, std::vector<StoredParticleProperty> &list_of_particles);

  enum class ReactionConstraint { NONE, CYL_Z, SLAB_Z };
  ReactionConstraint m_reaction_constraint = ReactionConstraint::NONE;
  double m_cyl_radius = -10.0;
  double m_cyl_x = -10.0;
  double m_cyl_y = -10.0;
  double m_slab_start_z = -10.0;
  double m_slab_end_z = -10.0;

protected:
  Utils::Vector3d get_random_position_in_box();
};

} // namespace ReactionMethods
#endif
