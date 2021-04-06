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

#include "reaction_methods/SingleReaction.hpp"

#include "random.hpp"

#include <utils/Vector.hpp>

#include <map>
#include <random>
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
  double temperature = -10.0;
  double exclusion_radius =
      0.0; // this is used as a kind of hard sphere radius, if
           // particles are closer than that it is assumed that
           // their interaction energy gets approximately
           // infinite => these configurations do not contribute
           // to the partition function and ensemble averages.
  double volume = -10.0;
  bool box_is_cylindric_around_z_axis = false;
  double cyl_radius = -10.0;
  double cyl_x = -10.0;
  double cyl_y = -10.0;
  bool box_has_wall_constraints = false;
  double slab_start_z = -10.0;
  double slab_end_z = -10.0;
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

  int delete_particle(int p_id);
  void add_reaction(double gamma, const std::vector<int> &reactant_types,
                    const std::vector<int> &reactant_coefficients,
                    const std::vector<int> &product_types,
                    const std::vector<int> &product_coefficients);
  void delete_reaction(int reaction_id) {
    reactions.erase(reactions.begin() + reaction_id);
  }

  bool do_global_mc_move_for_particles_of_type(int type,
                                               int particle_number_of_type,
                                               bool use_wang_landau);

  bool particle_inside_exclusion_radius_touched;

protected:
  std::vector<int> m_empty_p_ids_smaller_than_max_seen_particle;
  void generic_oneway_reaction(int reaction_id);
  virtual void on_reaction_entry(int &old_state_index) {}
  virtual void
  on_reaction_rejection_directly_after_entry(int &old_state_index) {}
  virtual void on_attempted_reaction(int &new_state_index) {}
  virtual void on_end_reaction(int &accepted_state) {}

  virtual void on_mc_rejection_directly_after_entry(int &old_state_index){};
  virtual void on_mc_accept(int &new_state_index){};
  virtual void on_mc_reject(int &old_state_index){};
  virtual int on_mc_use_WL_get_new_state() { return -10; }

  void make_reaction_attempt(
      SingleReaction const &current_reaction,
      std::vector<StoredParticleProperty> &changed_particles_properties,
      std::vector<int> &p_ids_created_particles,
      std::vector<StoredParticleProperty> &hidden_particles_properties);
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
  bool all_reactant_particles_exist(int reaction_id) const;

protected:
  virtual double calculate_acceptance_probability(
      SingleReaction const &current_reaction, double E_pot_old,
      double E_pot_new, std::map<int, int> const &old_particle_numbers,
      int old_state_index, int new_state_index,
      bool only_make_configuration_changing_move) const {
    return -10;
  }

private:
  std::mt19937 m_generator;
  std::normal_distribution<double> m_normal_distribution;
  std::uniform_real_distribution<double> m_uniform_real_distribution;

  std::map<int, int> save_old_particle_numbers(int reaction_id);

  void replace_particle(int p_id, int desired_type);
  int create_particle(int desired_type);
  void hide_particle(int p_id, int previous_type);

  void append_particle_property_of_random_particle(
      int type, std::vector<StoredParticleProperty> &list_of_particles);

  Utils::Vector3d get_random_position_in_box();
};

} // namespace ReactionMethods
#endif
