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
#ifndef REACTION_ENSEMBLE_H
#define REACTION_ENSEMBLE_H

#include "energy.hpp"
#include "particle_data.hpp"
#include "random.hpp"

#include <map>
#include <string>
#include <utils/Accumulator.hpp>

namespace ReactionEnsemble {

struct SingleReaction {
  // strict input to the algorithm
  std::vector<int> reactant_types;
  std::vector<int> reactant_coefficients;
  std::vector<int> product_types;
  std::vector<int> product_coefficients;
  double gamma = {};
  // calculated values that are stored for performance reasons
  int nu_bar = {};
  Utils::Accumulator accumulator_exponentials = Utils::Accumulator(1);
  int tried_moves = 0;
  int accepted_moves = 0;
  double get_acceptance_rate() {
    return static_cast<double>(accepted_moves) /
           static_cast<double>(tried_moves);
  };
};

struct StoredParticleProperty {
  int p_id;
  double charge;
  int type;
};

struct CollectiveVariable {
  double CV_minimum = {};
  double CV_maximum = {};
  double delta_CV = {};
  virtual double determine_current_state() = 0; // use pure virtual, otherwise
                                                // this will be used in vector
                                                // of collective variables
  virtual ~CollectiveVariable() = default;
};

class WangLandauReactionEnsemble;

struct EnergyCollectiveVariable : public CollectiveVariable {
  std::string energy_boundaries_filename;
  double determine_current_state() override {
    return calculate_current_potential_energy_of_system();
  }
  void
  load_CV_boundaries(WangLandauReactionEnsemble &m_current_wang_landau_system);
};

struct DegreeOfAssociationCollectiveVariable : public CollectiveVariable {
  std::vector<int> corresponding_acid_types;
  int associated_type;
  double determine_current_state() override {
    return calculate_degree_of_association();
  }

private:
  /**
   * Returns the degree of association for the current collective variable. This
   * is needed since you may use multiple degrees of association as collective
   * variable for the Wang-Landau algorithm.
   */
  double calculate_degree_of_association() {
    int total_number_of_corresponding_acid = 0;
    for (int corresponding_acid_type : corresponding_acid_types) {
      int num_of_current_type =
          number_of_particles_with_type(corresponding_acid_type);
      total_number_of_corresponding_acid += num_of_current_type;
    }
    if (total_number_of_corresponding_acid == 0) {
      throw std::runtime_error("Have you forgotten to specify all "
                               "corresponding acid types? Total particle "
                               "number of corresponding acid type is zero\n");
    }
    int num_of_associated_acid = number_of_particles_with_type(associated_type);
    double degree_of_association =
        static_cast<double>(num_of_associated_acid) /
        total_number_of_corresponding_acid; // cast to double because otherwise
                                            // any fractional part is lost
    return degree_of_association;
  }
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
  double get_acceptance_rate_configurational_moves() {
    return static_cast<double>(m_accepted_configurational_MC_moves) /
           static_cast<double>(m_tried_configurational_MC_moves);
  }

  void set_cuboid_reaction_ensemble_volume();
  virtual int do_reaction(int reaction_steps);
  void check_reaction_ensemble();

  int delete_particle(int p_id);
  void add_reaction(double gamma, const std::vector<int> &_reactant_types,
                    const std::vector<int> &_reactant_coefficients,
                    const std::vector<int> &_product_types,
                    const std::vector<int> &_product_coefficients);
  void delete_reaction(int reaction_id) {
    reactions.erase(reactions.begin() + reaction_id);
  }

  bool do_global_mc_move_for_particles_of_type(int type,
                                               int particle_number_of_type,
                                               bool use_wang_landau);

  bool particle_inside_exclusion_radius_touched;

protected:
  std::vector<int> m_empty_p_ids_smaller_than_max_seen_particle;
  bool generic_oneway_reaction(int reaction_id);
  virtual void on_reaction_entry(int &old_state_index) {}
  virtual void
  on_reaction_rejection_directly_after_entry(int &old_state_index) {}
  virtual void on_attempted_reaction(int &new_state_index) {}
  virtual void on_end_reaction(int &accepted_state) {}

  virtual void on_mc_rejection_directly_after_entry(int &old_state_index){};
  virtual void on_mc_accept(int &new_state_index){};
  virtual void on_mc_reject(int &old_state_index){};
  virtual int on_mc_use_WL_get_new_state() { return -10; };

  void make_reaction_attempt(
      SingleReaction &current_reaction,
      std::vector<StoredParticleProperty> &changed_particles_properties,
      std::vector<int> &p_ids_created_particles,
      std::vector<StoredParticleProperty> &hidden_particles_properties);
  void restore_properties(std::vector<StoredParticleProperty> &property_list,
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
  bool all_reactant_particles_exist(int reaction_id);

private:
  std::mt19937 m_generator;
  std::normal_distribution<double> m_normal_distribution;
  std::uniform_real_distribution<double> m_uniform_real_distribution;

  std::map<int, int> save_old_particle_numbers(int reaction_id);

  int calculate_nu_bar(
      std::vector<int> &reactant_coefficients,
      std::vector<int> &product_coefficients); // should only be used when
                                               // defining a new reaction
  int m_invalid_charge =
      -10000; // this is the default charge which is assigned to a type which
              // occurs in a reaction. this charge has to be overwritten. if it
              // is not overwritten the reaction ensemble will complain.
  void replace_particle(int p_id, int desired_type);
  int create_particle(int desired_type);
  void hide_particle(int p_id, int previous_type);

  void append_particle_property_of_random_particle(
      int type, std::vector<StoredParticleProperty> &list_of_particles);

  virtual double calculate_acceptance_probability(
      SingleReaction &current_reaction, double E_pot_old, double E_pot_new,
      std::map<int, int> &old_particle_numbers, int old_state_index,
      int new_state_index, bool only_make_configuration_changing_move) {
    return -10;
  };

  void add_types_to_index(std::vector<int> &type_list);
  Utils::Vector3d get_random_position_in_box();
  std::vector<double>
  get_random_position_in_box_enhanced_proposal_of_small_radii();
};

///////////////////////////// actual declaration of specific reaction algorithms

/** Reaction ensemble method.
 *  Works for the reaction ensemble at constant volume and temperature. For the
 *  reaction ensemble at constant pressure, additionally employ a barostat!
 *  NOTE: a chemical reaction consists of a forward and backward reaction.
 *  Here both reactions have to be defined separately. The extent of the
 *  reaction is here chosen to be +1. If the reaction trial move for a
 *  dissociation of HA is accepted then there is one more dissociated ion
 *  pair H+ and A-. Implementation of @cite smith94c.
 */
class ReactionEnsemble : public ReactionAlgorithm {
public:
  ReactionEnsemble(int seed) : ReactionAlgorithm(seed) {}

private:
  double calculate_acceptance_probability(
      SingleReaction &current_reaction, double E_pot_old, double E_pot_new,
      std::map<int, int> &old_particle_numbers, int dummy_old_state_index,
      int dummy_new_state_index,
      bool dummy_only_make_configuration_changing_move) override;
};

/** Wang-Landau reaction ensemble method */
class WangLandauReactionEnsemble : public ReactionAlgorithm {
public:
  WangLandauReactionEnsemble(int seed) : ReactionAlgorithm(seed) {}
  bool do_energy_reweighting = false;
  bool do_not_sample_reaction_partition_function = false;
  double final_wang_landau_parameter = 0.00001;

  void add_new_CV_degree_of_association(
      int associated_type, double CV_minimum, double CV_maximum,
      const std::vector<int> &_corresponding_acid_types);
  void add_new_CV_potential_energy(const std::string &filename,
                                   double delta_CV);
  std::vector<std::shared_ptr<CollectiveVariable>> collective_variables;

  std::string output_filename = "";

  std::vector<double> min_boundaries_energies;
  std::vector<double> max_boundaries_energies;

  std::vector<double>
      minimum_energies_at_flat_index; // only present in energy preparation run
  std::vector<double>
      maximum_energies_at_flat_index; // only present in energy preparation run

  int update_maximum_and_minimum_energies_at_current_state(); // use for
                                                              // preliminary
                                                              // energy
                                                              // reweighting
                                                              // runs
  int do_reaction(int reaction_steps) override;
  void write_out_preliminary_energy_run_results(const std::string &filename);

  // checkpointing, only designed to reassign values of a previous simulation to
  // a new simulation with the same initialization process
  int load_wang_landau_checkpoint(const std::string &identifier);
  int write_wang_landau_checkpoint(const std::string &identifier);
  void write_wang_landau_results_to_file(
      const std::string &full_path_to_output_filename);

private:
  void on_reaction_entry(int &old_state_index) override;
  void
  on_reaction_rejection_directly_after_entry(int &old_state_index) override;
  void on_attempted_reaction(int &new_state_index) override;
  void on_end_reaction(int &accepted_state) override;
  double calculate_acceptance_probability(
      SingleReaction &current_reaction, double E_pot_old, double E_pot_new,
      std::map<int, int> &old_particle_numbers, int old_state_index,
      int new_state_index, bool only_make_configuration_changing_move) override;
  void on_mc_rejection_directly_after_entry(int &old_state_index) override;
  void on_mc_accept(int &new_state_index) override;
  void on_mc_reject(int &old_state_index) override;
  int on_mc_use_WL_get_new_state() override;

  std::vector<int> histogram;
  std::vector<double> wang_landau_potential; // equals the logarithm to basis e
                                             // of the degeneracy of the states

  std::vector<int> nr_subindices_of_collective_variable;
  double wang_landau_parameter = 1.0; // equals the logarithm to basis e of the
  // modification factor of the degeneracy of
  // states when the state is visited
  double initial_wang_landau_parameter = 1.0;

  int int_fill_value = -10;
  double double_fill_value = -10.0;

  int used_bins = -10;             // for 1/t algorithm
  int monte_carlo_trial_moves = 0; // for 1/t algorithm

  int get_flattened_index_wang_landau_without_energy_collective_variable(
      int flattened_index_with_EnergyCollectiveVariable,
      int collective_variable_index_energy_observable); // needed for energy

  int get_flattened_index_wang_landau(
      std::vector<double> &current_state,
      std::vector<double> &collective_variables_minimum_values,
      std::vector<double> &collective_variables_maximum_values,
      std::vector<double> &delta_collective_variables_values,
      int nr_collective_variables); // collective variable
  int get_flattened_index_wang_landau_of_current_state();

  void update_wang_landau_potential_and_histogram(
      int index_of_state_after_acceptance_or_rejection);
  int m_WL_accepted_moves = 0;
  int m_WL_tries = 0;
  bool can_refine_wang_landau_one_over_t();
  bool m_system_is_in_1_over_t_regime = false;
  bool achieved_desired_number_of_refinements_one_over_t();
  void refine_wang_landau_parameter_one_over_t();

  int initialize_wang_landau(); // has to be called (at least) after the last
                                // collective variable is added
  double calculate_delta_degree_of_association(
      DegreeOfAssociationCollectiveVariable &current_collective_variable);
  int get_num_needed_bins();
  void invalidate_bins();
  void remove_bins_that_have_not_been_sampled();
  void reset_histogram();
  double get_minimum_CV_value_on_delta_CV_spaced_grid(double min_CV_value,
                                                      double delta_CV);
};

/**
 * Constant-pH Ensemble, for derivation see @cite reed92a.
 * For the constant pH reactions you need to provide the deprotonation and
 * afterwards the corresponding protonation reaction (in this order). If you
 * want to deal with multiple reactions do it multiple times. Note that there is
 * a difference in the usecase of the constant pH reactions and the above
 * reaction ensemble. For the constant pH simulation directily the
 * **apparent equilibrium constant which carries a unit** needs to be provided
 * -- this is equivalent to the gamma of the reaction ensemble above, where the
 * dimensionless reaction constant needs to be provided. Again: For the
 * constant-pH algorithm not the dimensionless reaction constant needs to be
 * provided here, but the apparent reaction constant.
 */
class ConstantpHEnsemble : public ReactionAlgorithm {
public:
  ConstantpHEnsemble(int seed) : ReactionAlgorithm(seed) {}
  double m_constant_pH = -10;
  int do_reaction(int reaction_steps) override;

private:
  double calculate_acceptance_probability(
      SingleReaction &current_reaction, double E_pot_old, double E_pot_new,
      std::map<int, int> &dummy_old_particle_numbers, int dummy_old_state_index,
      int dummy_new_state_index,
      bool dummy_only_make_configuration_changing_move) override;
  int get_random_valid_p_id();
};

/** Widom insertion method */
class WidomInsertion : public ReactionAlgorithm {
public:
  WidomInsertion(int seed) : ReactionAlgorithm(seed) {}
  std::pair<double, double> measure_excess_chemical_potential(int reaction_id);
};

//////////////////////////////////////////////////////////////////free functions
double calculate_factorial_expression(SingleReaction &current_reaction,
                                      std::map<int, int> &old_particle_numbers);
double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i);
} // namespace ReactionEnsemble
#endif
