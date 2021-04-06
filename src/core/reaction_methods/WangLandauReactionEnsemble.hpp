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
#ifndef REACTION_METHODS_WANG_LANDAU_REACTION_ENSEMBLE_HPP
#define REACTION_METHODS_WANG_LANDAU_REACTION_ENSEMBLE_HPP

#include "reaction_methods/ReactionAlgorithm.hpp"

#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace ReactionMethods {

struct CollectiveVariable {
  double CV_minimum = {};
  double CV_maximum = {};
  double delta_CV = {};
  // use pure virtual, otherwise this will be used in vector of collective
  // variables
  virtual double determine_current_state() const = 0;
  virtual ~CollectiveVariable() = default;
};

class WangLandauReactionEnsemble;

struct EnergyCollectiveVariable : public CollectiveVariable {
  double determine_current_state() const override;
  void load_CV_boundaries(WangLandauReactionEnsemble &wl_system,
                          std::istream &infile);
};

/** Measure the degree of association of a chemical species.
 *  As an example, consider a polybasic acid A which can be protonated
 *  into species AH and AH2. The degree of association of species A is
 *  equal to n(A) / (n(A) + n(AH) + n(AH2)).
 */
struct DegreeOfAssociationCollectiveVariable : public CollectiveVariable {
  /** List of all conjugated species */
  std::vector<int> corresponding_acid_types;
  /** Reference species from which the degree of association is measured */
  int associated_type;
  /** Calculate the degree of association of the reference species */
  double determine_current_state() const override {
    return calculate_degree_of_association();
  }

private:
  /**
   * Return the degree of association for the current collective variable. This
   * is needed since you may use multiple degrees of association as collective
   * variable for the Wang-Landau algorithm.
   */
  double calculate_degree_of_association() const;
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
  void add_new_CV_potential_energy(std::istream &infile, double delta_CV);
  std::vector<std::shared_ptr<CollectiveVariable>> collective_variables;

  std::string output_filename = "";

  std::vector<double> min_boundaries_energies;
  std::vector<double> max_boundaries_energies;

  // only present in energy preparation run
  std::vector<double> minimum_energies_at_flat_index;
  // only present in energy preparation run
  std::vector<double> maximum_energies_at_flat_index;

  // use for preliminary energy reweighting runs
  int update_maximum_and_minimum_energies_at_current_state();
  int do_reaction(int reaction_steps) override;
  void write_out_preliminary_energy_run_results(const std::string &filename);

  // checkpointing, only designed to reassign values of a previous simulation to
  // a new simulation with the same initialization process
  void load_wang_landau_checkpoint(const std::string &identifier);
  void write_wang_landau_checkpoint(const std::string &identifier);
  void write_wang_landau_results_to_file(const std::string &filename);

protected:
  double calculate_acceptance_probability(
      SingleReaction const &current_reaction, double E_pot_old,
      double E_pot_new, std::map<int, int> const &old_particle_numbers,
      int old_state_index, int new_state_index,
      bool only_make_configuration_changing_move) const override;
  void format_wang_landau_results(std::ostream &out);

private:
  void on_reaction_entry(int &old_state_index) override;
  void
  on_reaction_rejection_directly_after_entry(int &old_state_index) override;
  void on_attempted_reaction(int &new_state_index) override;
  void on_end_reaction(int &accepted_state) override;
  void on_mc_rejection_directly_after_entry(int &old_state_index) override;
  void on_mc_accept(int &new_state_index) override;
  void on_mc_reject(int &old_state_index) override;
  int on_mc_use_WL_get_new_state() override;

  std::vector<int> histogram;

  // logarithm to basis e of the degeneracy of the states
  std::vector<double> wang_landau_potential;

  std::vector<int> nr_subindices_of_collective_variable;

  // logarithm to basis e of the modification factor of the degeneracy of
  // states when the state is visited
  double wang_landau_parameter = 1.0;

  int int_fill_value = -10;
  double double_fill_value = -10.0;

  long used_bins = -10;            // for 1/t algorithm
  int monte_carlo_trial_moves = 0; // for 1/t algorithm

  int get_flattened_index_wang_landau_without_energy_collective_variable(
      int flattened_index_with_EnergyCollectiveVariable,
      int collective_variable_index_energy_observable); // needed for energy

  int get_flattened_index_wang_landau(
      std::vector<double> const &current_state,
      std::vector<double> const &collective_variables_minimum_values,
      std::vector<double> const &collective_variables_maximum_values,
      std::vector<double> const &delta_collective_variables_values,
      int nr_collective_variables); // collective variable
  int get_flattened_index_wang_landau_of_current_state();

  void update_wang_landau_potential_and_histogram(
      int index_of_state_after_acceptance_or_rejection);
  int m_WL_tries = 0;
  bool can_refine_wang_landau_one_over_t() const;
  bool m_system_is_in_1_over_t_regime = false;
  bool achieved_desired_number_of_refinements_one_over_t() const;
  void refine_wang_landau_parameter_one_over_t();

  /**
   * Initialize the current Wang-Landau system.
   * Has to be called (at least) after the last collective variable is added.
   */
  int initialize_wang_landau();
  double calculate_delta_degree_of_association(
      DegreeOfAssociationCollectiveVariable &current_collective_variable);
  void invalidate_bins();
  void reset_histogram();
  double get_minimum_CV_value_on_delta_CV_spaced_grid(double min_CV_value,
                                                      double delta_CV) const;
};

} // namespace ReactionMethods
#endif
