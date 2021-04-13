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

#include "reaction_methods/WangLandauReactionEnsemble.hpp"

#include "utils.hpp"

#include "energy.hpp"
#include "particle_data.hpp"

#include <boost/range/algorithm.hpp>

#include <utils/index.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ReactionMethods {

double EnergyCollectiveVariable::determine_current_state() const {
  return calculate_current_potential_energy_of_system();
}

double
DegreeOfAssociationCollectiveVariable::calculate_degree_of_association() const {
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
      static_cast<double>(total_number_of_corresponding_acid);
  return degree_of_association;
}

/** Load minimum and maximum energies as a function of the other collective
 *  variables.
 */
void EnergyCollectiveVariable::load_CV_boundaries(
    WangLandauReactionEnsemble &wl_system, std::istream &infile) {

  wl_system.do_energy_reweighting = true;

  // Note that you cannot change the other collective variables in the
  // pre-production run and the production run

  // The data is formatted as space-separated floating point values
  // (the first line is a header). The min and max energies are stored in
  // the last two columns. The first N columns are the collective variables.
  std::string line = "";
  std::getline(infile, line); // dummy read to throw away header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::vector<double> values;
    double value = -1.0;
    while (iss >> value) {
      values.push_back(value);
    }
    assert(values.size() >= 2);
    wl_system.max_boundaries_energies.emplace_back(values.back());
    wl_system.min_boundaries_energies.emplace_back(values[values.size() - 2]);
  }

  CV_minimum = *boost::range::min_element(wl_system.min_boundaries_energies);
  CV_maximum = *boost::range::max_element(wl_system.max_boundaries_energies);
}

void WangLandauReactionEnsemble::on_reaction_entry(int &old_state_index) {
  old_state_index = get_flattened_index_wang_landau_of_current_state();
  if (old_state_index >= 0) {
    if (histogram[old_state_index] >= 0)
      monte_carlo_trial_moves += 1;
  }
}

void WangLandauReactionEnsemble::on_reaction_rejection_directly_after_entry(
    int &old_state_index) {
  // increase the Wang-Landau potential and histogram at the current nbar
  // (this case covers the cases nbar=0 or nbar=1)
  update_wang_landau_potential_and_histogram(old_state_index);
}

void WangLandauReactionEnsemble::on_attempted_reaction(int &new_state_index) {
  new_state_index = get_flattened_index_wang_landau_of_current_state();
}

void WangLandauReactionEnsemble::on_end_reaction(int &accepted_state) {
  update_wang_landau_potential_and_histogram(accepted_state);
}

void WangLandauReactionEnsemble::on_mc_rejection_directly_after_entry(
    int &old_state_index) {
  if (do_energy_reweighting)
    update_wang_landau_potential_and_histogram(old_state_index);
}

void WangLandauReactionEnsemble::on_mc_accept(int &new_state_index) {
  if (do_energy_reweighting) {
    // modify wang_landau histogram and potential
    update_wang_landau_potential_and_histogram(new_state_index);
  }
}

void WangLandauReactionEnsemble::on_mc_reject(int &old_state_index) {
  if (do_energy_reweighting)
    update_wang_landau_potential_and_histogram(old_state_index);
}

int WangLandauReactionEnsemble::on_mc_use_WL_get_new_state() {
  return get_flattened_index_wang_landau_of_current_state();
}

/**
 * Adds a new collective variable (CV) of the type degree of association to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_degree_of_association(
    int associated_type, double CV_minimum, double CV_maximum,
    const std::vector<int> &corresponding_acid_types) {
  std::shared_ptr<DegreeOfAssociationCollectiveVariable>
      new_collective_variable =
          std::make_shared<DegreeOfAssociationCollectiveVariable>();
  new_collective_variable->associated_type = associated_type;
  new_collective_variable->CV_minimum = CV_minimum;
  new_collective_variable->CV_maximum = CV_maximum;
  new_collective_variable->corresponding_acid_types = corresponding_acid_types;
  new_collective_variable->delta_CV =
      calculate_delta_degree_of_association(*new_collective_variable);
  collective_variables.push_back(new_collective_variable);
  initialize_wang_landau();
}

/**
 * Adds a new collective variable (CV) of the type potential energy to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_potential_energy(
    std::istream &infile, double delta_CV) {
  std::shared_ptr<EnergyCollectiveVariable> new_collective_variable =
      std::make_shared<EnergyCollectiveVariable>();
  new_collective_variable->delta_CV = delta_CV;
  new_collective_variable->load_CV_boundaries(*this, infile);
  collective_variables.emplace_back(new_collective_variable);
  initialize_wang_landau();
}

/**
 * Adds a new collective variable (CV) of the type potential energy to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_potential_energy(
    const std::string &filename, double delta_CV) {
  std::ifstream infile;
  infile.open(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Cannot read " + filename);
  }
  add_new_CV_potential_energy(infile, delta_CV);
  infile.close();
}

/**
 * Returns the flattened index of the multidimensional Wang-Landau histogram
 */
int WangLandauReactionEnsemble::get_flattened_index_wang_landau(
    std::vector<double> const &current_state,
    std::vector<double> const &collective_variables_minimum_values,
    std::vector<double> const &collective_variables_maximum_values,
    std::vector<double> const &delta_collective_variables_values,
    int nr_collective_variables) {

  // check for the current state to be an allowed state in the range
  // [collective_variables_minimum_values:collective_variables_maximum_values],
  // else return a negative index (to indicate error)
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    if (current_state[CV_i] >
            collective_variables_maximum_values[CV_i] +
                delta_collective_variables_values[CV_i] * 0.98 ||
        current_state[CV_i] <
            collective_variables_minimum_values[CV_i] -
                delta_collective_variables_values[CV_i] * 0.01) {
      return -10;
    }
  }

  std::vector<int> individual_indices(nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    auto const value =
        (current_state[CV_i] - collective_variables_minimum_values[CV_i]) /
        delta_collective_variables_values[CV_i];
    int rounded_value;
    if (CV_i == collective_variables.size() - 1 && do_energy_reweighting) {
      // for energy collective variable (simple truncating conversion desired)
      rounded_value = static_cast<int>(value);
    } else {
      // for degree of association collective variables (rounding conversion)
      rounded_value = static_cast<int>(std::floor(value));
    }
    if (rounded_value < 0 or
        rounded_value >= nr_subindices_of_collective_variable[CV_i]) {
      return -10;
    }
    individual_indices[CV_i] = rounded_value;
  }
  // get flattened index from individual_indices
  int index = 0;
  // this is already part of the algorithm to find the correct index
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    int factor = 1;
    for (int j = CV_i + 1; j < nr_collective_variables; j++) {
      factor *= nr_subindices_of_collective_variable[j];
    }
    index += factor * individual_indices[CV_i];
  }
  return index;
}

/**
 * Returns the flattened index of the multidimensional Wang-Landau histogram for
 * the current state of the simulation.
 */
int WangLandauReactionEnsemble::
    get_flattened_index_wang_landau_of_current_state() {
  auto const nr_collective_variables =
      static_cast<int>(collective_variables.size());
  // get current state
  std::vector<double> current_state(nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++)
    current_state[CV_i] = collective_variables[CV_i]->determine_current_state();
  // get collective_variables_minimum_values
  std::vector<double> collective_variables_minimum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_minimum_values[CV_i] =
        collective_variables[CV_i]->CV_minimum;
  }
  // get collective_variables_maximum_values
  std::vector<double> collective_variables_maximum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_maximum_values[CV_i] =
        collective_variables[CV_i]->CV_maximum;
  }
  // get delta_collective_variables_values
  std::vector<double> delta_collective_variables_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    delta_collective_variables_values[CV_i] =
        collective_variables[CV_i]->delta_CV;
  }
  int index = get_flattened_index_wang_landau(
      current_state, collective_variables_minimum_values,
      collective_variables_maximum_values, delta_collective_variables_values,
      nr_collective_variables);
  return index;
}

/**
 * Returns the minimum value of the collective variable on a delta_CV spaced
 * grid which starts at 0
 */
double WangLandauReactionEnsemble::get_minimum_CV_value_on_delta_CV_spaced_grid(
    double min_CV_value, double delta_CV) const {
  // assume grid has it s origin at 0
  double minimum_CV_value_on_delta_CV_spaced_grid =
      floor(min_CV_value / delta_CV) * delta_CV;
  return minimum_CV_value_on_delta_CV_spaced_grid;
}

/**
 * Calculates the smallest difference in the degree of association which can be
 * observed when changing the degree of association by one single reaction.
 */
double WangLandauReactionEnsemble::calculate_delta_degree_of_association(
    DegreeOfAssociationCollectiveVariable &current_collective_variable) {
  // calculate Delta in the degree of association so that EVERY reaction step is
  // driven.
  int total_number_of_corresponding_acid = 0;
  for (int corresponding_acid_type :
       current_collective_variable.corresponding_acid_types) {
    int num_of_current_type =
        number_of_particles_with_type(corresponding_acid_type);
    total_number_of_corresponding_acid += num_of_current_type;
  }
  double delta = 1.0 / total_number_of_corresponding_acid;
  // now modify the minimum value of the CV to lie on the grid
  current_collective_variable.CV_minimum =
      get_minimum_CV_value_on_delta_CV_spaced_grid(
          current_collective_variable.CV_minimum, delta);
  return delta;
}

void WangLandauReactionEnsemble::invalidate_bins() {
  // make values in histogram and Wang-Landau potential negative if they are not
  // allowed at the given degree of association, because the energy boundaries
  // prohibit them

  long empty_bins_in_memory = 0;
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    std::vector<std::size_t> unraveled_index(
        nr_subindices_of_collective_variable.size());
    Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                         nr_subindices_of_collective_variable.end(),
                         unraveled_index.begin(), flattened_index);
    // use unraveled index
    int EnergyCollectiveVariable_index = 0;
    if (collective_variables.size() > 1) {
      // assume the energy collective variable to be the last added
      // collective variable
      EnergyCollectiveVariable_index =
          static_cast<int>(collective_variables.size()) - 1;
    }
    double current_energy =
        static_cast<double>(unraveled_index[EnergyCollectiveVariable_index]) *
            collective_variables[EnergyCollectiveVariable_index]->delta_CV +
        collective_variables[EnergyCollectiveVariable_index]->CV_minimum;
    int flat_index_without_energy_CV =
        get_flattened_index_wang_landau_without_energy_collective_variable(
            static_cast<int>(flattened_index), EnergyCollectiveVariable_index);
    std::shared_ptr<CollectiveVariable> energy_CV =
        collective_variables[EnergyCollectiveVariable_index];
    if (current_energy >
            max_boundaries_energies[flat_index_without_energy_CV] ||
        current_energy < min_boundaries_energies[flat_index_without_energy_CV] -
                             energy_CV->delta_CV) {
      histogram[flattened_index] = int_fill_value;
      wang_landau_potential[flattened_index] = double_fill_value;
      empty_bins_in_memory++;
    }
  }

  used_bins =
      static_cast<long>(wang_landau_potential.size()) - empty_bins_in_memory;
}

int WangLandauReactionEnsemble::initialize_wang_landau() {

  nr_subindices_of_collective_variable.resize(collective_variables.size(), 0);
  auto const &cv = collective_variables.back();
  // add 1 for collective variables which are of type degree of association
  nr_subindices_of_collective_variable.back() =
      static_cast<int>((cv->CV_maximum - cv->CV_minimum) / cv->delta_CV) + 1;

  // construct (possibly higher dimensional) histogram and potential over
  // gamma (the space which should be equally sampled when the Wang-Landau
  // algorithm has converged)
  // add 1 for degrees of association-related part of histogram (think of only
  // one acid particle)
  auto const needed_bins = std::accumulate(
      collective_variables.begin(), collective_variables.end(), 1l,
      [](long acc, std::shared_ptr<CollectiveVariable> const &cv) {
        return acc * (static_cast<long>((cv->CV_maximum - cv->CV_minimum) /
                                        cv->delta_CV) +
                      1l);
      });
  assert(needed_bins >= 0);
  histogram.resize(needed_bins, 0);
  wang_landau_potential.resize(needed_bins, 0);

  used_bins = needed_bins; // initialize for 1/t wang_landau algorithm

  if (do_energy_reweighting) {
    invalidate_bins();
  }
  return ES_OK;
}

/** Calculate the expression in the acceptance probability of the Wang-Landau
 *  reaction ensemble.
 *  Modify Boltzmann factor according to Wang-Landau algorithm in @cite yan02b.
 */
double WangLandauReactionEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &old_particle_numbers, int old_state_index,
    int new_state_index, bool only_make_configuration_changing_move) const {

  double beta = 1.0 / temperature;
  double bf = 1.0;

  if (!(do_not_sample_reaction_partition_function ||
        only_make_configuration_changing_move)) {
    auto const factorial_expr =
        calculate_factorial_expression(current_reaction, old_particle_numbers);
    bf = std::pow(volume, current_reaction.nu_bar) * current_reaction.gamma *
         factorial_expr;
  }

  if (!do_energy_reweighting) {
    bf *= exp(-beta * (E_pot_new - E_pot_old));
  }

  // Check whether the proposed state lies in the reaction coordinate space
  // gamma and add the Wang-Landau modification factor, this is a bit nasty
  // due to the energy collective variable case (memory layout of storage
  // array of the histogram and the wang_landau_potential values is "cuboid").
  if (old_state_index >= 0 && new_state_index >= 0) {
    if (histogram[new_state_index] >= 0 && histogram[old_state_index] >= 0) {
      // Modify Boltzmann factor (bf) according to Wang-Landau algorithm
      // @cite yan02b.
      // This makes the new state being accepted with the conditional
      // probability bf (bf is a transition probability = conditional
      // probability from the old state to move to the new state).
      bf = std::min(1.0, bf * exp(wang_landau_potential[old_state_index] -
                                  wang_landau_potential[new_state_index]));
    } else {
      if (histogram[old_state_index] < 0)
        bf = 10; // accept the reaction if we found a state in gamma
                 // (histogram[new_state_index] >= 0) or to sample new configs
                 // which might lie in gamma(histogram[new_state_index] < 0)
      else if (histogram[new_state_index] < 0)
        bf = -10; // reject the reaction, since the new state
                  // is not in gamma while the old sate was in gamma
    }
  } else if (old_state_index < 0 && new_state_index >= 0) {
    bf = 10; // accept the reaction if we found a new state in gamma
             // (new_state_index >= 0) or to sample new configs which
             // might lie in gamma (new_state_index < 0)
  } else if (old_state_index >= 0 && new_state_index < 0) {
    bf = -10; // reject the reaction, since the new state is
              // not in gamma while the old sate was in gamma
  }
  return bf;
}

/** Perform a randomly selected reaction using the Wang-Landau algorithm.
 *
 *  Make sure to perform additional configuration changing steps, after the
 *  reaction step! Like in @cite yan02b. This can be done with MD in the case
 *  of the no-energy-reweighting case, or with the function
 *  @ref ReactionAlgorithm::do_global_mc_move_for_particles_of_type.
 *
 *  Perform additional Monte Carlo moves to sample configurational
 *  partition function according to @cite yan02b. Do as many steps
 *  as needed to get to a new conformation.
 */
int WangLandauReactionEnsemble::do_reaction(int reaction_steps) {
  m_WL_tries += reaction_steps;
  for (int step = 0; step < reaction_steps; step++) {
    int reaction_id = i_random(static_cast<int>(reactions.size()));
    generic_oneway_reaction(reaction_id);
    if (can_refine_wang_landau_one_over_t() && m_WL_tries % 10000 == 0) {
      // check for convergence
      if (achieved_desired_number_of_refinements_one_over_t()) {
        write_wang_landau_results_to_file(output_filename);
        return -10; // return negative value to indicate that the Wang-Landau
                    // algorithm has converged
      }
      refine_wang_landau_parameter_one_over_t();
    }
  }
  // shift Wang-Landau potential minimum to zero
  if (m_WL_tries % (std::max(90000, 9 * reaction_steps)) == 0) {
    // for numerical stability here we also subtract the minimum positive value
    // of the wang_landau_potential from the wang_landau potential, allowed
    // since only the difference in the wang_landau potential is of interest.
    double minimum_wang_landau_potential =
        find_minimum_non_negative_value(wang_landau_potential);
    for (double &i : wang_landau_potential) {
      if (i >= 0) // check for whether we are in the
                  // valid range of the collective variable
        i -= minimum_wang_landau_potential;
    }
    // write out preliminary Wang-Landau potential results
    write_wang_landau_results_to_file(output_filename);
  }
  return 0;
}

/** Increase the Wang-Landau potential and histogram at the current nbar */
void WangLandauReactionEnsemble::update_wang_landau_potential_and_histogram(
    int index_of_state_after_acceptance_or_rejection) {
  if (index_of_state_after_acceptance_or_rejection >= 0) {
    if (histogram[index_of_state_after_acceptance_or_rejection] >= 0) {
      histogram[index_of_state_after_acceptance_or_rejection] += 1;
      wang_landau_potential[index_of_state_after_acceptance_or_rejection] +=
          wang_landau_parameter;
    }
  }
}

/**
 *Determines whether we can reduce the Wang-Landau parameter
 */
bool WangLandauReactionEnsemble::can_refine_wang_landau_one_over_t() const {
  double minimum_required_value =
      0.80 * average_list_of_allowed_entries(
                 histogram); // This is an additional constraint to sample
                             // configuration space better. Use flatness
                             // criterion according to 1/t algorithm as long as
                             // you are not in 1/t regime.
  if (do_energy_reweighting)
    minimum_required_value = 20; // get faster in energy reweighting case

  return *(std::min_element(histogram.begin(), histogram.end())) >
             minimum_required_value ||
         m_system_is_in_1_over_t_regime;
}

/**
 *Reset the Wang-Landau histogram.
 */
void WangLandauReactionEnsemble::reset_histogram() {
  printf("Histogram is flat. Refining. Previous Wang-Landau modification "
         "parameter was %f.\n",
         wang_landau_parameter);
  fflush(stdout);

  for (int i = 0; i < wang_landau_potential.size(); i++) {
    if (histogram[i] >= 0) { // checks for validity of index i (think of energy
                             // collective variables, in a cubic memory layout
                             // there will be indices which are not allowed by
                             // the energy boundaries. These values will be
                             // initialized with a negative fill value)
      histogram[i] = 0;
    }
  }
}

/**
 *Refine the Wang-Landau parameter using the 1/t rule.
 */
void WangLandauReactionEnsemble::refine_wang_landau_parameter_one_over_t() {
  double monte_carlo_time = static_cast<double>(monte_carlo_trial_moves) /
                            static_cast<double>(used_bins);
  if (wang_landau_parameter / 2.0 <= 1.0 / monte_carlo_time ||
      m_system_is_in_1_over_t_regime) {
    wang_landau_parameter = 1.0 / monte_carlo_time;
    if (!m_system_is_in_1_over_t_regime) {
      m_system_is_in_1_over_t_regime = true;
      printf("Refining: Wang-Landau parameter is now 1/t.\n");
    }
  } else {
    reset_histogram();
    wang_landau_parameter = wang_landau_parameter / 2.0;
  }
}

/**
 *Determine whether the desired number of refinements was achieved.
 */
bool WangLandauReactionEnsemble::
    achieved_desired_number_of_refinements_one_over_t() const {
  if (wang_landau_parameter < final_wang_landau_parameter) {
    printf("Achieved desired number of refinements\n");
    return true;
  }
  return false;
}

void WangLandauReactionEnsemble::format_wang_landau_results(std::ostream &out) {
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    // unravel index
    if (std::abs(wang_landau_potential[flattened_index] - double_fill_value) >
        1) {
      // only output data if they are not equal to double_fill_value.
      // This ensures that for the energy observable not allowed energies
      // (energies in the interval [global_E_min, global_E_max]) in the
      // multidimensional Wang-Landau potential are printed out, since
      // the range [E_min(nbar), E_max(nbar)] for each nbar may be a
      // different one
      std::vector<std::size_t> unraveled_index(
          nr_subindices_of_collective_variable.size());
      Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                           nr_subindices_of_collective_variable.end(),
                           unraveled_index.begin(), flattened_index);
      // use unraveled index
      for (std::size_t i = 0; i < collective_variables.size(); i++) {
        auto const value = static_cast<double>(unraveled_index[i]) *
                               collective_variables[i]->delta_CV +
                           collective_variables[i]->CV_minimum;
        out << value << " ";
      }
      out << wang_landau_potential[flattened_index] << " \n";
    }
  }
}

/**
 *Writes the Wang-Landau potential to file.
 */
void WangLandauReactionEnsemble::write_wang_landau_results_to_file(
    const std::string &filename) {
  std::ofstream outfile;

  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }

  format_wang_landau_results(outfile);
  outfile.close();
}

/**
 *Update the minimum and maximum observed energies using the current state.
 *Needed for preliminary energy reweighting runs.
 */
int WangLandauReactionEnsemble::
    update_maximum_and_minimum_energies_at_current_state() {
  if (minimum_energies_at_flat_index.empty() ||
      maximum_energies_at_flat_index.empty()) {
    minimum_energies_at_flat_index.resize(wang_landau_potential.size(),
                                          double_fill_value);
    maximum_energies_at_flat_index.resize(wang_landau_potential.size(),
                                          double_fill_value);
  }

  const double E_pot_current = calculate_current_potential_energy_of_system();
  int index = get_flattened_index_wang_landau_of_current_state();

  // update stored energy values
  if (((E_pot_current < minimum_energies_at_flat_index[index]) ||
       std::abs(minimum_energies_at_flat_index[index] - double_fill_value) <
           std::numeric_limits<double>::epsilon())) {
    minimum_energies_at_flat_index[index] = E_pot_current;
  }
  if (((E_pot_current > maximum_energies_at_flat_index[index]) ||
       std::abs(maximum_energies_at_flat_index[index] - double_fill_value) <
           std::numeric_limits<double>::epsilon())) {
    maximum_energies_at_flat_index[index] = E_pot_current;
  }

  return 0;
}

/**
 *Write out an energy boundary file using the energy boundaries observed in a
 *preliminary energy reweighting run.
 */
void WangLandauReactionEnsemble::write_out_preliminary_energy_run_results(
    const std::string &filename) {

  std::ofstream outfile;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  outfile << "#nbar E_min E_max\n";
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    // unravel index
    std::vector<std::size_t> unraveled_index(
        nr_subindices_of_collective_variable.size());
    Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                         nr_subindices_of_collective_variable.end(),
                         unraveled_index.begin(), flattened_index);
    // use unraveled index
    for (std::size_t i = 0; i < collective_variables.size(); i++) {
      auto const value = static_cast<double>(unraveled_index[i]) *
                             collective_variables[i]->delta_CV +
                         collective_variables[i]->CV_minimum;
      outfile << value << " ";
    }
    outfile << minimum_energies_at_flat_index[flattened_index] << " "
            << maximum_energies_at_flat_index[flattened_index] << " \n";
  }
  outfile.close();
}

/**
 *Returns the flattened index of a given flattened index without the energy
 *collective variable.
 */
int WangLandauReactionEnsemble::
    get_flattened_index_wang_landau_without_energy_collective_variable(
        int flattened_index_with_EnergyCollectiveVariable,
        int CV_index_energy_observable) {
  std::vector<std::size_t> unraveled_index(
      nr_subindices_of_collective_variable.size());
  Utils::unravel_index(
      nr_subindices_of_collective_variable.begin(),
      nr_subindices_of_collective_variable.end(), unraveled_index.begin(),
      static_cast<std::size_t>(flattened_index_with_EnergyCollectiveVariable));
  // use unraveled index, but skip last collective variable
  // (the energy collective variable)
  auto const nr_collective_variables =
      static_cast<int>(collective_variables.size()) - 1;
  std::vector<double> current_state(nr_collective_variables);
  for (int i = 0; i < nr_collective_variables; i++) {
    current_state[i] = static_cast<double>(unraveled_index[i]) *
                           collective_variables[i]->delta_CV +
                       collective_variables[i]->CV_minimum;
  }

  // get collective_variables_minimum_values
  std::vector<double> collective_variables_minimum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_minimum_values[CV_i] =
        collective_variables[CV_i]->CV_minimum;
  }
  // get collective_variables_maximum_values
  std::vector<double> collective_variables_maximum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_maximum_values[CV_i] =
        collective_variables[CV_i]->CV_maximum;
  }
  // get delta_collective_variables_values
  std::vector<double> delta_collective_variables_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    delta_collective_variables_values[CV_i] =
        collective_variables[CV_i]->delta_CV;
  }
  int index = get_flattened_index_wang_landau(
      current_state, collective_variables_minimum_values,
      collective_variables_maximum_values, delta_collective_variables_values,
      nr_collective_variables);
  return index;
}

/**
 * Writes the Wang-Landau parameter, the histogram and the potential to a file.
 * You can restart a Wang-Landau simulation using this information.
 * Additionally you should store the positions of the particles.
 * Not storing them introduces small, small statistical errors.
 */
void WangLandauReactionEnsemble::write_wang_landau_checkpoint(
    const std::string &identifier) {
  std::string filename;
  std::ofstream outfile;

  // write current Wang-Landau parameters (wang_landau_parameter,
  // monte_carlo_trial_moves, flat_index_of_current_state)
  filename = std::string("checkpoint_wang_landau_parameters_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  outfile << wang_landau_parameter << " " << monte_carlo_trial_moves << " "
          << get_flattened_index_wang_landau_of_current_state() << "\n";
  outfile.close();

  // write histogram
  filename = std::string("checkpoint_wang_landau_histogram_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  for (int i = 0; i < wang_landau_potential.size(); i++) {
    outfile << histogram[i] << "\n";
  }
  outfile.close();
  // write Wang-Landau potential
  filename = std::string("checkpoint_wang_landau_potential_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  for (double i : wang_landau_potential) {
    outfile << i << "\n";
  }
  outfile.close();
}

/**
 *Loads the Wang-Landau checkpoint
 */
void WangLandauReactionEnsemble::load_wang_landau_checkpoint(
    const std::string &identifier) {
  std::string filename;
  std::ifstream infile;

  // restore Wang-Landau parameters
  filename = std::string("checkpoint_wang_landau_parameters_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    double wang_landau_parameter_entry;
    int wang_landau_monte_carlo_trial_moves_entry;
    int flat_index_of_state_at_checkpointing;
    int line = 0;
    while (infile >> wang_landau_parameter_entry >>
           wang_landau_monte_carlo_trial_moves_entry >>
           flat_index_of_state_at_checkpointing) {
      wang_landau_parameter = wang_landau_parameter_entry;
      monte_carlo_trial_moves = wang_landau_monte_carlo_trial_moves_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // restore histogram
  filename = std::string("checkpoint_wang_landau_histogram_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    int hist_entry;
    int line = 0;
    while (infile >> hist_entry) {
      histogram[line] = hist_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // restore Wang-Landau potential
  filename = std::string("checkpoint_wang_landau_potential_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    double wang_landau_potential_entry;
    int line = 0;
    while (infile >> wang_landau_potential_entry) {
      wang_landau_potential[line] = wang_landau_potential_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // possible task: restore state in which the system was when the checkpoint
  // was written. However as long as checkpointing and restoring the system form
  // the checkpoint is rare this should not matter statistically.
}

} // namespace ReactionMethods
