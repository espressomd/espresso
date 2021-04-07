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
#ifndef REACTION_METHODS_UTILS_HPP
#define REACTION_METHODS_UTILS_HPP

#include "reaction_methods/ReactionAlgorithm.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

namespace ReactionMethods {

/**
 * Calculates the whole product of factorial expressions which occur in the
 * reaction ensemble acceptance probability.
 *
 * See @cite smith94c.
 */
double
calculate_factorial_expression(SingleReaction const &current_reaction,
                               std::map<int, int> const &old_particle_numbers);

/**
 * Calculates the factorial expression which occurs in the reaction ensemble
 * acceptance probability
 */
double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i);

/**
 * Calculate the average of an array (used for the histogram of the
 * Wang-Landau algorithm). It excludes values which are initialized to be
 * negative. Those values indicate that the Wang-Landau algorithm should not
 * sample those values. The values still occur in the list because we can only
 * store "rectangular" value ranges.
 */
template <typename T>
double average_list_of_allowed_entries(const std::vector<T> &rng) {
  T result = 0;
  int counter_allowed_entries = 0;
  for (auto &val : rng) {
    if (val >= 0) { // checks for validity of index i (think of energy
                    // collective variables, in a cubic memory layout
                    // there will be indices which are not allowed by
                    // the energy boundaries. These values will be
                    // initialized with a negative fill value)
      result += val;
      counter_allowed_entries += 1;
    }
  }
  if (counter_allowed_entries) {
    return static_cast<double>(result) /
           static_cast<double>(counter_allowed_entries);
  }
  return 0.0;
}

/**
 * Finds the minimum non negative value in the provided range and returns
 * this value.
 */
inline double find_minimum_non_negative_value(std::vector<double> const &rng) {
  if (rng.empty())
    throw std::runtime_error("range is empty\n");
  // think of negative histogram values that indicate not
  // allowed energies in the case of an energy observable
  auto const it = std::min_element(rng.begin(), rng.end(),
                                   [](double const &a, double const &b) {
                                     if (a <= 0)
                                       return false;
                                     if (b <= 0)
                                       return true;
                                     return a < b;
                                   });
  if (it == rng.end() or *it < 0) {
    return rng.back();
  }
  return *it;
}

} // namespace ReactionMethods
#endif
