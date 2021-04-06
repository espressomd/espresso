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
#ifndef REACTION_METHODS_SINGLE_REACTION_HPP
#define REACTION_METHODS_SINGLE_REACTION_HPP

#include <utils/Accumulator.hpp>

#include <numeric>
#include <vector>

namespace ReactionMethods {

struct SingleReaction {
  SingleReaction() = default;
  SingleReaction(double gamma, std::vector<int> const &reactant_types,
                 std::vector<int> const &reactant_coefficients,
                 std::vector<int> const &product_types,
                 std::vector<int> const &product_coefficients) {
    this->reactant_types = reactant_types;
    this->reactant_coefficients = reactant_coefficients;
    this->product_types = product_types;
    this->product_coefficients = product_coefficients;
    this->gamma = gamma;
    nu_bar = std::accumulate(product_coefficients.begin(),
                             product_coefficients.end(), 0) -
             std::accumulate(reactant_coefficients.begin(),
                             reactant_coefficients.end(), 0);
  }

  // strict input to the algorithm
  std::vector<int> reactant_types;
  std::vector<int> reactant_coefficients;
  std::vector<int> product_types;
  std::vector<int> product_coefficients;
  double gamma = {};
  // calculated values that are stored for performance reasons
  int nu_bar = {}; ///< change in particle numbers for the reaction
  Utils::Accumulator accumulator_exponentials = Utils::Accumulator(1);
  int tried_moves = 0;
  int accepted_moves = 0;
  double get_acceptance_rate() const {
    return static_cast<double>(accepted_moves) /
           static_cast<double>(tried_moves);
  }
};

} // namespace ReactionMethods
#endif
