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
#ifndef REACTION_METHODS_WIDOM_INSERTION_HPP
#define REACTION_METHODS_WIDOM_INSERTION_HPP

#include "ReactionAlgorithm.hpp"

#include <utility>

namespace ReactionMethods {

/** Widom insertion method */
class WidomInsertion : public ReactionAlgorithm {
public:
  WidomInsertion(int seed) : ReactionAlgorithm(seed) {}
  double calculate_particle_insertion_potential_energy(
      SingleReaction &current_reaction);
};

} // namespace ReactionMethods
#endif
