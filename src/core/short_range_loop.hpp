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
#ifndef CORE_SHORT_RANGE_HPP
#define CORE_SHORT_RANGE_HPP

#include "cells.hpp"
#include "grid.hpp"
#include "integrate.hpp"

#include <profiler/profiler.hpp>

#include <utility>

namespace detail {
/**
 * @brief Functor that returns true for
 *        any arguments.
 */
struct True {
  template <class... T> bool operator()(T...) const { return true; }
};
} // namespace detail

template <class PairKernel, class VerletCriterion = detail::True>
void short_range_loop(PairKernel pair_kernel,
                      const VerletCriterion &verlet_criterion) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (interaction_range() != INACTIVE_CUTOFF) {
    cell_structure.pair_loop(pair_kernel, verlet_criterion);
  }
}

template <class BondKernel, class PairKernel,
          class VerletCriterion = detail::True>
void short_range_loop(BondKernel bond_kernel, PairKernel pair_kernel,
                      const VerletCriterion &verlet_criterion) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  cell_structure.bond_loop(bond_kernel);

  short_range_loop(pair_kernel, verlet_criterion);
}
#endif
