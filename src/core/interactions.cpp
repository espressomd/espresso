/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "electrostatics/coulomb.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/mpi/collectives/broadcast.hpp>

#include <algorithm>

static double recalc_long_range_cutoff() {
  auto max_cut_long_range = INACTIVE_CUTOFF;
#ifdef ELECTROSTATICS
  max_cut_long_range = std::max(max_cut_long_range, Coulomb::cutoff());
#endif

#ifdef DIPOLES
  max_cut_long_range = std::max(max_cut_long_range, Dipoles::cutoff());
#endif

  return max_cut_long_range;
}

double maximal_cutoff(bool single_node) {
  auto max_cut = get_min_global_cut();
  auto const max_cut_long_range = recalc_long_range_cutoff();
  auto const max_cut_bonded = maximal_cutoff_bonded();
  auto const max_cut_nonbonded = maximal_cutoff_nonbonded();

  max_cut = std::max(max_cut, max_cut_long_range);
  if (not single_node) {
    // If there is just one node, the bonded cutoff can be omitted
    // because bond partners are always on the local node.
    max_cut = std::max(max_cut, max_cut_bonded);
  }
  max_cut = std::max(max_cut, max_cut_nonbonded);
  max_cut = std::max(max_cut, collision_detection_cutoff());
  return max_cut;
}

bool long_range_interactions_sanity_checks() {
  try {
#ifdef ELECTROSTATICS
    Coulomb::sanity_checks();
#endif
#ifdef DIPOLES
    Dipoles::sanity_checks();
#endif
  } catch (std::runtime_error const &err) {
    runtimeErrorMsg() << err.what();
    return true;
  }
  return false;
}
