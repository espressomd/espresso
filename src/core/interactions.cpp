/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#include "communication.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "collision.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "event.hpp"
#include "grid.hpp"

#include "serialization/IA_parameters.hpp"

#include <boost/mpi.hpp>

#include <algorithm>

static double recalc_long_range_cutoff() {
  auto max_cut_long_range = INACTIVE_CUTOFF;
#ifdef ELECTROSTATICS
  max_cut_long_range =
      std::max(max_cut_long_range, Coulomb::cutoff(box_geo.length()));
#endif

#ifdef DIPOLES
  max_cut_long_range =
      std::max(max_cut_long_range, Dipole::cutoff(box_geo.length()));
#endif

  return max_cut_long_range;
}

double maximal_cutoff() {
  auto max_cut = min_global_cut;
  auto const max_cut_long_range = recalc_long_range_cutoff();
  auto const max_cut_bonded = maximal_cutoff_bonded();
  auto const max_cut_nonbonded = maximal_cutoff_nonbonded();

  max_cut = std::max(max_cut, max_cut_long_range);
  max_cut = std::max(max_cut, max_cut_bonded);
  max_cut = std::max(max_cut, max_cut_nonbonded);
  max_cut = std::max(max_cut, collision_detection_cutoff());

  return max_cut;
}

bool long_range_interactions_sanity_checks() {
  /* set to zero if initialization was not successful. */
  bool failed = false;

#ifdef ELECTROSTATICS
  failed |= Coulomb::sanity_checks();
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  failed |= Dipole::sanity_checks();
#endif /* ifdef DIPOLES */

  return failed;
}

void mpi_bcast_ia_params_local(int i, int j) {
  boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);
  on_short_range_ia_change();
}

REGISTER_CALLBACK(mpi_bcast_ia_params_local)

void mpi_bcast_ia_params(int i, int j) {
  mpi_call_all(mpi_bcast_ia_params_local, i, j);
}
