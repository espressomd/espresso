/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/** \file
 *  Implementation of nonbonded_interaction_data.hpp
 */
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include "serialization/IA_parameters.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

/****************************************
 * variables
 *****************************************/
int max_seen_particle_type = 0;
std::vector<IA_parameters> ia_params;

double min_global_cut = INACTIVE_CUTOFF;

/*****************************************
 * function prototypes
 *****************************************/

/*****************************************
 * general low-level functions
 *****************************************/

IA_parameters *get_ia_param_safe(int i, int j) {
  make_particle_type_exist(std::max(i, j));
  return get_ia_param(i, j);
}

std::string ia_params_get_state() {
  std::stringstream out;
  boost::archive::binary_oarchive oa(out);
  oa << ia_params;
  oa << max_seen_particle_type;
  return out.str();
}

void ia_params_set_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);
  ia_params.clear();
  ia >> ia_params;
  ia >> max_seen_particle_type;
  mpi_bcast_max_seen_particle_type(max_seen_particle_type);
  mpi_bcast_all_ia_params();
}

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

static double recalc_maximal_cutoff(const IA_parameters &data) {
  auto max_cut_current = INACTIVE_CUTOFF;

#ifdef LENNARD_JONES
  max_cut_current = std::max(max_cut_current, (data.lj.cut + data.lj.offset));
#endif

#ifdef WCA
  max_cut_current = std::max(max_cut_current, data.wca.cut);
#endif

#ifdef DPD
  max_cut_current = std::max(
      max_cut_current, std::max(data.dpd_radial.cutoff, data.dpd_trans.cutoff));
#endif

#ifdef LENNARD_JONES_GENERIC
  max_cut_current =
      std::max(max_cut_current, (data.ljgen.cut + data.ljgen.offset));
#endif

#ifdef SMOOTH_STEP
  max_cut_current = std::max(max_cut_current, data.smooth_step.cut);
#endif

#ifdef HERTZIAN
  max_cut_current = std::max(max_cut_current, data.hertzian.sig);
#endif

#ifdef GAUSSIAN
  max_cut_current = std::max(max_cut_current, data.gaussian.cut);
#endif

#ifdef BMHTF_NACL
  max_cut_current = std::max(max_cut_current, data.bmhtf.cut);
#endif

#ifdef MORSE
  max_cut_current = std::max(max_cut_current, data.morse.cut);
#endif

#ifdef BUCKINGHAM
  max_cut_current = std::max(max_cut_current, data.buckingham.cut);
#endif

#ifdef SOFT_SPHERE
  max_cut_current = std::max(max_cut_current,
                             (data.soft_sphere.cut + data.soft_sphere.offset));
#endif

#ifdef HAT
  max_cut_current = std::max(max_cut_current, data.hat.r);
#endif

#ifdef LJCOS
  max_cut_current =
      std::max(max_cut_current, (data.ljcos.cut + data.ljcos.offset));
#endif

#ifdef LJCOS2
  max_cut_current =
      std::max(max_cut_current, (data.ljcos2.cut + data.ljcos2.offset));
#endif

#ifdef GAY_BERNE
  max_cut_current = std::max(max_cut_current, data.gay_berne.cut);
#endif

#ifdef TABULATED
  max_cut_current = std::max(max_cut_current, data.tab.cutoff());
#endif

#ifdef THOLE
  // If THOLE is active, use p3m cutoff
  if (data.thole.scaling_coeff != 0)
    max_cut_current =
        std::max(max_cut_current, Coulomb::cutoff(box_geo.length()));
#endif

  return max_cut_current;
}

double maximal_cutoff_nonbonded() {
  auto max_cut_nonbonded = INACTIVE_CUTOFF;

  for (auto &data : ia_params) {
    data.max_cut = recalc_maximal_cutoff(data);
    max_cut_nonbonded = std::max(max_cut_nonbonded, data.max_cut);
  }

  return max_cut_nonbonded;
}

double maximal_cutoff() {
  auto max_cut = min_global_cut;
  auto const max_cut_long_range = recalc_long_range_cutoff();
  auto const max_cut_bonded = maximal_cutoff_bonded();
  auto const max_cut_nonbonded = maximal_cutoff_nonbonded();

  max_cut = std::max(max_cut, max_cut_long_range);
  max_cut = std::max(max_cut, max_cut_bonded);
  max_cut = std::max(max_cut, max_cut_nonbonded);

  return max_cut;
}

/** This function increases the LOCAL ia_params field for non-bonded
   interactions
    to the given size. This function is not exported
    since it does not do this on all nodes. Use
    make_particle_type_exist for that.
*/
void realloc_ia_params(int nsize) {
  if (nsize <= max_seen_particle_type)
    return;

  auto new_params = std::vector<IA_parameters>(nsize * (nsize + 1) / 2);

  /* if there is an old field, move entries */
  for (int i = 0; i < max_seen_particle_type; i++)
    for (int j = i; j < max_seen_particle_type; j++) {
      new_params.at(Utils::upper_triangular(i, j, nsize)) =
          std::move(*get_ia_param(i, j));
    }

  max_seen_particle_type = nsize;
  std::swap(ia_params, new_params);
}

void reset_ia_params() {
  boost::fill(ia_params, IA_parameters{});
  mpi_bcast_all_ia_params();
}

bool is_new_particle_type(int type) {
  return (type + 1) > max_seen_particle_type;
}

void make_particle_type_exist(int type) {
  if (is_new_particle_type(type))
    mpi_bcast_max_seen_particle_type(type + 1);
}

void make_particle_type_exist_local(int type) {
  if (is_new_particle_type(type))
    realloc_ia_params(type + 1);
}

int interactions_sanity_checks() {
  /* set to zero if initialization was not successful. */
  int state = 1;

#ifdef ELECTROSTATICS
  Coulomb::sanity_checks(state);
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::nonbonded_sanity_check(state);
#endif /* ifdef  DIPOLES */

  return state;
}
