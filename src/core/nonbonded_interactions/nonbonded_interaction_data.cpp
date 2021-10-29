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

#include "communication.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "serialization/IA_parameters.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <utils/index.hpp>

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

/****************************************
 * variables
 *****************************************/
int max_seen_particle_type = 0;
std::vector<IA_parameters> ia_params;

double min_global_cut = INACTIVE_CUTOFF;

/*****************************************
 * general low-level functions
 *****************************************/

static void mpi_realloc_ia_params_local(int new_size) {
  if (new_size <= max_seen_particle_type)
    return;

  auto new_params = std::vector<IA_parameters>(new_size * (new_size + 1) / 2);

  /* if there is an old field, move entries */
  for (int i = 0; i < max_seen_particle_type; i++)
    for (int j = i; j < max_seen_particle_type; j++) {
      new_params.at(Utils::upper_triangular(i, j, new_size)) =
          std::move(*get_ia_param(i, j));
    }

  max_seen_particle_type = new_size;
  std::swap(ia_params, new_params);
}

REGISTER_CALLBACK(mpi_realloc_ia_params_local)

/** Increase the size of the @ref ia_params vector. */
inline void mpi_realloc_ia_params(int new_size) {
  mpi_call_all(mpi_realloc_ia_params_local, new_size);
}

static void mpi_bcast_all_ia_params_local() {
  boost::mpi::broadcast(comm_cart, ia_params, 0);
}

REGISTER_CALLBACK(mpi_bcast_all_ia_params_local)

/** Broadcast @ref ia_params to all nodes. */
inline void mpi_bcast_all_ia_params() {
  mpi_call_all(mpi_bcast_all_ia_params_local);
}

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
  mpi_realloc_ia_params(max_seen_particle_type);
  mpi_bcast_all_ia_params();
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

void reset_ia_params() {
  boost::fill(ia_params, IA_parameters{});
  mpi_bcast_all_ia_params();
}

bool is_new_particle_type(int type) {
  return (type + 1) > max_seen_particle_type;
}

void make_particle_type_exist(int type) {
  if (is_new_particle_type(type))
    mpi_realloc_ia_params(type + 1);
}

void make_particle_type_exist_local(int type) {
  if (is_new_particle_type(type))
    mpi_realloc_ia_params_local(type + 1);
}

void mpi_set_min_global_cut_local(double min_global_cut) {
  ::min_global_cut = min_global_cut;
  on_skin_change();
}

REGISTER_CALLBACK(mpi_set_min_global_cut_local)

void mpi_set_min_global_cut(double min_global_cut) {
  mpi_call_all(mpi_set_min_global_cut_local, min_global_cut);
}
