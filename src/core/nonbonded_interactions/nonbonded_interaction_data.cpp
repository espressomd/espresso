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
/** \file
 *  Implementation of nonbonded_interaction_data.hpp
 */
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "communication.hpp"
#include "electrostatics/coulomb.hpp"
#include "event.hpp"
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
std::vector<IA_parameters> old_nonbonded_ia_params;
std::vector<std::shared_ptr<IA_parameters>> nonbonded_ia_params;

/** Minimal global interaction cutoff. Particles with a distance
 *  smaller than this are guaranteed to be available on the same node
 *  (through ghosts).
 */
static double min_global_cut = INACTIVE_CUTOFF;

/*****************************************
 * general low-level functions
 *****************************************/

void mpi_realloc_ia_params_local(int new_size) {
  auto const old_size = ::max_seen_particle_type;
  if (new_size <= old_size)
    return;

  auto const n_pairs = new_size * (new_size + 1) / 2;
  auto new_params_ = std::vector<IA_parameters>(n_pairs);
  auto new_params = std::vector<std::shared_ptr<IA_parameters>>(n_pairs);

  /* if there is an old field, move entries */
  for (int i = 0; i < old_size; i++) {
    for (int j = i; j < old_size; j++) {
      auto const new_key = Utils::upper_triangular(i, j, new_size);
      auto const old_key = Utils::upper_triangular(i, j, old_size);
      new_params_.at(new_key) = std::move(old_nonbonded_ia_params[old_key]);
      new_params[new_key] = std::move(nonbonded_ia_params[old_key]);
    }
  }
  for (auto &ia_params : new_params) {
    if (ia_params == nullptr) {
      ia_params = std::make_shared<IA_parameters>();
    }
  }

  ::max_seen_particle_type = new_size;
  std::swap(::old_nonbonded_ia_params, new_params_);
  std::swap(::nonbonded_ia_params, new_params);
}

REGISTER_CALLBACK(mpi_realloc_ia_params_local)

/** Increase the size of the @ref nonbonded_ia_params vector. */
static void mpi_realloc_ia_params(int new_size) {
  mpi_call_all(mpi_realloc_ia_params_local, new_size);
}

static void mpi_bcast_all_ia_params_local() {
  boost::mpi::broadcast(comm_cart, old_nonbonded_ia_params, 0);
}

REGISTER_CALLBACK(mpi_bcast_all_ia_params_local)

/** Broadcast @ref old_nonbonded_ia_params to all nodes. */
static void mpi_bcast_all_ia_params() {
  mpi_call_all(mpi_bcast_all_ia_params_local);
}

IA_parameters *get_ia_param_safe(int i, int j) {
  make_particle_type_exist(std::max(i, j));
  return &get_ia_param(i, j);
}

std::string ia_params_get_state() {
  std::stringstream out;
  boost::archive::binary_oarchive oa(out);
  oa << old_nonbonded_ia_params;
  oa << max_seen_particle_type;
  return out.str();
}

void ia_params_set_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);
  old_nonbonded_ia_params.clear();
  ia >> old_nonbonded_ia_params;
  ia >> max_seen_particle_type;
  mpi_realloc_ia_params(max_seen_particle_type);
  mpi_bcast_all_ia_params();
}

static double recalc_maximal_cutoff(const IA_parameters &data) {
  auto max_cut_current = INACTIVE_CUTOFF;

#ifdef LENNARD_JONES
  max_cut_current = std::max(max_cut_current, data.lj.max_cutoff());
#endif

#ifdef WCA
  max_cut_current = std::max(max_cut_current, data.wca.max_cutoff());
#endif

#ifdef DPD
  max_cut_current = std::max(max_cut_current, data.dpd.max_cutoff());
#endif

#ifdef LENNARD_JONES_GENERIC
  max_cut_current = std::max(max_cut_current, data.ljgen.max_cutoff());
#endif

#ifdef SMOOTH_STEP
  max_cut_current = std::max(max_cut_current, data.smooth_step.cut);
#endif

#ifdef HERTZIAN
  max_cut_current = std::max(max_cut_current, data.hertzian.max_cutoff());
#endif

#ifdef GAUSSIAN
  max_cut_current = std::max(max_cut_current, data.gaussian.max_cutoff());
#endif

#ifdef BMHTF_NACL
  max_cut_current = std::max(max_cut_current, data.bmhtf.max_cutoff());
#endif

#ifdef MORSE
  max_cut_current = std::max(max_cut_current, data.morse.max_cutoff());
#endif

#ifdef BUCKINGHAM
  max_cut_current = std::max(max_cut_current, data.buckingham.max_cutoff());
#endif

#ifdef SOFT_SPHERE
  max_cut_current = std::max(max_cut_current, data.soft_sphere.max_cutoff());
#endif

#ifdef HAT
  max_cut_current = std::max(max_cut_current, data.hat.max_cutoff());
#endif

#ifdef LJCOS
  max_cut_current = std::max(max_cut_current, data.ljcos.max_cutoff());
#endif

#ifdef LJCOS2
  max_cut_current = std::max(max_cut_current, data.ljcos2.max_cutoff());
#endif

#ifdef GAY_BERNE
  max_cut_current = std::max(max_cut_current, data.gay_berne.max_cutoff());
#endif

#ifdef TABULATED
  max_cut_current = std::max(max_cut_current, data.tab.cutoff());
#endif

#ifdef THOLE
  // If THOLE is active, use p3m cutoff
  if (data.thole.scaling_coeff != 0.)
    max_cut_current = std::max(max_cut_current, Coulomb::cutoff());
#endif

  return max_cut_current;
}

double maximal_cutoff_nonbonded() {
  auto max_cut_nonbonded = INACTIVE_CUTOFF;

  for (auto &data : old_nonbonded_ia_params) {
    data.max_cut = recalc_maximal_cutoff(data);
    max_cut_nonbonded = std::max(max_cut_nonbonded, data.max_cut);
  }

  for (auto &data : nonbonded_ia_params) {
    data->max_cut = recalc_maximal_cutoff(*data);
    max_cut_nonbonded = std::max(max_cut_nonbonded, data->max_cut);
  }

  return max_cut_nonbonded;
}

void reset_ia_params() {
  boost::fill(old_nonbonded_ia_params, IA_parameters{});
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

void set_min_global_cut(double min_global_cut) {
  ::min_global_cut = min_global_cut;
  on_skin_change();
}

double get_min_global_cut() { return ::min_global_cut; }
