/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
    Implementation of nonbonded_interaction_data.hpp
 */
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "actor/DipolarDirectSum.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "dpd.hpp"
#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"
#include "electrostatics_magnetostatics/mdlc_correction.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "nonbonded_interaction_data.hpp"
#include "nonbonded_interactions/buckingham.hpp"
#include "nonbonded_interactions/cos2.hpp"
#include "nonbonded_interactions/gaussian.hpp"
#include "nonbonded_interactions/gb.hpp"
#include "nonbonded_interactions/hat.hpp"
#include "nonbonded_interactions/hertzian.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "nonbonded_interactions/ljcos.hpp"
#include "nonbonded_interactions/ljcos2.hpp"
#include "nonbonded_interactions/ljgen.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/steppot.hpp"
#ifdef DIPOLAR_BARNES_HUT
#include "actor/DipolarBarnesHut.hpp"
#endif
#include "layered.hpp"
#include "object-in-fluid/affinity.hpp"
#include "object-in-fluid/membrane_collision.hpp"
#include "pressure.hpp"
#include "rattle.hpp"
#include "serialization/IA_parameters.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <cstdlib>
#include <cstring>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"

/****************************************
 * variables
 *****************************************/
int max_seen_particle_type = 0;
std::vector<IA_parameters> ia_params;

double min_global_cut = 0.0;

double max_cut;
double max_cut_nonbonded;
double max_cut_bonded;
/** maximal cutoff of type-independent short range ia, mainly
    electrostatics and DPD*/
double max_cut_global;
/** Everything which is in the global cutoff except real space cutoffs
    of dipolar and Coulomb methods */
double max_cut_global_without_coulomb_and_dipolar;

/*****************************************
 * function prototypes
 *****************************************/

/** calculates and returns the maximal global nonbonded cutoff that is
    required.  Currently, this are just the cutoffs from the
    electrostatics method and some dpd cutoffs. */
static void recalc_global_maximal_nonbonded_cutoff();

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

static void recalc_global_maximal_nonbonded_and_long_range_cutoff() {
  /* user defined minimal global cut. This makes sure that data of
   pairs of particles with a distance smaller than this are always
   available on the same node (through ghosts). Required for example
   for the relative virtual sites algorithm. */
  max_cut_global = min_global_cut;

  // global cutoff without dipolar and Coulomb methods is needed
  // for more selective addition of particle pairs to Verlet lists
  max_cut_global_without_coulomb_and_dipolar = max_cut_global;

  // Electrostatics and magnetostatics

  /* Coulomb::Cutoff:
     Cutoff for the real space electrostatics.
     Note that the box length may have changed,
     but the method not yet reinitialized.
   */
#ifdef ELECTROSTATICS
  max_cut_global = std::max(max_cut_global, Coulomb::cutoff(box_l));
#endif

#ifdef DIPOLES
  max_cut_global = std::max(max_cut_global, Dipole::cutoff(box_l));
#endif
}

static void recalc_maximal_cutoff_nonbonded() {
  int i, j;

  recalc_global_maximal_nonbonded_and_long_range_cutoff();

  max_cut_nonbonded = max_cut_global;

  for (i = 0; i < max_seen_particle_type; i++)
    for (j = i; j < max_seen_particle_type; j++) {
      double max_cut_current = INACTIVE_CUTOFF;

      IA_parameters *data = get_ia_param(i, j);

#ifdef LENNARD_JONES
      if (max_cut_current < (data->LJ_cut + data->LJ_offset))
        max_cut_current = (data->LJ_cut + data->LJ_offset);
#endif

#ifdef WCA
      max_cut_current = std::max(max_cut_current, data->WCA_cut);
#endif

#ifdef DPD
      max_cut_current = std::max(max_cut_current,
                                 std::max(data->dpd_r_cut, data->dpd_tr_cut));
#endif

#ifdef LENNARD_JONES_GENERIC
      if (max_cut_current < (data->LJGEN_cut + data->LJGEN_offset))
        max_cut_current = (data->LJGEN_cut + data->LJGEN_offset);
#endif

#ifdef SMOOTH_STEP
      if (max_cut_current < data->SmSt_cut)
        max_cut_current = data->SmSt_cut;
#endif

#ifdef HERTZIAN
      if (max_cut_current < data->Hertzian_sig)
        max_cut_current = data->Hertzian_sig;
#endif

#ifdef GAUSSIAN
      if (max_cut_current < data->Gaussian_cut)
        max_cut_current = data->Gaussian_cut;
#endif

#ifdef BMHTF_NACL
      if (max_cut_current < data->BMHTF_cut)
        max_cut_current = data->BMHTF_cut;
#endif

#ifdef MORSE
      if (max_cut_current < data->MORSE_cut)
        max_cut_current = data->MORSE_cut;
#endif

#ifdef BUCKINGHAM
      if (max_cut_current < data->BUCK_cut)
        max_cut_current = data->BUCK_cut;
#endif

#ifdef SOFT_SPHERE
      if (max_cut_current < (data->soft_cut + data->soft_offset))
        max_cut_current = (data->soft_cut + data->soft_offset);
#endif

#ifdef AFFINITY
      if (max_cut_current < data->affinity_cut)
        max_cut_current = data->affinity_cut;
#endif

#ifdef MEMBRANE_COLLISION
      if (max_cut_current < data->membrane_cut)
        max_cut_current = data->membrane_cut;
#endif

#ifdef HAT
      if (max_cut_current < data->HAT_r)
        max_cut_current = data->HAT_r;
#endif

#ifdef LJCOS
      {
        double max_cut_tmp = data->LJCOS_cut + data->LJCOS_offset;
        if (max_cut_current < max_cut_tmp)
          max_cut_current = max_cut_tmp;
      }
#endif

#ifdef LJCOS2
      {
        double max_cut_tmp = data->LJCOS2_cut + data->LJCOS2_offset;
        if (max_cut_current < max_cut_tmp)
          max_cut_current = max_cut_tmp;
      }
#endif

#ifdef COS2
      {
        double max_cut_tmp = data->COS2_cut + data->COS2_offset;
        if (max_cut_current < max_cut_tmp)
          max_cut_current = max_cut_tmp;
      }
#endif

#ifdef GAY_BERNE
      if (max_cut_current < data->GB_cut)
        max_cut_current = data->GB_cut;
#endif

#ifdef TABULATED
      max_cut_current = std::max(max_cut_current, data->TAB.cutoff());
#endif

#ifdef SWIMMER_REACTIONS
      if (max_cut_current < data->REACTION_range)
        max_cut_current = data->REACTION_range;
#endif

#ifdef THOLE
      // If THOLE is active, use p3m cutoff
      if (data->THOLE_scaling_coeff != 0)
        max_cut_current = std::max(max_cut_current, Coulomb::cutoff(box_l));
#endif

      IA_parameters *data_sym = get_ia_param(j, i);

      /* no interaction ever touched it, at least no real
         short-ranged one (that writes to the nonbonded energy) */
      data_sym->particlesInteract = data->particlesInteract =
          (max_cut_current > 0.0);

      data_sym->max_cut = data->max_cut = max_cut_current;

      if (max_cut_current > max_cut_nonbonded)
        max_cut_nonbonded = max_cut_current;
    }
}

void recalc_maximal_cutoff() {
  recalc_maximal_cutoff_bonded();
  recalc_maximal_cutoff_nonbonded();

  /* make max_cut the maximal cutoff of both bonded and non-bonded
     interactions */
  if (max_cut_nonbonded > max_cut_bonded)
    max_cut = max_cut_nonbonded;
  else
    max_cut = max_cut_bonded;
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

  auto new_params = std::vector<IA_parameters>(nsize * nsize);

  /* if there is an old field, move entries */
  for (int i = 0; i < max_seen_particle_type; i++)
    for (int j = 0; j < max_seen_particle_type; j++) {
      new_params[i * nsize + j] =
          std::move(ia_params[i * max_seen_particle_type + j]);
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
