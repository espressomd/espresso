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
#ifndef _INTERACTION_DATA_H
#define _INTERACTION_DATA_H
/** \file
 *  Various procedures concerning interactions between particles.
 */

#include "TabulatedPotential.hpp"
#include "config.hpp"

#include <utils/index.hpp>

#include <MpiCallbacks.hpp>
#include <algorithm>
#include <cassert>
#include <communication.hpp>
#include <event.hpp>
#include <string>
#include <vector>

/** Cutoff for deactivated interactions. Must be negative, so that even
 *  particles on top of each other don't interact by chance.
 */
constexpr double INACTIVE_CUTOFF = -1.;

/** Lennard-Jones with shift */
struct LJ_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double shift = 0.0;
  double offset = 0.0;
  double min = 0.0;
};

/** WCA potential */
struct WCA_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
};

/** Generic Lennard-Jones with shift */
struct LJGen_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double shift = 0.0;
  double offset = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double lambda1 = 1.0;
  double softrad = 0.0;
};

/** smooth step potential */
struct SmoothStep_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double d = 0.0;
  int n = 0;
  double k0 = 0.0;
};

/** Hertzian potential */
struct Hertzian_Parameters {
  double eps = 0.0;
  double sig = INACTIVE_CUTOFF;
};

/** Gaussian potential */
struct Gaussian_Parameters {
  double eps = 0.0;
  double sig = 1.0;
  double cut = INACTIVE_CUTOFF;
};

/** BMHTF NaCl potential */
struct BMHTF_Parameters {
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double computed_shift = 0.0;
};

/** Morse potential */
struct Morse_Parameters {
  double eps = INACTIVE_CUTOFF;
  double alpha = INACTIVE_CUTOFF;
  double rmin = INACTIVE_CUTOFF;
  double cut = INACTIVE_CUTOFF;
  double rest = INACTIVE_CUTOFF;
};

/** Buckingham potential */
struct Buckingham_Parameters {
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double cut = INACTIVE_CUTOFF;
  double discont = 0.0;
  double shift = 0.0;
  double F1 = 0.0;
  double F2 = 0.0;
};

/** soft-sphere potential */
struct SoftSphere_Parameters {
  double a = 0.0;
  double n = 0.0;
  double cut = INACTIVE_CUTOFF;
  double offset = 0.0;
};

/** hat potential */
struct Hat_Parameters {
  double Fmax = 0.0;
  double r = INACTIVE_CUTOFF;
};

/** Lennard-Jones+Cos potential */
struct LJcos_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double offset = 0.0;
  double alfa = 0.0;
  double beta = 0.0;
  double rmin = 0.0;
};

/** Lennard-Jones with a different Cos potential */
struct LJcos2_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double offset = 0.0;
  double w = 0.0;
  double rchange = 0.0;
};

/** Gay-Berne potential */
struct GayBerne_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double k1 = 0.0;
  double k2 = 0.0;
  double mu = 0.0;
  double nu = 0.0;
  double chi1 = 0.0;
  double chi2 = 0.0;
};

/** Thole potential */
struct Thole_Parameters {
  double scaling_coeff;
  double q1q2;
};

/** DPD potential */
struct DPDParameters {
  double gamma = 0.;
  double k = 1.;
  double cutoff = -1.;
  int wf = 0;
  double pref = 0.0;
};

/** Data structure containing the interaction parameters for non-bonded
 *  interactions.
 *  Access via <tt>get_ia_param(i, j)</tt> with
 *  <tt>i</tt>, <tt>j</tt> \< \ref max_seen_particle_type
 */
struct IA_parameters {
  /** maximal cutoff for this pair of particle types. This contains
   *  contributions from the short-ranged interactions, plus any
   *  cutoffs from global interactions like electrostatics.
   */
  double max_cut = INACTIVE_CUTOFF;

#ifdef LENNARD_JONES
  LJ_Parameters lj;
#endif

#ifdef WCA
  WCA_Parameters wca;
#endif

#ifdef LENNARD_JONES_GENERIC
  LJGen_Parameters ljgen;
#endif

#ifdef SMOOTH_STEP
  SmoothStep_Parameters smooth_step;
#endif

#ifdef HERTZIAN
  Hertzian_Parameters hertzian;
#endif

#ifdef GAUSSIAN
  Gaussian_Parameters gaussian;
#endif

#ifdef BMHTF_NACL
  BMHTF_Parameters bmhtf;
#endif

#ifdef MORSE
  Morse_Parameters morse;
#endif

#ifdef BUCKINGHAM
  Buckingham_Parameters buckingham;
#endif

#ifdef SOFT_SPHERE
  SoftSphere_Parameters soft_sphere;
#endif

#ifdef HAT
  Hat_Parameters hat;
#endif

#ifdef LJCOS
  LJcos_Parameters ljcos;
#endif

#ifdef LJCOS2
  LJcos2_Parameters ljcos2;
#endif

#ifdef GAY_BERNE
  GayBerne_Parameters gay_berne;
#endif

#ifdef TABULATED
  TabulatedPotential tab;
#endif

#ifdef DPD
  /** \name DPD as interaction */
  /**@{*/
  DPDParameters dpd_radial;
  DPDParameters dpd_trans;
  /**@}*/
#endif

#ifdef THOLE
  Thole_Parameters thole;
#endif
};

extern std::vector<IA_parameters> ia_params;

/** Maximal particle type seen so far. */
extern int max_seen_particle_type;

/** Maximal interaction cutoff (real space/short range non-bonded
 *  interactions).
 */
double maximal_cutoff_nonbonded();

/** Minimal global interaction cutoff. Particles with a distance
 *  smaller than this are guaranteed to be available on the same node
 *  (through ghosts).
 */
extern double min_global_cut;

/**
 * @brief Get interaction parameters between particle types i and j
 *
 * This is symmetric, e.g. it holds that get_ia_param(i, j) and
 * get_ia_param(j, i) point to the same data.
 *
 * @param i First type, has to be smaller than @ref max_seen_particle_type.
 * @param j Second type, has to be smaller than @ref max_seen_particle_type.
 *
 * @return Pointer to interaction parameters for the type pair.
 */
inline IA_parameters *get_ia_param(int i, int j) {
  assert(i >= 0 && i < max_seen_particle_type);
  assert(j >= 0 && j < max_seen_particle_type);

  return &ia_params[Utils::upper_triangular(std::min(i, j), std::max(i, j),
                                            max_seen_particle_type)];
}

/** Get interaction parameters between particle types i and j.
 *  Slower than @ref get_ia_param, but can also be used on not
 *  yet present particle types
 */
IA_parameters *get_ia_param_safe(int i, int j);

/** @brief Get the state of all non-bonded interactions.
 */
std::string ia_params_get_state();

/** @brief Set the state of all non-bonded interactions.
 */
void ia_params_set_state(std::string const &);

bool is_new_particle_type(int type);
/** Make sure that ia_params is large enough to cover interactions
 *  for this particle type. The interactions are initialized with values
 *  such that no physical interaction occurs.
 */
void make_particle_type_exist(int type);

void make_particle_type_exist_local(int type);

/**
 * @brief Reset all interaction parameters to their defaults.
 */
void reset_ia_params();

/** Check if a non-bonded interaction is defined */
inline bool checkIfInteraction(IA_parameters const &data) {
  return data.max_cut != INACTIVE_CUTOFF;
}

/** @brief Set and broadcast the min_global_cut
 *  @param min_global_cut minimum global cutoff
 */
void mpi_set_min_global_cut(double min_global_cut);
#endif
