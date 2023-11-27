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

#pragma once

/** \file
 *  Various procedures concerning interactions between particles.
 */

#include "TabulatedPotential.hpp"
#include "config/config.hpp"

#include <utils/index.hpp>
#include <utils/math/int_pow.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
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
  LJ_Parameters() = default;
  LJ_Parameters(double epsilon, double sigma, double cutoff, double offset,
                double min, double shift);
  double get_auto_shift() const {
    auto auto_shift = 0.;
    if (cut != 0.) {
      auto_shift = Utils::int_pow<6>(sig / cut) - Utils::int_pow<12>(sig / cut);
    }
    return auto_shift;
  }
  double max_cutoff() const { return cut + offset; }
  double min_cutoff() const { return min + offset; }
};

/** WCA potential */
struct WCA_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  WCA_Parameters() = default;
  WCA_Parameters(double epsilon, double sigma);
  double max_cutoff() const { return cut; }
};

/** Generic Lennard-Jones with shift */
struct LJGen_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double shift = 0.0;
  double offset = 0.0;
  double lambda = 1.0;
  double softrad = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  LJGen_Parameters() = default;
  LJGen_Parameters(double epsilon, double sigma, double cutoff, double shift,
                   double offset,
#ifdef LJGEN_SOFTCORE
                   double lam, double delta,
#endif
                   double e1, double e2, double b1, double b2);
  double get_auto_shift() const {
    auto auto_shift = 0.;
    if (cut != 0.) {
      auto_shift = b2 * std::pow(sig / cut, a2) - b1 * std::pow(sig / cut, a1);
    }
    return auto_shift;
  }
  double max_cutoff() const { return cut + offset; }
};

/** smooth step potential */
struct SmoothStep_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double d = 0.0;
  int n = 0;
  double k0 = 0.0;
  SmoothStep_Parameters() = default;
  SmoothStep_Parameters(double eps, double sig, double cutoff, double d, int n,
                        double k0);
  double max_cutoff() const { return cut; }
};

/** Hertzian potential */
struct Hertzian_Parameters {
  double eps = 0.0;
  double sig = INACTIVE_CUTOFF;
  Hertzian_Parameters() = default;
  Hertzian_Parameters(double eps, double sig);
  double max_cutoff() const { return sig; }
};

/** Gaussian potential */
struct Gaussian_Parameters {
  double eps = 0.0;
  double sig = 1.0;
  double cut = INACTIVE_CUTOFF;
  Gaussian_Parameters() = default;
  Gaussian_Parameters(double eps, double sig, double cutoff);
  double max_cutoff() const { return cut; }
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
  BMHTF_Parameters() = default;
  BMHTF_Parameters(double A, double B, double C, double D, double sig,
                   double cut);
  double max_cutoff() const { return cut; }
};

/** Morse potential */
struct Morse_Parameters {
  double eps = 0.;
  double alpha = INACTIVE_CUTOFF;
  double rmin = INACTIVE_CUTOFF;
  double cut = INACTIVE_CUTOFF;
  double rest = INACTIVE_CUTOFF;
  Morse_Parameters() = default;
  Morse_Parameters(double eps, double alpha, double rmin, double cutoff);
  double max_cutoff() const { return cut; }
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
  Buckingham_Parameters() = default;
  Buckingham_Parameters(double a, double b, double c, double d, double cutoff,
                        double discont, double shift);
  double max_cutoff() const { return cut; }
};

/** soft-sphere potential */
struct SoftSphere_Parameters {
  double a = 0.0;
  double n = 0.0;
  double cut = INACTIVE_CUTOFF;
  double offset = 0.0;
  SoftSphere_Parameters() = default;
  SoftSphere_Parameters(double a, double n, double cutoff, double offset);
  double max_cutoff() const { return cut + offset; }
};

/** hat potential */
struct Hat_Parameters {
  double Fmax = 0.0;
  double r = INACTIVE_CUTOFF;
  Hat_Parameters() = default;
  Hat_Parameters(double F_max, double cutoff);
  double max_cutoff() const { return r; }
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
  LJcos_Parameters() = default;
  LJcos_Parameters(double epsilon, double sigma, double cutoff, double offset);
  double max_cutoff() const { return cut + offset; }
};

/** Lennard-Jones with a different Cos potential */
struct LJcos2_Parameters {
  double eps = 0.0;
  double sig = 0.0;
  double cut = INACTIVE_CUTOFF;
  double offset = 0.0;
  double w = 0.0;
  double rchange = 0.0;
  LJcos2_Parameters() = default;
  LJcos2_Parameters(double epsilon, double sigma, double offset, double width);
  double max_cutoff() const { return cut + offset; }
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
  GayBerne_Parameters() = default;
  GayBerne_Parameters(double eps, double sig, double cut, double k1, double k2,
                      double mu, double nu);
  double max_cutoff() const { return cut; }
};

/** Thole potential */
struct Thole_Parameters {
  double scaling_coeff = 0.; // inactive cutoff is 0
  double q1q2 = 0.;
  Thole_Parameters() = default;
  Thole_Parameters(double scaling_coeff, double q1q2)
      : scaling_coeff{scaling_coeff}, q1q2{q1q2} {}
};

/** DPD potential */
struct DPDParameters {
  double gamma = 0.;
  double k = 1.;
  double cutoff = INACTIVE_CUTOFF;
  int wf = 0;
  double pref = 0.0;
};

struct DPD_Parameters {
  DPDParameters radial;
  DPDParameters trans;
  DPD_Parameters() = default;
  DPD_Parameters(double gamma, double k, double r_c, int wf, double tgamma,
                 double tr_c, int twf) {
    radial = DPDParameters{gamma, k, r_c, wf, -1.};
    trans = DPDParameters{tgamma, k, tr_c, twf, -1.};
  }
  double max_cutoff() const { return std::max(radial.cutoff, trans.cutoff); }
};

/** @brief Parameters for non-bonded interactions. */
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
  DPD_Parameters dpd;
#endif

#ifdef THOLE
  Thole_Parameters thole;
#endif
};

class InteractionsNonBonded {
  /** @brief List of pairwise interactions. */
  std::vector<std::shared_ptr<IA_parameters>> m_nonbonded_ia_params{};
  /** @brief Maximal particle type seen so far. */
  int max_seen_particle_type = -1;

  void realloc_ia_params(int type) {
    assert(type >= 0);
    auto const old_size = m_nonbonded_ia_params.size();
    m_nonbonded_ia_params.resize(Utils::lower_triangular(type, type) + 1);
    auto const new_size = m_nonbonded_ia_params.size();
    if (new_size > old_size) {
      for (auto &data : m_nonbonded_ia_params) {
        if (data == nullptr) {
          data = std::make_shared<IA_parameters>();
        }
      }
    }
  }

public:
  InteractionsNonBonded() {
    /* make sure interaction 0<->0 always exists */
    make_particle_type_exist(0);
  }

  /**
   * @brief Make sure the interaction parameter list is large enough to cover
   * interactions for this particle type.
   * New interactions are initialized with values such that no physical
   * interaction occurs.
   */
  void make_particle_type_exist(int type) {
    assert(type >= 0);
    if (type > max_seen_particle_type) {
      realloc_ia_params(type);
      max_seen_particle_type = type;
    }
  }

  auto get_ia_param_key(int i, int j) const {
    assert(i >= 0 and i <= max_seen_particle_type);
    assert(j >= 0 and j <= max_seen_particle_type);
    auto const key = static_cast<std::size_t>(
        Utils::lower_triangular(std::max(i, j), std::min(i, j)));
    assert(key < m_nonbonded_ia_params.size());
    return key;
  }

  /**
   * @brief Get interaction parameters between particle types i and j
   *
   * This is symmetric, e.g. it holds that `get_ia_param(i, j)` and
   * `get_ia_param(j, i)` point to the same data.
   *
   * @param i First type, must exist
   * @param j Second type, must exist
   *
   * @return Reference to interaction parameters for the type pair.
   */
  auto &get_ia_param(int i, int j) {
    return *m_nonbonded_ia_params[get_ia_param_key(i, j)];
  }

  auto const &get_ia_param(int i, int j) const {
    return *m_nonbonded_ia_params[get_ia_param_key(i, j)];
  }

  auto get_ia_param_ref_counted(int i, int j) const {
    return m_nonbonded_ia_params[get_ia_param_key(i, j)];
  }

  void set_ia_param(int i, int j, std::shared_ptr<IA_parameters> const &ia) {
    m_nonbonded_ia_params[get_ia_param_key(i, j)] = ia;
  }

  auto get_max_seen_particle_type() const { return max_seen_particle_type; }

  /** @brief Recalculate cutoff of each interaction struct. */
  void recalc_maximal_cutoffs();

  /** @brief Get maximal cutoff. */
  double maximal_cutoff() const;

  double maximal_cutoff(std::vector<int> types) const {
    auto max_cut_nonbonded = INACTIVE_CUTOFF;
    auto const n_types = static_cast<int>(types.size());
    for (int i = 0; i < n_types; i++) {
      for (int j = i; j < n_types; j++) {
        if (types[i] > max_seen_particle_type or
            types[j] > max_seen_particle_type) {
          continue;
        }
        auto const data = get_ia_param(types[i], types[j]);
        max_cut_nonbonded = std::max(max_cut_nonbonded, data.max_cut);
      }
    }
    return max_cut_nonbonded;
  }

  double maximal_cutoff_all_except(std::vector<int> types) const {
    auto max_cut_nonbonded = INACTIVE_CUTOFF;
    for (int i = 0; i < max_seen_particle_type; i++) {
      if (std::find(types.begin(), types.end(), i) != types.end()) {
        // i in types
        continue;
      }
      for (int j = i; j < max_seen_particle_type; j++) {
        if (std::find(types.begin(), types.end(), j) != types.end()) {
          // j in types
          continue;
        }
        auto const data = get_ia_param(i, j);
        max_cut_nonbonded = std::max(max_cut_nonbonded, data.max_cut);
      }
    }
    return max_cut_nonbonded;
  }

  bool only_central_forces(std::vector<int> types) const {
    // gay-berne is the only non-central force
    bool ret = true;
#ifdef GAY_BERNE
    auto const n_types = static_cast<int>(types.size());
    for (int i = 0; i < n_types; i++) {
      for (int j = i; j < n_types; j++) {
        auto const data = get_ia_param(i, j);
        if (data.gay_berne.cut != INACTIVE_CUTOFF) {
          ret = false;
          break;
        }
      }
    }
#endif
    return ret;
  }
};
