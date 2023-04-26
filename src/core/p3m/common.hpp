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
 *  Common functions for dipolar and charge P3M.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  @cite hockney88a and @cite deserno98a @cite deserno98b. The file p3m
 *  contains only the Particle-Mesh part.
 *
 *  Further reading: @cite ewald21a, @cite hockney88a, @cite deserno98a,
 *  @cite deserno98b, @cite deserno00e, @cite deserno00b, @cite cerda08d
 *
 */

#ifndef ESPRESSO_SRC_CORE_P3M_COMMON_HPP
#define ESPRESSO_SRC_CORE_P3M_COMMON_HPP

#include "config/config.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <vector>

/** This value indicates metallic boundary conditions. */
auto constexpr P3M_EPSILON_METALLIC = 0.0;

#if defined(P3M) || defined(DP3M)

#include "LocalBox.hpp"

#include <array>
#include <stdexcept>
#include <vector>

namespace detail {
/** @brief Index helpers for direct and reciprocal space.
 *  After the FFT the data is in order YZX, which
 *  means that Y is the slowest changing index.
 */
namespace FFT_indexing {
enum FFT_REAL_VECTOR : int { RX = 0, RY = 1, RZ = 2 };
enum FFT_WAVE_VECTOR : int { KY = 0, KZ = 1, KX = 2 };
} // namespace FFT_indexing
} // namespace detail

/** Structure to hold P3M parameters and some dependent variables. */
struct P3MParameters {
  /** tuning or production? */
  bool tuning;
  /** Ewald splitting parameter (0<alpha<1), rescaled to
   *  @p alpha_L = @p alpha * @p box_l. */
  double alpha_L;
  /** cutoff radius for real space electrostatics (>0), rescaled to
   *  @p r_cut_iL = @p r_cut * @p box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  Utils::Vector3i mesh = {};
  /** offset of the first mesh point (lower left corner) from the
   *  coordinate origin ([0,1[). */
  Utils::Vector3d mesh_off;
  /** charge assignment order ([0,7]). */
  int cao;
  /** accuracy of the actual parameter set. */
  double accuracy;

  /** epsilon of the "surrounding dielectric". */
  double epsilon;
  /** cutoff for charge assignment. */
  Utils::Vector3d cao_cut;
  /** mesh constant. */
  Utils::Vector3d a;
  /** inverse mesh constant. */
  Utils::Vector3d ai;
  /** unscaled @ref P3MParameters::alpha_L "alpha_L" for use with fast
   *  inline functions only */
  double alpha;
  /** unscaled @ref P3MParameters::r_cut_iL "r_cut_iL" for use with fast
   *  inline functions only */
  double r_cut;
  /** number of points unto which a single charge is interpolated, i.e.
   *  @ref P3MParameters::cao "cao" cubed */
  int cao3;

  P3MParameters(bool tuning, double epsilon, double r_cut,
                Utils::Vector3i const &mesh, Utils::Vector3d const &mesh_off,
                int cao, double alpha, double accuracy)
      : tuning{tuning}, alpha_L{0.}, r_cut_iL{0.}, mesh{mesh},
        mesh_off{mesh_off}, cao{cao}, accuracy{accuracy}, epsilon{epsilon},
        cao_cut{}, a{}, ai{}, alpha{alpha}, r_cut{r_cut}, cao3{-1} {

    auto constexpr value_to_tune = -1.;

    if (epsilon < 0.) {
      throw std::domain_error("Parameter 'epsilon' must be >= 0");
    }

    if (accuracy <= 0.) {
      throw std::domain_error("Parameter 'accuracy' must be > 0");
    }

    if (r_cut <= 0.) {
      if (tuning and r_cut == value_to_tune) {
        this->r_cut = 0.;
      } else {
        throw std::domain_error("Parameter 'r_cut' must be > 0");
      }
    }

    if (alpha <= 0.) {
      if (tuning and alpha == value_to_tune) {
        this->alpha = 0.;
      } else {
        throw std::domain_error("Parameter 'alpha' must be > 0");
      }
    }

    if (not(mesh >= Utils::Vector3i::broadcast(1) or
            ((mesh[0] >= 1) and (mesh == Utils::Vector3i{{mesh[0], -1, -1}})) or
            (tuning and mesh == Utils::Vector3i::broadcast(-1)))) {
      throw std::domain_error("Parameter 'mesh' must be > 0");
    }

    if (not(mesh_off >= Utils::Vector3d::broadcast(0.) and
            mesh_off <= Utils::Vector3d::broadcast(1.))) {
      if (mesh_off == Utils::Vector3d::broadcast(-1.)) {
        this->mesh_off = Utils::Vector3d::broadcast(P3M_MESHOFF);
      } else {
        throw std::domain_error("Parameter 'mesh_off' must be >= 0 and <= 1");
      }
    }

    if ((cao < 1 or cao > 7) and (not tuning or cao != -1)) {
      throw std::domain_error("Parameter 'cao' must be >= 1 and <= 7");
    }

    if (not tuning and (Utils::Vector3i::broadcast(cao) > mesh)) {
      throw std::domain_error("Parameter 'cao' cannot be larger than 'mesh'");
    }
  }

  /**
   * @brief Recalculate quantities derived from the mesh and box length:
   * @ref P3MParameters::a "a",
   * @ref P3MParameters::ai "ai" and
   * @ref P3MParameters::cao_cut "cao_cut".
   */
  void recalc_a_ai_cao_cut(Utils::Vector3d const &box_l) {
    ai = Utils::hadamard_division(mesh, box_l);
    a = Utils::hadamard_division(Utils::Vector3d::broadcast(1.), ai);
    cao_cut = (static_cast<double>(cao) / 2.) * a;
  }
};

/** Structure for local mesh parameters. */
struct P3MLocalMesh {
  /* local mesh characterization. */
  /** dimension (size) of local mesh. */
  Utils::Vector3i dim;
  /** number of local mesh points. */
  int size;
  /** index of lower left corner of the
      local mesh in the global mesh. */
  int ld_ind[3];
  /** position of the first local mesh point. */
  double ld_pos[3];
  /** dimension of mesh inside node domain. */
  int inner[3];
  /** inner left down grid point */
  int in_ld[3];
  /** inner up right grid point + (1,1,1) */
  int in_ur[3];
  /** number of margin mesh points. */
  int margin[6];
  /** number of margin mesh points from neighbour nodes */
  int r_margin[6];
  /** offset between mesh lines of the last dimension */
  int q_2_off;
  /** offset between mesh lines of the two last dimensions */
  int q_21_off;

  /**
   * @brief Recalculate quantities derived from the mesh and box length:
   * @ref P3MLocalMesh::ld_pos "ld_pos" (position of the left down mesh).
   */
  void recalc_ld_pos(P3MParameters const &params) {
    // spatial position of left down mesh point
    for (unsigned int i = 0; i < 3; i++) {
      ld_pos[i] = (ld_ind[i] + params.mesh_off[i]) * params.a[i];
    }
  }

  /**
   * @brief Calculate properties of the local FFT mesh
   * for the charge assignment process.
   */
  void calc_local_ca_mesh(P3MParameters const &params,
                          LocalBox<double> const &local_geo, double skin,
                          double space_layer);
};

/** One of the aliasing sums used to compute k-space errors.
 *  Fortunately the one which is most important (because it converges
 *  most slowly, since it is not damped exponentially) can be
 *  calculated analytically. The result (which depends on the order of
 *  the spline interpolation) can be written as an even trigonometric
 *  polynomial. The results are tabulated here (the employed formula
 *  is eq. (7.66) in @cite hockney88a).
 */
double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao);

#endif /* P3M || DP3M */

namespace detail {
/** Calculate indices that shift @ref P3MParameters::mesh "mesh" by `mesh/2`.
 *  For each mesh size @f$ n @f$ in @c mesh_size, create a sequence of integer
 *  values @f$ \left( 0, \ldots, \lfloor n/2 \rfloor, -\lfloor n/2 \rfloor,
 *  \ldots, -1\right) @f$ if @c zero_out_midpoint is false, otherwise
 *  @f$ \left( 0, \ldots, \lfloor n/2 - 1 \rfloor, 0, -\lfloor n/2 \rfloor,
 *  \ldots, -1\right) @f$.
 */
std::array<std::vector<int>, 3> inline calc_meshift(
    Utils::Vector3i const &mesh_size, bool zero_out_midpoint = false) {
  std::array<std::vector<int>, 3> ret{};

  for (unsigned int i = 0; i < 3; i++) {
    ret[i] = std::vector<int>(mesh_size[i]);

    for (int j = 1; j <= mesh_size[i] / 2; j++) {
      ret[i][j] = j;
      ret[i][mesh_size[i] - j] = -j;
    }
    if (zero_out_midpoint)
      ret[i][mesh_size[i] / 2] = 0;
  }

  return ret;
}
} // namespace detail

#endif
