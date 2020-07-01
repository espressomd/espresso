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
#ifndef _P3M_COMMON_H
#define _P3M_COMMON_H
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
#include "config.hpp"

#if defined(P3M) || defined(DP3M)

#include "LocalBox.hpp"

#include <utils/Vector.hpp>

/** Error Codes for p3m tuning (version 2) */
enum P3M_TUNE_ERROR {
  /** force evaluation failed */
  P3M_TUNE_FAIL = 1,
  /** could not find a valid realspace cutoff radius */
  P3M_TUNE_NOCUTOFF = 2,
  /** charge assignment order too large for mesh size */
  P3M_TUNE_CAO_TOO_LARGE = 4,
  /** conflict with ELC gap size */
  P3M_TUNE_ELCTEST = 8,
  P3M_TUNE_CUTOFF_TOO_LARGE = 16,
  /** could not achieve target accuracy */
  P3M_TUNE_ACCURACY_TOO_LARGE = 32
};

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/** increment size of charge assignment fields. */
#define CA_INCREMENT 32
/** precision limit for the r_cut zero */
#define P3M_RCUT_PREC 1e-3
/** granularity of the time measurement */
#define P3M_TIME_GRAN 2

/************************************************
 * data types
 ************************************************/

/** Structure for local mesh parameters. */
typedef struct {
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
} p3m_local_mesh;

/** Structure to hold P3M parameters and some dependent variables. */
typedef struct {
  /** tuning or production? */
  bool tuning = false;
  /** Ewald splitting parameter (0<alpha<1), rescaled to
   *  @p alpha_L = @p alpha * @p box_l. */
  double alpha_L = 0.0;
  /** cutoff radius for real space electrostatics (>0), rescaled to
   *  @p r_cut_iL = @p r_cut * @p box_l_i. */
  double r_cut_iL = 0.0;
  /** number of mesh points per coordinate direction (>0). */
  int mesh[3] = {};
  /** offset of the first mesh point (lower left corner) from the
   *  coordinate origin ([0,1[). */
  double mesh_off[3] = {P3M_MESHOFF, P3M_MESHOFF, P3M_MESHOFF};
  /** charge assignment order ([0,7]). */
  int cao = 0;
  /** accuracy of the actual parameter set. */
  double accuracy = 0.0;

  /** epsilon of the "surrounding dielectric". */
  double epsilon = P3M_EPSILON;
  /** cutoff for charge assignment. */
  double cao_cut[3] = {};
  /** mesh constant. */
  double a[3] = {};
  /** inverse mesh constant. */
  Utils::Vector3d ai = {};
  /** unscaled @ref P3MParameters::alpha_L "alpha_L" for use with fast
   *  inline functions only */
  double alpha = 0.0;
  /** unscaled @ref P3MParameters::r_cut_iL "r_cut_iL" for use with fast
   *  inline functions only */
  double r_cut = -1.;
  /** number of points unto which a single charge is interpolated, i.e.
   *  p3m.cao^3 */
  int cao3 = 0;
  /** additional points around the charge assignment mesh, for method like
   *  dielectric ELC creating virtual charges. */
  double additional_mesh[3] = {};

  template <typename Archive> void serialize(Archive &ar, long int) {
    ar &tuning &alpha_L &r_cut_iL &mesh;
    ar &mesh_off &cao &accuracy &epsilon &cao_cut;
    ar &a &ai &alpha &r_cut &cao3 &additional_mesh;
  }

} P3MParameters;

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  output array with dimension dim[3] at start position start[3].
 *
 *  \param in          Pointer to first element of input block data.
 *  \param out         Pointer to first element of output grid.
 *  \param start       Start position of block in output grid.
 *  \param size        Dimensions of the block
 *  \param dim         Dimensions of the output grid.
 */
void p3m_add_block(double const *in, double *out, int const start[3],
                   int const size[3], int const dim[3]);

/** One of the aliasing sums used by \ref p3m_k_space_error.
 *  Fortunately the one which is most important (because it converges
 *  most slowly, since it is not damped exponentially) can be
 *  calculated analytically. The result (which depends on the order of
 *  the spline interpolation) can be written as an even trigonometric
 *  polynomial. The results are tabulated here (the employed formula
 *  is eq. (7.66) in @cite hockney88a).
 */
double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao);

/** Calculate properties of the local FFT mesh for the
 *   charge assignment process.
 */
void p3m_calc_local_ca_mesh(p3m_local_mesh &local_mesh,
                            const P3MParameters &params,
                            const LocalBox<double> &local_geo, double skin);

/** Calculate the spatial position of the left down mesh
 *  point of the local mesh, to be stored in
 *  @ref p3m_local_mesh::ld_pos "ld_pos".
 *
 *  Function called by @ref p3m_calc_local_ca_mesh() once and by
 *  @ref p3m_scaleby_box_l() whenever the box size changes.
 */
void p3m_calc_lm_ld_pos(p3m_local_mesh &local_mesh,
                        const P3MParameters &params);

#endif /* P3M || DP3M */

#endif /* _P3M_COMMON_H */
