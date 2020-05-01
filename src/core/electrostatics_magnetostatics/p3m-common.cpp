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
 *  P3M main file.
 */
#include "p3m-common.hpp"

#if defined(P3M) || defined(DP3M)
#include "errorhandling.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <cstdio>

/* For debug messages */
extern int this_node;

void p3m_add_block(double const *in, double *out, int const start[3],
                   int const size[3], int const dim[3]) {
  /* fast,mid and slow changing indices */
  int f, m, s;
  /* linear index of in grid, linear index of out grid */
  int li_in = 0, li_out = 0;
  /* offsets for indices in output grid */
  int m_out_offset, s_out_offset;

  li_out = start[2] + (dim[2] * (start[1] + (dim[1] * start[0])));
  m_out_offset = dim[2] - size[2];
  s_out_offset = (dim[2] * (dim[1] - size[1]));

  for (s = 0; s < size[0]; s++) {
    for (m = 0; m < size[1]; m++) {
      for (f = 0; f < size[2]; f++) {
        out[li_out++] += in[li_in++];
      }
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao) {
  double c, res = 0.0;
  c = Utils::sqr(cos(Utils::pi() * mesh_i * (double)n));

  switch (cao) {
  case 1: {
    res = 1;
    break;
  }
  case 2: {
    res = (1.0 + c * 2.0) / 3.0;
    break;
  }
  case 3: {
    res = (2.0 + c * (11.0 + c * 2.0)) / 15.0;
    break;
  }
  case 4: {
    res = (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
    break;
  }
  case 5: {
    res = (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0)))) / 2835.0;
    break;
  }
  case 6: {
    res = (1382.0 +
           c * (35396.0 +
                c * (83021.0 + c * (34096.0 + c * (2026.0 + c * 4.0))))) /
          155925.0;
    break;
  }
  case 7: {
    res =
        (21844.0 +
         c * (776661.0 + c * (2801040.0 +
                              c * (2123860.0 +
                                   c * (349500.0 + c * (8166.0 + c * 4.0)))))) /
        6081075.0;
    break;
  }
  default: {
    fprintf(stderr,
            "%d: INTERNAL_ERROR: The value %d for the interpolation order "
            "should not occur!\n",
            this_node, cao);
    errexit();
  }
  }

  return res;
}

void p3m_calc_local_ca_mesh(p3m_local_mesh &local_mesh,
                            const P3MParameters &params,
                            const LocalBox<double> &local_geo, double skin) {
  int i;
  int ind[3];
  /* total skin size */
  double full_skin[3];

  for (i = 0; i < 3; i++)
    full_skin[i] = params.cao_cut[i] + skin + params.additional_mesh[i];

  /* inner left down grid point (global index) */
  for (i = 0; i < 3; i++)
    local_mesh.in_ld[i] =
        (int)ceil(local_geo.my_left()[i] * params.ai[i] - params.mesh_off[i]);
  /* inner up right grid point (global index) */
  for (i = 0; i < 3; i++)
    local_mesh.in_ur[i] =
        (int)floor(local_geo.my_right()[i] * params.ai[i] - params.mesh_off[i]);

  /* correct roundoff errors at boundary */
  for (i = 0; i < 3; i++) {
    if ((local_geo.my_right()[i] * params.ai[i] - params.mesh_off[i]) -
            local_mesh.in_ur[i] <
        ROUND_ERROR_PREC)
      local_mesh.in_ur[i]--;
    if (1.0 + (local_geo.my_left()[i] * params.ai[i] - params.mesh_off[i]) -
            local_mesh.in_ld[i] <
        ROUND_ERROR_PREC)
      local_mesh.in_ld[i]--;
  }
  /* inner grid dimensions */
  for (i = 0; i < 3; i++)
    local_mesh.inner[i] = local_mesh.in_ur[i] - local_mesh.in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for (i = 0; i < 3; i++)
    local_mesh.ld_ind[i] =
        (int)ceil((local_geo.my_left()[i] - full_skin[i]) * params.ai[i] -
                  params.mesh_off[i]);
  /* left down margin */
  for (i = 0; i < 3; i++)
    local_mesh.margin[i * 2] = local_mesh.in_ld[i] - local_mesh.ld_ind[i];
  /* up right grid point */
  for (i = 0; i < 3; i++)
    ind[i] =
        (int)floor((local_geo.my_right()[i] + full_skin[i]) * params.ai[i] -
                   params.mesh_off[i]);
  /* correct roundoff errors at up right boundary */
  for (i = 0; i < 3; i++)
    if (((local_geo.my_right()[i] + full_skin[i]) * params.ai[i] -
         params.mesh_off[i]) -
            ind[i] ==
        0)
      ind[i]--;
  /* up right margin */
  for (i = 0; i < 3; i++)
    local_mesh.margin[(i * 2) + 1] = ind[i] - local_mesh.in_ur[i];

  /* grid dimension */
  local_mesh.size = 1;
  for (i = 0; i < 3; i++) {
    local_mesh.dim[i] = ind[i] - local_mesh.ld_ind[i] + 1;
    local_mesh.size *= local_mesh.dim[i];
  }
  /* reduce inner grid indices from global to local */
  for (i = 0; i < 3; i++)
    local_mesh.in_ld[i] = local_mesh.margin[i * 2];
  for (i = 0; i < 3; i++)
    local_mesh.in_ur[i] = local_mesh.margin[i * 2] + local_mesh.inner[i];

  local_mesh.q_2_off = local_mesh.dim[2] - params.cao;
  local_mesh.q_21_off = local_mesh.dim[2] * (local_mesh.dim[1] - params.cao);
}

void p3m_calc_lm_ld_pos(p3m_local_mesh &local_mesh,
                        const P3MParameters &params) {
  /* spatial position of left down mesh point */
  for (int i = 0; i < 3; i++) {
    local_mesh.ld_pos[i] =
        (local_mesh.ld_ind[i] + params.mesh_off[i]) * params.a[i];
  }
}

#endif /* defined(P3M) || defined(DP3M) */
