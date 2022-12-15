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

#include "config/config.hpp"

#if defined(P3M) || defined(DP3M)

#include "common.hpp"

#include "LocalBox.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <stdexcept>

double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao) {
  auto const c =
      Utils::sqr(std::cos(Utils::pi() * mesh_i * static_cast<double>(n)));

  switch (cao) {
  case 1: {
    return 1.0;
  }
  case 2: {
    return (1.0 + c * 2.0) / 3.0;
  }
  case 3: {
    return (2.0 + c * (11.0 + c * 2.0)) / 15.0;
  }
  case 4: {
    return (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
  }
  case 5: {
    return (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0)))) /
           2835.0;
  }
  case 6: {
    return (1382.0 +
            c * (35396.0 +
                 c * (83021.0 + c * (34096.0 + c * (2026.0 + c * 4.0))))) /
           155925.0;
  }
  case 7: {
    return (21844.0 +
            c * (776661.0 +
                 c * (2801040.0 +
                      c * (2123860.0 +
                           c * (349500.0 + c * (8166.0 + c * 4.0)))))) /
           6081075.0;
  }
  default: {
    throw std::logic_error("Invalid value cao=" + std::to_string(cao));
  }
  }
}

void P3MLocalMesh::calc_local_ca_mesh(P3MParameters const &params,
                                      LocalBox<double> const &local_geo,
                                      double skin, double space_layer) {
  int i;
  int ind[3];
  // total skin size
  auto const full_skin = Utils::Vector3d{params.cao_cut} +
                         Utils::Vector3d::broadcast(skin) +
                         Utils::Vector3d{0., 0., space_layer};
  // inner left down corner
  auto const &inner_ld_pos = local_geo.my_left();
  // inner up right corner
  auto const &inner_ur_pos = local_geo.my_right();
  // outer up right corner
  auto const outer_ur_pos = inner_ur_pos + full_skin;
  // outer left down corner
  auto const outer_ld_pos = inner_ld_pos - full_skin;
  // convert spatial positions to grid indices
  auto const calc_grid_pos = [&params](Utils::Vector3d const &pos, int i) {
    return pos[i] * params.ai[i] - params.mesh_off[i];
  };

  /* inner left down grid point (global index) */
  for (i = 0; i < 3; i++)
    in_ld[i] = static_cast<int>(std::ceil(calc_grid_pos(inner_ld_pos, i)));
  /* inner up right grid point (global index) */
  for (i = 0; i < 3; i++)
    in_ur[i] = static_cast<int>(std::floor(calc_grid_pos(inner_ur_pos, i)));

  /* correct roundoff errors at boundary */
  for (i = 0; i < 3; i++) {
    if (calc_grid_pos(inner_ur_pos, i) - in_ur[i] < ROUND_ERROR_PREC)
      in_ur[i]--;
    if (calc_grid_pos(inner_ld_pos, i) - in_ld[i] + 1. < ROUND_ERROR_PREC)
      in_ld[i]--;
  }
  /* inner grid dimensions */
  for (i = 0; i < 3; i++)
    inner[i] = in_ur[i] - in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for (i = 0; i < 3; i++)
    ld_ind[i] = static_cast<int>(std::ceil(calc_grid_pos(outer_ld_pos, i)));
  /* left down margin */
  for (i = 0; i < 3; i++)
    margin[i * 2] = in_ld[i] - ld_ind[i];
  /* up right grid point */
  for (i = 0; i < 3; i++)
    ind[i] = static_cast<int>(std::floor(calc_grid_pos(outer_ur_pos, i)));
  /* correct roundoff errors at up right boundary */
  for (i = 0; i < 3; i++)
    if (calc_grid_pos(outer_ur_pos, i) - ind[i] == 0.)
      ind[i]--;
  /* up right margin */
  for (i = 0; i < 3; i++)
    margin[(i * 2) + 1] = ind[i] - in_ur[i];

  /* grid dimension */
  size = 1;
  for (i = 0; i < 3; i++) {
    dim[i] = ind[i] - ld_ind[i] + 1;
    size *= dim[i];
  }

  /* reduce inner grid indices from global to local */
  for (i = 0; i < 3; i++)
    in_ld[i] = margin[i * 2];
  for (i = 0; i < 3; i++)
    in_ur[i] = margin[i * 2] + inner[i];

  q_2_off = dim[2] - params.cao;
  q_21_off = dim[2] * (dim[1] - params.cao);
}

#endif /* defined(P3M) || defined(DP3M) */
