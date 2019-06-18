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
 *  This file contains everything needed to create a start-up configuration
 *  of (partially charged) diamond structure polymer chains with counterions
 *  and salt molecules, assigning velocities to the particles and
 *  cross-linking the polymers if necessary.
 *
 *  The corresponding header file is diamond.hpp.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "PartCfg.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"
#include "debug.hpp"
#include "diamond.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "polymer.hpp"
#include "random.hpp"

int create_counterions(PartCfg &partCfg, int const N_CI, int part_id,
                       int const mode, double const shield, int const max_try,
                       double const val_CI, int const type_CI) {
  Utils::Vector3d pos;

  int cnt1 = 0;
  int max_cnt = 0;
  for (int n = 0; n < N_CI; ++n) {
    for (cnt1 = 0; cnt1 < max_try; ++cnt1) {
      pos[0] = box_l[0] * d_random();
      pos[1] = box_l[1] * d_random();
      pos[2] = box_l[2] * d_random();
      if ((mode != 0) or (mindist(partCfg, pos) > shield))
        break;
    }
    if (cnt1 >= max_try)
      throw std::runtime_error(
          "Too many failed attempts finding valid position.");
    if (place_particle(part_id, pos.data()) == ES_PART_ERROR)
      throw std::runtime_error("Failed to place particle.");
    set_particle_q(part_id, val_CI);
    set_particle_type(part_id, type_CI);

    part_id++;
    max_cnt = std::max(cnt1, max_cnt);
  }
  POLY_TRACE(printf(" %d->%d \n", cnt1, max_cnt));

  return max_cnt;
}

int create_diamond(PartCfg &partCfg, double const a, double const bond_length,
                   int const MPC, int const N_CI, double const val_nodes,
                   double const val_cM, double const val_CI, int const cM_dist,
                   int const nonet) {
  Utils::Vector3d pos;
  double const off = bond_length / sqrt(3);
  auto const dnodes =
      a / 4. * Utils::Vector<Utils::Vector3d, 8>{{0, 0, 0}, {1, 1, 1},
                                                 {2, 2, 0}, {0, 2, 2},
                                                 {2, 0, 2}, {3, 3, 1},
                                                 {1, 3, 3}, {3, 1, 3}};
  static constexpr int dchain[16][5] = {
      {0, 1, +1, +1, +1}, {1, 2, +1, +1, -1}, {1, 3, -1, +1, +1},
      {1, 4, +1, -1, +1}, {2, 5, +1, +1, +1}, {3, 6, +1, +1, +1},
      {4, 7, +1, +1, +1}, {5, 0, +1, +1, -1}, {5, 3, +1, -1, +1},
      {5, 4, -1, +1, +1}, {6, 0, -1, +1, +1}, {6, 2, +1, -1, +1},
      {6, 4, +1, +1, -1}, {7, 0, +1, -1, +1}, {7, 2, -1, +1, +1},
      {7, 3, +1, +1, -1}};
  int const type_bond = 0;
  int const type_node = 0;
  int const type_cM = 1;
  int const type_nM = 1;
  int const type_CI = 2;

  int part_id = 0;
  /* place 8 tetra-functional nodes */
  for (auto &dnode : dnodes) {
    for (int i = 0; i < 3; ++i) {
      pos[i] = dnode[i];
    }
    if (place_particle(part_id, pos.data()) == ES_PART_ERROR)
      return -3;
    set_particle_q(part_id, val_nodes);
    set_particle_type(part_id, type_node);

    part_id++;
  }

  /* place intermediate monomers on chains connecting the nodes */
  for (auto const &chain : dchain) {
    for (int k = 1; k <= MPC; ++k) {
      for (int j = 0; j < 3; ++j)
        pos[j] = dnodes[chain[0]][j] + k * chain[2 + j] * off;
      if (place_particle(part_id, pos.data()) == ES_PART_ERROR)
        throw std::runtime_error("Failed to place particle.");
      set_particle_q(part_id, (k % cM_dist == 0) ? val_cM : 0.0);
      set_particle_type(part_id, (k % cM_dist == 0) ? type_cM : type_nM);

      int bond[2];
      bond[0] = type_bond;
      if (k == 1) {
        if (nonet != 1) {
          bond[1] = chain[0];
          add_particle_bond(part_id, bond);
        }
      } else {
        bond[1] = part_id - 1;
        add_particle_bond(part_id, bond);
      }
      if ((k == MPC) && (nonet != 1)) {
        bond[1] = chain[1];
        add_particle_bond(part_id, bond);
      }
      part_id++;
    }
  }

  /* place counterions (if any) */
  if (N_CI > 0)
    create_counterions(partCfg, N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);

  return 0;
}
