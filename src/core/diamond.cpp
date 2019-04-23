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
    This file contains everything needed to create a start-up configuration
    of (partially charged) diamond structure polymer chains with counterions
    and salt molecules, assigning velocities to the particles and
    cross-linking the polymers if necessary.

    The corresponding header file is diamond.hpp.
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
#include "global.hpp"
#include "integrate.hpp"
#include "diamond.hpp"
#include "random.hpp"
#include "utils.hpp"

#include "utils/vec_rotate.hpp"
#include <vector>
using Utils::vec_rotate;


/*************************************************************
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

double mindist4(PartCfg &partCfg, double pos[3]) {
  if (partCfg.size() == 0) {
    return std::min(std::min(box_l[0], box_l[1]), box_l[2]);
  }

  auto const mindist = std::accumulate(
      partCfg.begin(), partCfg.end(), std::numeric_limits<double>::infinity(),
      [&pos](double mindist, Particle const &p) {
        return std::min(mindist, get_mi_vector(pos, p.r.p).norm2());
      });

  if (mindist < std::numeric_limits<double>::infinity())
    return std::sqrt(mindist);
  return -1.0;
}

double buf_mindist4(double const pos[3], int n_add, double const *const add) {
  double mindist = 30000.0, dx, dy, dz;
  int i;

  if (n_add == 0)
    return (std::min(std::min(box_l[0], box_l[1]), box_l[2]));
  for (i = 0; i < n_add; i++) {
    dx = pos[0] - add[3 * i + 0];
    dx -= std::round(dx / box_l[0]) * box_l[0];
    dy = pos[1] - add[3 * i + 1];
    dy -= std::round(dy / box_l[1]) * box_l[1];
    dz = pos[2] - add[3 * i + 2];
    dz -= std::round(dz / box_l[2]) * box_l[2];
    mindist =
        std::min(mindist, Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz));
  }
  if (mindist < 30000.0)
    return (sqrt(mindist));
  return (-1.0);
}

int collision(PartCfg &partCfg, double pos[3], double shield, int n_add,
              double *add) {
  if (mindist4(partCfg, pos) > shield && buf_mindist4(pos, n_add, add) > shield)
    return (0);
  return (1);
}

int constraint_collision(double *p1, double *p2) {
  Vector3d folded_pos1 = folded_position({p1, p1 + 3});
  Vector3d folded_pos2 = folded_position({p2, p2 + 3});

  for (auto &c : Constraints::constraints) {
    auto cs =
        std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
    if (cs) {
      double d1, d2;
      double v[3];

      cs->calc_dist(folded_pos1, &d1, v);
      cs->calc_dist(folded_pos2, &d2, v);

      if (d1 * d2 < 0.0)
        return 1;
    }
  }
  return 0;
}

int counterionsC(PartCfg &partCfg, int N_CI, int part_id, int mode,
                 double shield, int max_try, double val_CI, int type_CI) {
  int n, cnt1, max_cnt;
  double pos[3];

  cnt1 = max_cnt = 0;
  for (n = 0; n < N_CI; n++) {
    for (cnt1 = 0; cnt1 < max_try; cnt1++) {
      pos[0] = box_l[0] * d_random();
      pos[1] = box_l[1] * d_random();
      pos[2] = box_l[2] * d_random();
      if ((mode != 0) || (collision(partCfg, pos, shield, 0, nullptr) == 0))
        break;
      POLY_TRACE(printf("c"); fflush(nullptr));
    }
    if (cnt1 >= max_try)
      return (-1);
    if (place_particle(part_id, pos) == ES_PART_ERROR)
      return (-3);
    set_particle_q(part_id, val_CI);
    set_particle_type(part_id, type_CI);

    part_id++;
    max_cnt = std::max(cnt1, max_cnt);

    POLY_TRACE(printf("C"); fflush(nullptr));
  }
  POLY_TRACE(printf(" %d->%d \n", cnt1, max_cnt));
  if (cnt1 >= max_try)
    return (-1);

  return (std::max(max_cnt, cnt1));
}

int diamondC(PartCfg &partCfg, double a, double bond_length, int MPC, int N_CI,
             double val_nodes, double val_cM, double val_CI, int cM_dist,
             int nonet) {
  int i, j, k, part_id, bond[2], type_bond = 0, type_node = 0, type_cM = 1,
                                 type_nM = 1, type_CI = 2;
  double pos[3], off = bond_length / sqrt(3);
  double dnodes[8][3] = {{0, 0, 0}, {1, 1, 1}, {2, 2, 0}, {0, 2, 2},
                         {2, 0, 2}, {3, 3, 1}, {1, 3, 3}, {3, 1, 3}};
  int dchain[16][5] = {
      {0, 1, +1, +1, +1}, {1, 2, +1, +1, -1}, {1, 3, -1, +1, +1},
      {1, 4, +1, -1, +1}, {2, 5, +1, +1, +1}, {3, 6, +1, +1, +1},
      {4, 7, +1, +1, +1}, {5, 0, +1, +1, -1}, {5, 3, +1, -1, +1},
      {5, 4, -1, +1, +1}, {6, 0, -1, +1, +1}, {6, 2, +1, -1, +1},
      {6, 4, +1, +1, -1}, {7, 0, +1, -1, +1}, {7, 2, -1, +1, +1},
      {7, 3, +1, +1, -1}};

  part_id = 0;
  /* place 8 tetra-functional nodes */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 3; j++) {
      dnodes[i][j] *= a / 4.;
      pos[j] = dnodes[i][j];
    }
    if (place_particle(part_id, pos) == ES_PART_ERROR)
      set_particle_q(part_id, val_nodes);
    set_particle_type(part_id, type_node);

    part_id++;
  }

  /* place intermediate monomers on chains connecting the nodes */
  for (i = 0; i < 2 * 8; i++) {
    for (k = 1; k <= MPC; k++) {
      for (j = 0; j < 3; j++)
        pos[j] = dnodes[dchain[i][0]][j] + k * dchain[i][2 + j] * off;
      if (place_particle(part_id, pos) == ES_PART_ERROR)
        return (-3);
      set_particle_q(part_id, (k % cM_dist == 0) ? val_cM : 0.0);
      set_particle_type(part_id, (k % cM_dist == 0) ? type_cM : type_nM);

      bond[0] = type_bond;
      if (k == 1) {
        if (nonet != 1) {
          bond[1] = dchain[i][0];
          add_particle_bond(part_id, bond);
        }
      } else {
        bond[1] = part_id - 1;
        add_particle_bond(part_id, bond);
      }
      if ((k == MPC) && (nonet != 1)) {
        bond[1] = dchain[i][1];
        add_particle_bond(part_id, bond);
      }
      part_id++;
    }
  }

  /* place counterions (if any) */
  if (N_CI > 0)
    counterionsC(partCfg, N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);

  return (0);
}
