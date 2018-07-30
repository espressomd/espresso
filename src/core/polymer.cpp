/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file polymer.cpp
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and crosslinking the polymers if
   necessary.

    The corresponding header file is polymer.hpp.
*/

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"
#include "debug.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "PartCfg.hpp"
#include "polymer.hpp"
#include "random.hpp"
#include "utils.hpp"

/*************************************************************
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

int mindist3(PartCfg & partCfg, int part_id, double r_catch, int *ids) {
  int caught = 0;

  auto const r_catch2 = r_catch * r_catch;
  auto const &part = partCfg[part_id];

  for (auto const &p : partCfg) {
    if (p != part) {
      if (get_mi_vector(part.r.p, p.r.p).norm2() < r_catch2)
        ids[caught++] = p.p.identity;
    }
  }

  return caught;
}

double mindist4(PartCfg & partCfg, double pos[3]) {
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

double buf_mindist4(double pos[3], int n_add, double *add) {
  double mindist = 30000.0, dx, dy, dz;
  int i;

  if (n_add == 0)
    return (std::min(std::min(box_l[0], box_l[1]), box_l[2]));
  for (i = 0; i < n_add; i++) {
    dx = pos[0] - add[3 * i + 0];
    dx -= dround(dx / box_l[0]) * box_l[0];
    dy = pos[1] - add[3 * i + 1];
    dy -= dround(dy / box_l[1]) * box_l[1];
    dz = pos[2] - add[3 * i + 2];
    dz -= dround(dz / box_l[2]) * box_l[2];
    mindist = std::min(mindist, Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz));
  }
  if (mindist < 30000.0)
    return (sqrt(mindist));
  return (-1.0);
}

int collision(PartCfg & partCfg, double pos[3], double shield, int n_add, double *add) {
  if (mindist4(partCfg, pos) > shield && buf_mindist4(pos, n_add, add) > shield)
    return (0);
  return (1);
}

int constraint_collision(double *p1, double *p2) {
  Particle part1, part2;
  double d1, d2, v[3];
  double folded_pos1[3];
  double folded_pos2[3];
  int img[3];

  memmove(folded_pos1, p1, 3 * sizeof(double));
  fold_position(folded_pos1, img);

  memmove(folded_pos2, p2, 3 * sizeof(double));
  fold_position(folded_pos2, img);

  for (auto &c : Constraints::constraints) {
    auto cs =
        std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
    if (cs) {
      cs->calc_dist(folded_pos1, &d1, v);
      cs->calc_dist(folded_pos2, &d2, v);

      if (d1 * d2 < 0.0)
        return 1;
    }
  }
  return 0;
}

int polymerC(PartCfg & partCfg, int N_P, int MPC, double bond_length, int part_id, double *posed,
             int mode, double shield, int max_try, double val_cM, int cM_dist,
             int type_nM, int type_cM, int type_bond, double angle,
             double angle2, double *posed2, int constr) {
  int p, n, cnt1, cnt2, max_cnt, bond_size, i;
  double phi, zz, rr;
  double pos[3];
  double poz[3];
  double poy[3] = {0, 0, 0};
  double pox[3] = {0, 0, 0};
  double a[3] = {0, 0, 0};
  double b[3], c[3] = {0., 0., 0.}, d[3];
  double absc;

  std::vector<double> poly(3 * MPC);

  bond_size = bonded_ia_params[type_bond].num;
  std::vector<int> bond(bond_size + 1);
  bond[0] = type_bond;

  cnt1 = cnt2 = max_cnt = 0;
  for (p = 0; p < N_P; p++) {
    if (p > 0) posed = nullptr;

    for (cnt2 = 0; cnt2 < max_try; cnt2++) {
      /* place start monomer */
      if (posed != nullptr) {
        /* if position of 1st monomer is given */
        pos[0] = posed[0];
        pos[1] = posed[1];
        pos[2] = posed[2];
      } else {
        /* randomly set position */
        for (cnt1 = 0; cnt1 < max_try; cnt1++) {
          pos[0] = box_l[0] * d_random();
          pos[1] = box_l[1] * d_random();
          pos[2] = box_l[2] * d_random();
          if ((mode == 1) || (collision(partCfg, pos, shield, 0, nullptr) == 0))
            break;
          POLY_TRACE(printf("s"); fflush(nullptr));
        }
        if (cnt1 >= max_try) {
          return (-1);
        }
      }
      poly[0] = pos[0];
      poly[1] = pos[1];
      poly[2] = pos[2];

      max_cnt = std::max(cnt1, max_cnt);
      POLY_TRACE(printf("S"); fflush(nullptr));

      poz[0] = pos[0];
      poz[1] = pos[1];
      poz[2] = pos[2];

      /* place 2nd monomer */
      n = 1;
      if (posed2 != nullptr && posed != nullptr && angle2 > -1.0) {
        /* if position of 2nd monomer is given */
        pos[0] = posed2[0];
        pos[1] = posed2[1];
        pos[2] = posed2[2];
        /* calculate preceding monomer so that bond_length is correct */
        absc = sqrt(Utils::sqr(pos[0] - poz[0]) + Utils::sqr(pos[1] - poz[1]) +
                    Utils::sqr(pos[2] - poz[2]));
        poz[0] = pos[0] + (poz[0] - pos[0]) * bond_length / absc;
        poz[1] = pos[1] + (poz[1] - pos[1]) * bond_length / absc;
        poz[2] = pos[2] + (poz[2] - pos[2]) * bond_length / absc;
        // POLY_TRACE(/* printf("virtually shifted position of first monomer to
        // (%f,%f,%f)\n",poz[0],poz[1],poz[2]) */);
      } else {
        /* randomly place 2nd monomer */
        for (cnt1 = 0; cnt1 < max_try; cnt1++) {
          zz = (2.0 * d_random() - 1.0) * bond_length;
          rr = sqrt(Utils::sqr(bond_length) - Utils::sqr(zz));
          phi = 2.0 * PI * d_random();
          pos[0] = poz[0] + rr * cos(phi);
          pos[1] = poz[1] + rr * sin(phi);
          pos[2] = poz[2] + zz;
          if (constr == 0 ||
              constraint_collision(pos, poly.data() + 3 * (n - 1)) == 0) {

            if (mode == 1 || collision(partCfg, pos, shield, n, poly.data()) == 0)
              break;
            if (mode == 0) {
              cnt1 = -1;
              break;
            }
          }
          POLY_TRACE(printf("m"); fflush(nullptr));
        }
        if (cnt1 >= max_try) {
          fprintf(stderr, "\nWarning! Attempt #%d to build polymer %d failed "
                          "while placing monomer 2!\n",
                  cnt2 + 1, p);
          fprintf(stderr, "         Retrying by re-setting the start-monomer "
                          "of current chain...\n");
        }
        if (cnt1 == -1 || cnt1 >= max_try) {
          continue; /* continue the main loop */
        }
      }
      if (posed2 != nullptr && p > 0) {
        posed2 = nullptr;
      }
      poly[3 * n] = pos[0];
      poly[3 * n + 1] = pos[1];
      poly[3 * n + 2] = pos[2];

      max_cnt = std::max(cnt1, max_cnt);
      POLY_TRACE(printf("M"); fflush(nullptr));

      /* place remaining monomers */
      for (n = 2; n < MPC; n++) {
        if (angle2 > -1.0) {
          if (n == 2) { /* if the 2nd angle is set, construct preceding monomer
                          with resulting plane perpendicular on the xy-plane */
            poy[0] = 2 * poz[0] - pos[0];
            poy[1] = 2 * poz[1] - pos[1];
            if (pos[2] == poz[2])
              poy[2] = poz[2] + 1;
            else
              poy[2] = poz[2];
          } else {
            /* save 3rd last monomer */
            pox[0] = poy[0];
            pox[1] = poy[1];
            pox[2] = poy[2];
          }
        }
        if (angle > -1.0) {
          /* save one but last monomer */
          poy[0] = poz[0];
          poy[1] = poz[1];
          poy[2] = poz[2];
        }
        /* save last monomer */
        poz[0] = pos[0];
        poz[1] = pos[1];
        poz[2] = pos[2];

        if (angle > -1.0) {
          a[0] = poy[0] - poz[0];
          a[1] = poy[1] - poz[1];
          a[2] = poy[2] - poz[2];

          b[0] = pox[0] - poy[0];
          b[1] = pox[1] - poy[1];
          b[2] = pox[2] - poy[2];

          vector_product(a, b, c);
        }

        for (cnt1 = 0; cnt1 < max_try; cnt1++) {
          if (angle > -1.0) {
            if (sqrlen(c) < ROUND_ERROR_PREC) {
              fprintf(stderr, "WARNING: rotation axis is 0,0,0, check the "
                              "angles given to the polymer command\n");
              c[0] = 1;
              c[1] = 0;
              c[2] = 0;
            }
            if (angle2 > -1.0 && n > 2) {
              vec_rotate(a, angle2, c, d);
            } else {
              phi = 2.0 * PI * d_random();
              vec_rotate(a, phi, c, d);
            }

            vec_rotate(d, angle, a, b);

            pos[0] = poz[0] + b[0];
            pos[1] = poz[1] + b[1];
            pos[2] = poz[2] + b[2];

          } else {
            zz = (2.0 * d_random() - 1.0) * bond_length;
            rr = sqrt(Utils::sqr(bond_length) - Utils::sqr(zz));
            phi = 2.0 * PI * d_random();
            pos[0] = poz[0] + rr * cos(phi);
            pos[1] = poz[1] + rr * sin(phi);
            pos[2] = poz[2] + zz;
          }

// POLY_TRACE(/* printf("a=(%f,%f,%f) absa=%f M=(%f,%f,%f) c=(%f,%f,%f) absMc=%f
// a*c=%f)\n",a[0],a[1],a[2],sqrt(Utils::sqr(a[0])+Utils::sqr(a[1])+Utils::sqr(a[2])),M[0],M[1],M[2],c[0],c[1],c[2],sqrt(Utils::sqr(M[0]+c[0])+Utils::sqr(M[1]+c[1])+Utils::sqr(M[2]+c[2])),a[0]*c[0]+a[1]*c[1]+a[2]*c[2])
// */);
// POLY_TRACE(/* printf("placed Monomer %d at
// (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);

          if (constr == 0 ||
              constraint_collision(pos, poly.data() + 3 * (n - 1)) == 0) {
            if (mode == 1 || collision(partCfg, pos, shield, n, poly.data()) == 0)
              break;
            if (mode == 0) {
              cnt1 = -2;
              break;
            }
          }
          POLY_TRACE(printf("m"); fflush(nullptr));
        }
        if (cnt1 >= max_try) {
          fprintf(stderr, "\nWarning! Attempt #%d to build polymer %d failed "
                          "after %d unsuccessful trials to place monomer %d!\n",
                  cnt2 + 1, p, cnt1, n);
          fprintf(stderr, "         Retrying by re-setting the start-monomer "
                          "of current chain...\n");
        }
        if (cnt1 == -2 || cnt1 >= max_try) {
          n = 0;
          break;
        }
        poly[3 * n] = pos[0];
        poly[3 * n + 1] = pos[1];
        poly[3 * n + 2] = pos[2];

        max_cnt = std::max(cnt1, max_cnt);

        POLY_TRACE(printf("M"); fflush(nullptr));
      }
      if (n > 0)
        break;
    } /* cnt2 */
    POLY_TRACE(printf(" %d/%d->%d \n", cnt1, cnt2, max_cnt));
    if (cnt2 >= max_try) {
      return (-2);
    } else

      max_cnt = std::max(max_cnt, std::max(cnt1, cnt2));

    /* actually creating current polymer in ESPResSo */
    for (n = 0; n < MPC; n++) {

      pos[0] = poly[3 * n];
      pos[1] = poly[3 * n + 1];
      pos[2] = poly[3 * n + 2];
      if (place_particle(part_id, pos) == ES_PART_ERROR ||
          (set_particle_q(part_id, ((n % cM_dist == 0) ? val_cM : 0.0)) ==
           ES_ERROR) ||
          (set_particle_type(part_id,
                             ((n % cM_dist == 0) ? type_cM : type_nM)) ==
           ES_ERROR)) {
        return (-3);
      }

      if (n >= bond_size) {
        bond[1] = part_id - bond_size;
        for (i = 2; i <= bond_size; i++) {
          bond[i] = part_id - bond_size + i;
        }
        if (change_particle_bond(part_id - bond_size + 1, bond.data(), 0) ==
            ES_ERROR) {
          return (-3);
        }
      }
      part_id++;
      // POLY_TRACE(/* printf("placed Monomer %d at
      // (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);
    }
  }

  return (std::max(max_cnt, cnt2));
}

int counterionsC(PartCfg & partCfg, int N_CI, int part_id, int mode, double shield, int max_try,
                 double val_CI, int type_CI) {
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
    if (set_particle_q(part_id, val_CI) == ES_ERROR)
      return (-3);
    if (set_particle_type(part_id, type_CI) == ES_ERROR)
      return (-3);
    part_id++;
    max_cnt = std::max(cnt1, max_cnt);

    POLY_TRACE(printf("C"); fflush(nullptr));
  }
  POLY_TRACE(printf(" %d->%d \n", cnt1, max_cnt));
  if (cnt1 >= max_try)
    return (-1);

  return (std::max(max_cnt, cnt1));
}






int diamondC(PartCfg & partCfg, double a, double bond_length, int MPC, int N_CI, double val_nodes,
             double val_cM, double val_CI, int cM_dist, int nonet) {
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
      return (-3);
    if (set_particle_q(part_id, val_nodes) == ES_ERROR)
      return (-3);
    if (set_particle_type(part_id, type_node) == ES_ERROR)
      return (-3);
    part_id++;
  }

  /* place intermediate monomers on chains connecting the nodes */
  for (i = 0; i < 2 * 8; i++) {
    for (k = 1; k <= MPC; k++) {
      for (j = 0; j < 3; j++)
        pos[j] = dnodes[dchain[i][0]][j] + k * dchain[i][2 + j] * off;
      if (place_particle(part_id, pos) == ES_PART_ERROR)
        return (-3);
      if (set_particle_q(part_id, (k % cM_dist == 0) ? val_cM : 0.0) ==
          ES_ERROR)
        return (-3);
      if (set_particle_type(part_id, (k % cM_dist == 0) ? type_cM : type_nM) ==
          ES_ERROR)
        return (-3);
      bond[0] = type_bond;
      if (k == 1) {
        if (nonet != 1) {
          bond[1] = dchain[i][0];
          if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
            return (-3);
        }
      } else {
        bond[1] = part_id - 1;
        if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
          return (-3);
      }
      if ((k == MPC) && (nonet != 1)) {
        bond[1] = dchain[i][1];
        if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
          return (-3);
      }
      part_id++;
    }
  }

  /* place counterions (if any) */
  if (N_CI > 0)
    counterionsC(partCfg, N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);

  return (0);
}

int icosaederC(PartCfg & partCfg, double ico_a, int MPC, int N_CI, double val_cM, double val_CI,
               int cM_dist) {
  int i, j, k, l, part_id, bond[2], type_bond = 0, type_cM = 0, type_nM = 1,
                                    type_CI = 2;
  double pos[3], pos_shift[3], vec[3], e_vec[3], vec_l,
      bond_length = (2 * ico_a / 3.) / (1. * MPC);
  double ico_g = ico_a * (1 + sqrt(5)) / 2.0, shift = 0.0;
  double ico_coord[12][3] = {
      {0, +ico_a, +ico_g}, {0, +ico_a, -ico_g}, {0, -ico_a, +ico_g},
      {0, -ico_a, -ico_g}, {+ico_a, +ico_g, 0}, {+ico_a, -ico_g, 0},
      {-ico_a, +ico_g, 0}, {-ico_a, -ico_g, 0}, {+ico_g, 0, +ico_a},
      {-ico_g, 0, +ico_a}, {+ico_g, 0, -ico_a}, {-ico_g, 0, -ico_a}};
  int ico_NN[12][5] = {{2, 8, 4, 6, 9},   {3, 10, 4, 6, 11}, {0, 8, 5, 7, 9},
                       {1, 10, 5, 7, 11}, {0, 6, 1, 10, 8},  {2, 7, 3, 10, 8},
                       {0, 4, 1, 11, 9},  {2, 5, 3, 11, 9},  {0, 2, 5, 10, 4},
                       {0, 2, 7, 11, 6},  {1, 3, 5, 8, 4},   {1, 3, 7, 9, 6}};
  int ico_ind[12][10];

  /* make sure that the edges in ico_NN are sorted such that NearestNeighbours
   * are next to each other */
  /* int    ico_NN[12][5]    = {{2,4,6,8, 9}, {3,4,6,10,11}, {0,5,7,8, 9},
  {1,5,7,10,11}, {0,1,6,8,10}, {2,3,7,8,10},
                             {0,1,4,9,11}, {2,3,5, 9,11}, {0,2,4,5,10}, {0,2,6,
  7,11}, {1,3,4,5, 8}, {1,3,6,7, 9}};
  for(i=0; i<12; i++) {
    printf("%d: { ",i);
    for(j=0; j<5; j++) printf("%d ",ico_NN[i][j]);
    printf("} -> ");
    for(j=0; j<5; j++)
      for(l=0; l<5; l++)
        if(j!=l)
          for(k=0; k<5; k++)
            if(ico_NN[ico_NN[i][j]][k]==ico_NN[i][l]) printf("%d = %d (@%d)
  ",ico_NN[i][j],ico_NN[i][l],k);
    printf("\n");
  } */

  /* shift coordinates to not be centered around zero but rather be positive */
  if (ico_a > ico_g)
    shift = ico_a;
  else
    shift = ico_g;

  /* create fulleren & soccer-ball */
  part_id = 0;
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 5; j++) {
      /* place chains along the 5 edges around each of the 12 icosaeder's
       * vertices */
      if (j < 4)
        for (l = 0; l < 3; l++)
          vec[l] =
              (ico_coord[ico_NN[i][j + 1]][l] - ico_coord[ico_NN[i][j]][l]) /
              3.;
      else
        for (l = 0; l < 3; l++)
          vec[l] =
              (ico_coord[ico_NN[i][0]][l] - ico_coord[ico_NN[i][4]][l]) / 3.;
      vec_l = sqrt(Utils::sqr(vec[0]) + Utils::sqr(vec[1]) + Utils::sqr(vec[2]));
      for (l = 0; l < 3; l++)
        e_vec[l] = vec[l] / vec_l;

      ico_ind[i][j] = part_id;
      bond_length = vec_l / (1. * MPC);
      for (l = 0; l < 3; l++)
        pos[l] = ico_coord[i][l] +
                 (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l]) / 3.;
      for (k = 0; k < MPC; k++) {
        for (l = 0; l < 3; l++)
          pos_shift[l] = pos[l] + shift;
        if (place_particle(part_id, pos_shift) == ES_PART_ERROR)
          return (-3);
        if (set_particle_q(part_id, val_cM) == ES_ERROR)
          return (-3);
        if (set_particle_type(part_id, type_cM) == ES_ERROR)
          return (-3);
        bond[0] = type_bond;
        if (k > 0) {
          bond[1] = part_id - 1;
          if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
            return (-3);
        }
        part_id++;
        for (l = 0; l < 3; l++)
          pos[l] += bond_length * e_vec[l];
      }

      /* place chains along the 5 edges on the middle third of the connection
       * between two NN vertices */
      if (i < ico_NN[i][j]) {
        for (l = 0; l < 3; l++)
          vec[l] = (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l]) / 3.;
        vec_l = sqrt(Utils::sqr(vec[0]) + Utils::sqr(vec[1]) + Utils::sqr(vec[2]));
        for (l = 0; l < 3; l++)
          e_vec[l] = vec[l] / vec_l;

        ico_ind[i][j + 5] = part_id;
        bond_length = vec_l / (1. * MPC);
        for (l = 0; l < 3; l++)
          pos[l] = ico_coord[i][l] +
                   (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l]) / 3. +
                   bond_length * e_vec[l];
        for (k = 1; k < MPC; k++) {
          for (l = 0; l < 3; l++)
            pos_shift[l] = pos[l] + shift;
          if (place_particle(part_id, pos_shift) == ES_ERROR)
            return (-3);
          if (set_particle_q(part_id, 0.0) == ES_ERROR)
            return (-3);
          if (set_particle_type(part_id, type_nM) == ES_ERROR)
            return (-3);
          bond[0] = type_bond;
          if (k > 1) {
            bond[1] = part_id - 1;
            if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
              return (-3);
          } else {
            bond[1] = ico_ind[i][j];
            if (change_particle_bond(part_id, bond, 0) == ES_ERROR)
              return (-3);
          }
          part_id++;
          for (l = 0; l < 3; l++)
            pos[l] += bond_length * e_vec[l];
        }
      }
    }

    for (j = 0; j < 5; j++) {
      /* add bonds between the edges around the vertices */
      bond[0] = type_bond;
      //      if(j>0) bond[1] = ico_ind[i][j-1] + (MPC-1); else bond[1] =
      //      ico_ind[i][4] + (MPC-1);
      if (j > 0)
        bond[1] = ico_ind[i][j - 1] + (MPC - 1);
      else if (MPC > 0)
        bond[1] = ico_ind[i][4] + (MPC - 1);
      else
        bond[1] = ico_ind[i][4];
      if (change_particle_bond(ico_ind[i][j], bond, 0) == ES_ERROR)
        return (-2);

      /* connect loose edges around vertices with chains along the middle third
       * already created earlier */
      if (i > ico_NN[i][j]) {
        bond[0] = type_bond;
        for (l = 0; l < 5; l++)
          if (ico_NN[ico_NN[i][j]][l] == i)
            break;
        if (l == 5) {
          fprintf(stderr, "INTERNAL ERROR: Couldn't find my neighbouring edge "
                          "upon creating the icosaeder!\n");
          errexit();
        }
        bond[1] = ico_ind[ico_NN[i][j]][l + 5] + (MPC - 2);
        if (change_particle_bond(ico_ind[i][j], bond, 0) == ES_ERROR)
          return (-1);
      }
    }
  }

  /* place counterions (if any) */
  if (N_CI > 0)
    counterionsC(partCfg, N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);

  return (0);
}
