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
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and cross-linking the polymers if
   necessary.

    The corresponding header file is polymer.hpp.
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
#include "polymer.hpp"
#include "random.hpp"
#include "utils.hpp"

#include "utils/vec_rotate.hpp"
using Utils::vec_rotate;


bool is_valid_position(const Utils::Vector3d *pos, const std::vector<std::vector<Utils::Vector3d>> *positions, PartCfg &partCfg, const double shield, const int constr) {
    // check if constraint is violated
    if (constr == true) {
        Utils::Vector3d folded_pos = folded_position(*pos);
        for (auto &c : Constraints::constraints) {
            auto cs = std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
            if (cs) {
                double d;
                double v[3];

                cs->calc_dist(folded_pos, &d, v);

                if (d <= 0) {
//                    std::cout << "CONSTRAINT violated." << std::endl;
                    return false;
                }
            }
        }
    }

    if (shield > 0) {
        // check for collision with existing particles
        if (mindist5(partCfg, *pos) < shield) {
//                    std::cout << "EXISTING violated." << std::endl;
            return false;
        }
        // check for collision with buffered positions
        double buff_mindist = std::numeric_limits<double>::infinity();
        double h;
        for (auto p = positions->begin(); p != positions->end(); ++p) {
            for (auto m = p->begin(); m != p->end(); ++m) {
                h = (*pos - *m).norm2();
                buff_mindist = std::min(h, buff_mindist);
            }
        }
        if (std::sqrt(buff_mindist) < shield) {
//                    std::cout << "BUFFERED violated." << std::endl;
            return false;
        }
    }
    return true;
}

Utils::Vector3d random_position() {
    Utils::Vector3d v;
    for (int i=0; i<3; ++i)
        v[i] = box_l[i] * d_random();
    return v;
}

Utils::Vector3d random_unit_vector() {
    Utils::Vector3d v;
    const int phi = acos( 1. - 2. * d_random() );
    const int theta = 2. * Utils::pi() * d_random();
    v[0] = sin(phi) * cos(theta);
    v[1] = sin(phi) * sin(theta);
    v[2] = cos(phi);
    v /= v.norm();
    return v;
}

std::vector<std::vector<Utils::Vector3d>> draw_polymer_positions(PartCfg &partCfg, int N_P, int MPC, double bond_length,
              std::vector<Utils::Vector3d> &start_positions, int mode, double shield, int max_try, const int use_bond_angle,
              double angle, double angle2, int constr) {
    std::vector<std::vector<Utils::Vector3d>> positions(N_P, std::vector<Utils::Vector3d>(MPC));

    // make sure that if given, all starting positions are valid
    if ((not start_positions.empty())
            and std::any_of(start_positions.begin(), start_positions.end(),
                [&positions, &partCfg, shield, constr](Utils::Vector3d v)
                { return not is_valid_position(&v, &positions, partCfg, shield, constr); }))
        throw std::runtime_error("Invalid starting positions.");
    // else generate initial positions
    for (int p = 0; p < N_P; ++p) {
        int counter = 0;
        Utils::Vector3d trial_pos;
        // first monomer
        if (start_positions.empty()) {
            do {
                trial_pos = random_position();
                counter++;
            } while ((not is_valid_position(&trial_pos, &positions, partCfg, shield, constr))
                    and (counter < max_try));
            if (counter == max_try) {
                throw std::runtime_error("Failed to create polymer positions.");
            } else {
                positions[p][0] = trial_pos;
            }
        } else {
            positions[p][0] = start_positions[p];
        }
    }
    // remaining monomers
    for (int p = 0; p < N_P; ++p) {
        int counter = 0;
        Utils::Vector3d trial_pos;
        for (int m = 1; m < MPC; ++m) {
            do {
                if (not use_bond_angle or m < 2) {
                    // random step
                    trial_pos = positions[p][m-1] + bond_length * random_unit_vector();
                } else {
                    // use prescribed angle
                    Utils::Vector3d last_vec = positions[p][m-1] - positions[p][m-2];
                    trial_pos = positions[p][m-1] + vec_rotate(
                            vector_product(last_vec, random_unit_vector()),
                            angle,
                            last_vec);
                }
                counter++;
            } while ((not is_valid_position(&trial_pos, &positions, partCfg, shield, constr))
                    and (counter < max_try));
            if (counter == max_try) {
                throw std::runtime_error("22 Failed to create polymer positions.");
            }
            positions[p][m] = trial_pos;
            counter = 0;
        }
    }
    return positions;
}

/*************************************************************
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

double mindist5(PartCfg &partCfg, const Utils::Vector3d pos) {
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

int mindist3(PartCfg &partCfg, int part_id, double r_catch, int *ids) {
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
  Utils::Vector3d folded_pos1 = folded_position({p1, p1 + 3});
  Utils::Vector3d folded_pos2 = folded_position({p2, p2 + 3});

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

int polymerC(PartCfg &partCfg, int N_P, int MPC, double bond_length,
             int part_id, double const *posed, int mode, double shield,
             int max_try, double val_cM, int cM_dist, int type_nM, int type_cM,
             int type_bond, double angle, double angle2, double const *posed2,
             int constr) {
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
    if (p > 0)
      posed = nullptr;

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
          phi = 2.0 * Utils::pi() * d_random();
          pos[0] = poz[0] + rr * cos(phi);
          pos[1] = poz[1] + rr * sin(phi);
          pos[2] = poz[2] + zz;
          if (constr == 0 ||
              constraint_collision(pos, poly.data() + 3 * (n - 1)) == 0) {

            if (mode == 1 ||
                collision(partCfg, pos, shield, n, poly.data()) == 0)
              break;
            if (mode == 0) {
              cnt1 = -1;
              break;
            }
          }
          POLY_TRACE(printf("m"); fflush(nullptr));
        }
        if (cnt1 >= max_try) {
          fprintf(stderr,
                  "\nWarning! Attempt #%d to build polymer %d failed "
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
              phi = 2.0 * Utils::pi() * d_random();
              vec_rotate(a, phi, c, d);
            }

            vec_rotate(d, angle, a, b);

            pos[0] = poz[0] + b[0];
            pos[1] = poz[1] + b[1];
            pos[2] = poz[2] + b[2];

          } else {
            zz = (2.0 * d_random() - 1.0) * bond_length;
            rr = sqrt(Utils::sqr(bond_length) - Utils::sqr(zz));
            phi = 2.0 * Utils::pi() * d_random();
            pos[0] = poz[0] + rr * cos(phi);
            pos[1] = poz[1] + rr * sin(phi);
            pos[2] = poz[2] + zz;
          }

          // POLY_TRACE(/* printf("a=(%f,%f,%f) absa=%f M=(%f,%f,%f)
          // c=(%f,%f,%f) absMc=%f
          // a*c=%f)\n",a[0],a[1],a[2],sqrt(Utils::sqr(a[0])+Utils::sqr(a[1])+Utils::sqr(a[2])),M[0],M[1],M[2],c[0],c[1],c[2],sqrt(Utils::sqr(M[0]+c[0])+Utils::sqr(M[1]+c[1])+Utils::sqr(M[2]+c[2])),a[0]*c[0]+a[1]*c[1]+a[2]*c[2])
          // */);
          // POLY_TRACE(/* printf("placed Monomer %d at
          // (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);

          if (constr == 0 ||
              constraint_collision(pos, poly.data() + 3 * (n - 1)) == 0) {
            if (mode == 1 ||
                collision(partCfg, pos, shield, n, poly.data()) == 0)
              break;
            if (mode == 0) {
              cnt1 = -2;
              break;
            }
          }
          POLY_TRACE(printf("m"); fflush(nullptr));
        }
        if (cnt1 >= max_try) {
          fprintf(stderr,
                  "\nWarning! Attempt #%d to build polymer %d failed "
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
    }

    max_cnt = std::max(max_cnt, std::max(cnt1, cnt2));

    /* actually creating current polymer in ESPResSo */
    for (n = 0; n < MPC; n++) {

      pos[0] = poly[3 * n];
      pos[1] = poly[3 * n + 1];
      pos[2] = poly[3 * n + 2];
      if (place_particle(part_id, pos) == ES_PART_ERROR)
        return -3;

      set_particle_q(part_id, ((n % cM_dist == 0) ? val_cM : 0.0));
      set_particle_type(part_id, ((n % cM_dist == 0) ? type_cM : type_nM));

      if (n >= bond_size) {
        bond[1] = part_id - bond_size;
        for (i = 2; i <= bond_size; i++) {
          bond[i] = part_id - bond_size + i;
        }
        add_particle_bond(part_id - bond_size + 1, bond);
      }
      part_id++;
      // POLY_TRACE(/* printf("placed Monomer %d at
      // (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);
    }
  }

  return (std::max(max_cnt, cnt2));
}
