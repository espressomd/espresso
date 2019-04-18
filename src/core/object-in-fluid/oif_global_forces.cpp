/*
  Copyright (C) 2012-2018 The ESPResSo project

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
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
 */

#include "oif_global_forces.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "particle_data.hpp"

#include "utils/math/triangle_functions.hpp"
using Utils::angle_btw_triangles;
using Utils::area_triangle;
using Utils::get_n_triangle;

/** set parameters for the OIF_GLOBAL_FORCES potential.
 */
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g,
                                 double V0, double kv) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.oif_global_forces.ka_g = ka_g;
  bonded_ia_params[bond_type].p.oif_global_forces.A0_g = A0_g;
  bonded_ia_params[bond_type].p.oif_global_forces.V0 = V0;
  bonded_ia_params[bond_type].p.oif_global_forces.kv = kv;

  bonded_ia_params[bond_type].type = BONDED_IA_OIF_GLOBAL_FORCES;
  bonded_ia_params[bond_type].num = 2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

/** called in force_calc() from within forces.cpp
 *  calculates the global area and global volume for a cell before the forces
 * are handled
 *  sums up parts for area with mpi_reduce from local triangles
 *  synchronization with allreduce
 *
 *  !!! loop over particles from domain_decomposition !!!
 */

void calc_oif_global(double *area_volume,
                     int molType) { // first-fold-then-the-same approach
  double partArea = 0.0;
  double part_area_volume[2]; // added

  // z volume
  double VOL_partVol = 0.;

  /* loop over particles */
  Particle *p1, *p2, *p3;
  Utils::Vector3d p11, p22, p33;
  int img[3];
  double AA[3], BB[3];
  Bonded_ia_parameters *iaparams;
  int type_num, n_partners, id;
  BondedInteraction type;

  int test = 0;

  for (auto &p : local_cells.particles()) {
    int j = 0;
    p1 = &p;
    while (j < p1->bl.n) {
      /* bond type */
      type_num = p1->bl.e[j++];
      iaparams = &bonded_ia_params[type_num];
      type = iaparams->type;
      n_partners = iaparams->num;
      id = p1->p.mol_id;
      if (type == BONDED_IA_OIF_GLOBAL_FORCES &&
          id == molType) { // BONDED_IA_OIF_GLOBAL_FORCES with correct molType
        test++;
        /* fetch particle 2 */
        p2 = local_particles[p1->bl.e[j++]];
        if (!p2) {
          runtimeErrorMsg() << "oif global calc: bond broken between particles "
                            << p1->p.identity << " and " << p1->bl.e[j - 1]
                            << " (particles not stored on the same node - "
                               "oif_global_forces1); n "
                            << p1->bl.n << " max " << p1->bl.max;
          return;
        }
        /* fetch particle 3 */
        // if(n_partners>2){
        p3 = local_particles[p1->bl.e[j++]];
        if (!p3) {
          runtimeErrorMsg() << "oif global calc: bond broken between particles "
                            << p1->p.identity << ", " << p1->bl.e[j - 2]
                            << " and " << p1->bl.e[j - 1]
                            << " (particles not stored on the same node - "
                               "oif_global_forces1); n "
                            << p1->bl.n << " max " << p1->bl.max;
          return;
        }
        // remaining neighbors fetched

        // getting unfolded positions of all particles
        // first find out which particle out of p1, p2 (possibly p3, p4) is not
        // a ghost particle. In almost all cases it is p1, however, it might be
        // other one. we call this particle reference particle.
        if (p1->l.ghost != 1) {
          // unfold non-ghost particle using image, because for physical
          // particles, the structure p->l.i is correctly set
          p11 = unfolded_position(p1);
          // other coordinates are obtained from its relative positions to the
          // reference particle
          get_mi_vector(AA, p2->r.p, p11);
          get_mi_vector(BB, p3->r.p, p11);
          for (int i = 0; i < 3; i++) {
            p22[i] = p11[i] + AA[i];
            p33[i] = p11[i] + BB[i];
          }
        } else {
          // in case the first particle is a ghost particle
          if (p2->l.ghost != 1) {
            p22 = unfolded_position(p2);
            get_mi_vector(AA, p1->r.p, p22);
            get_mi_vector(BB, p3->r.p, p22);
            for (int i = 0; i < 3; i++) {
              p11[i] = p22[i] + AA[i];
              p33[i] = p22[i] + BB[i];
            }
          } else {
            // in case the first and the second particle are ghost particles
            if (p3->l.ghost != 1) {
              p33 = unfolded_position(p3);
              get_mi_vector(AA, p1->r.p, p33);
              get_mi_vector(BB, p2->r.p, p33);
              for (int i = 0; i < 3; i++) {
                p11[i] = p33[i] + AA[i];
                p22[i] = p33[i] + BB[i];
              }
            } else {
              printf("Something wrong in oif_global_forces.hpp: All particles "
                     "in a bond are ghost particles, impossible to unfold the "
                     "positions...");
              return;
            }
          }
        }

        // unfolded positions correct
        auto const VOL_A = area_triangle(p11, p22, p33);
        partArea += VOL_A;

        auto const VOL_norm = get_n_triangle(p11, p22, p33);
        auto const VOL_dn = VOL_norm.norm();
        auto const VOL_hz = 1.0 / 3.0 * (p11[2] + p22[2] + p33[2]);
        VOL_partVol += VOL_A * -1 * VOL_norm[2] / VOL_dn * VOL_hz;
      } else {
        j += n_partners;
      }
    }
  }

  part_area_volume[0] = partArea;
  part_area_volume[1] = VOL_partVol;

  MPI_Allreduce(part_area_volume, area_volume, 2, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

void add_oif_global_forces(double const *area_volume,
                           int molType) { // first-fold-then-the-same approach
  double area = area_volume[0];
  double VOL_volume = area_volume[1];

  int test = 0;

  for (auto &p : local_cells.particles()) {
    int j = 0;
    auto p1 = &p;
    while (j < p1->bl.n) {
      /* bond type */
      auto const type_num = p1->bl.e[j++];
      auto iaparams = &bonded_ia_params[type_num];
      auto const type = iaparams->type;
      auto const n_partners = iaparams->num;
      auto const id = p1->p.mol_id;
      if (type == BONDED_IA_OIF_GLOBAL_FORCES &&
          id == molType) { // BONDED_IA_OIF_GLOBAL_FORCES with correct molType
        test++;
        /* fetch particle 2 */
        auto p2 = local_particles[p1->bl.e[j++]];
        if (!p2) {
          runtimeErrorMsg() << "add area: bond broken between particles "
                            << p1->p.identity << " and " << p1->bl.e[j - 1]
                            << " (particles not stored on the same node - "
                               "oif_globalforce2); n "
                            << p1->bl.n << " max " << p1->bl.max;
          return;
        }
        /* fetch particle 3 */
        // if(n_partners>2){
        auto p3 = local_particles[p1->bl.e[j++]];
        if (!p3) {
          runtimeErrorMsg()
              << "add area: bond broken between particles " << p1->p.identity
              << ", " << p1->bl.e[j - 2] << " and " << p1->bl.e[j - 1]
              << " (particles not stored on the same node); n " << p1->bl.n
              << " max " << p1->bl.max;
          return;
        }

        auto const p11 = unfolded_position(*p1);
        auto const p22 = p11 + get_mi_vector(p2->r.p, p11);
        auto const p33 = p11 + get_mi_vector(p3->r.p, p11);

        // unfolded positions correct
        /// starting code from volume force
        auto const VOL_norm = get_n_triangle(p11, p22, p33).normalize();
        auto const VOL_A = area_triangle(p11, p22, p33);
        auto const VOL_vv = (VOL_volume - iaparams->p.oif_global_forces.V0) /
                            iaparams->p.oif_global_forces.V0;

        auto const VOL_force = (1.0 / 3.0) * iaparams->p.oif_global_forces.kv *
                               VOL_vv * VOL_A * VOL_norm;
        p1->f.f += VOL_force;
        p2->f.f += VOL_force;
        p3->f.f += VOL_force;
        ///  ending code from volume force

        auto const h = (1. / 3.) * (p11 + p22 + p33);

        auto const deltaA = (area - iaparams->p.oif_global_forces.A0_g) /
                            iaparams->p.oif_global_forces.A0_g;

        auto const m1 = h - p11;
        auto const m2 = h - p22;
        auto const m3 = h - p33;

        auto const m1_length = m1.norm();
        auto const m2_length = m2.norm();
        auto const m3_length = m3.norm();

        auto const fac = iaparams->p.oif_global_forces.ka_g * VOL_A * deltaA /
                         (m1_length * m1_length + m2_length * m2_length +
                          m3_length * m3_length);

        p1->f.f += fac * m1;
        p2->f.f += fac * m2;
        p3->f.f += fac * m3;
      } else {
        j += n_partners;
      }
    }
  }
}

int max_oif_objects = 0;
