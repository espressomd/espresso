/*
 * Copyright (C) 2012-2019 The ESPResSo project
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

#include "oif_global_forces.hpp"

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include <utils/math/triangle_functions.hpp>
using Utils::angle_btw_triangles;
using Utils::area_triangle;
using Utils::get_n_triangle;

#include <utils/constants.hpp>

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

void calc_oif_global(double *area_volume, int molType, CellStructure &cs) {
  // first-fold-then-the-same approach
  double partArea = 0.0;
  // z volume
  double VOL_partVol = 0.;

  double part_area_volume[2]; // added

  for (auto &p : cs.local_particles()) {
    if (p.p.mol_id != molType)
      continue;

    cs.execute_bond_handler(p, [&partArea, &VOL_partVol](
                                   Particle &p1, int bond_id,
                                   Utils::Span<Particle *> partners) {
      auto const &iaparams = bonded_ia_params[bond_id];

      if (iaparams.type == BONDED_IA_OIF_GLOBAL_FORCES) {
        // remaining neighbors fetched
        auto const p11 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
        auto const p22 = p11 + get_mi_vector(partners[0]->r.p, p11, box_geo);
        auto const p33 = p11 + get_mi_vector(partners[1]->r.p, p11, box_geo);

        // unfolded positions correct
        auto const VOL_A = area_triangle(p11, p22, p33);
        partArea += VOL_A;

        auto const VOL_norm = get_n_triangle(p11, p22, p33);
        auto const VOL_dn = VOL_norm.norm();
        auto const VOL_hz = 1.0 / 3.0 * (p11[2] + p22[2] + p33[2]);
        VOL_partVol += VOL_A * -1 * VOL_norm[2] / VOL_dn * VOL_hz;
      }

      return false;
    });
  }

  part_area_volume[0] = partArea;
  part_area_volume[1] = VOL_partVol;

  MPI_Allreduce(part_area_volume, area_volume, 2, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

void add_oif_global_forces(double const *area_volume, int molType,
                           CellStructure &cs) {
  // first-fold-then-the-same approach
  double area = area_volume[0];
  double VOL_volume = area_volume[1];

  for (auto &p : cs.local_particles()) {
    if (p.p.mol_id != molType)
      continue;

    cs.execute_bond_handler(p, [area,
                                VOL_volume](Particle &p1, int bond_id,
                                            Utils::Span<Particle *> partners) {
      auto const &iaparams = bonded_ia_params[bond_id];

      if (iaparams.type == BONDED_IA_OIF_GLOBAL_FORCES) {
        auto const p11 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
        auto const p22 = p11 + get_mi_vector(partners[0]->r.p, p11, box_geo);
        auto const p33 = p11 + get_mi_vector(partners[1]->r.p, p11, box_geo);

        // unfolded positions correct
        // starting code from volume force
        auto const VOL_norm = get_n_triangle(p11, p22, p33).normalize();
        auto const VOL_A = area_triangle(p11, p22, p33);
        auto const VOL_vv = (VOL_volume - iaparams.p.oif_global_forces.V0) /
                            iaparams.p.oif_global_forces.V0;

        auto const VOL_force = (1.0 / 3.0) * iaparams.p.oif_global_forces.kv *
                               VOL_vv * VOL_A * VOL_norm;

        auto const h = (1. / 3.) * (p11 + p22 + p33);

        auto const deltaA = (area - iaparams.p.oif_global_forces.A0_g) /
                            iaparams.p.oif_global_forces.A0_g;

        auto const m1 = h - p11;
        auto const m2 = h - p22;
        auto const m3 = h - p33;

        auto const m1_length = m1.norm();
        auto const m2_length = m2.norm();
        auto const m3_length = m3.norm();

        auto const fac = iaparams.p.oif_global_forces.ka_g * VOL_A * deltaA /
                         (m1_length * m1_length + m2_length * m2_length +
                          m3_length * m3_length);

        p1.f.f += fac * m1 + VOL_force;
        partners[0]->f.f += fac * m2 + VOL_force;
        partners[1]->f.f += fac * m3 + VOL_force;
      }

      return false;
    });
  }
}

int max_oif_objects = 0;
