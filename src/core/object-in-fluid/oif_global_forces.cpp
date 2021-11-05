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

#include "BoxGeometry.hpp"
#include "CellStructure.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "grid.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/triangle_functions.hpp>

#include <boost/mpi/collectives.hpp>

#include <functional>

using Utils::area_triangle;
using Utils::get_n_triangle;

Utils::Vector2d calc_oif_global(int molType, CellStructure &cs) {
  // first-fold-then-the-same approach
  double partArea = 0.0;
  // z volume
  double VOL_partVol = 0.;

  cs.bond_loop([&partArea, &VOL_partVol,
                molType](Particle &p1, int bond_id,
                         Utils::Span<Particle *> partners) {
    if (p1.p.mol_id != molType)
      return false;

    if (boost::get<OifGlobalForcesBond>(bonded_ia_params.at(bond_id).get()) !=
        nullptr) {
      // remaining neighbors fetched
      auto const p11 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
      auto const p22 = p11 + box_geo.get_mi_vector(partners[0]->r.p, p11);
      auto const p33 = p11 + box_geo.get_mi_vector(partners[1]->r.p, p11);

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

  auto const area_volume_local = Utils::Vector2d{{partArea, VOL_partVol}};
  return boost::mpi::all_reduce(comm_cart, area_volume_local, std::plus<>());
}

void add_oif_global_forces(Utils::Vector2d const &area_volume, int molType,
                           CellStructure &cs) {
  // first-fold-then-the-same approach
  double area = area_volume[0];
  double VOL_volume = area_volume[1];

  cs.bond_loop([area, VOL_volume, molType](Particle &p1, int bond_id,
                                           Utils::Span<Particle *> partners) {
    if (p1.p.mol_id != molType)
      return false;

    if (auto const *iaparams = boost::get<OifGlobalForcesBond>(
            bonded_ia_params.at(bond_id).get())) {
      auto const p11 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
      auto const p22 = p11 + box_geo.get_mi_vector(partners[0]->r.p, p11);
      auto const p33 = p11 + box_geo.get_mi_vector(partners[1]->r.p, p11);

      // unfolded positions correct
      // starting code from volume force
      auto const VOL_norm = get_n_triangle(p11, p22, p33).normalize();
      auto const VOL_A = area_triangle(p11, p22, p33);
      auto const VOL_vv = (VOL_volume - iaparams->V0) / iaparams->V0;

      auto const VOL_force =
          (1.0 / 3.0) * iaparams->kv * VOL_vv * VOL_A * VOL_norm;

      auto const h = (1. / 3.) * (p11 + p22 + p33);

      auto const deltaA = (area - iaparams->A0_g) / iaparams->A0_g;

      auto const m1 = h - p11;
      auto const m2 = h - p22;
      auto const m3 = h - p33;

      auto const m1_length = m1.norm();
      auto const m2_length = m2.norm();
      auto const m3_length = m3.norm();

      auto const fac = iaparams->ka_g * VOL_A * deltaA /
                       (m1_length * m1_length + m2_length * m2_length +
                        m3_length * m3_length);

      p1.f.f += fac * m1 + VOL_force;
      partners[0]->f.f += fac * m2 + VOL_force;
      partners[1]->f.f += fac * m3 + VOL_force;
    }

    return false;
  });
}

int max_oif_objects = 0;

void mpi_set_max_oif_objects_local(int max_oif_objects) {
  ::max_oif_objects = max_oif_objects;
}

REGISTER_CALLBACK(mpi_set_max_oif_objects_local)

void mpi_set_max_oif_objects(int max_oif_objects) {
  mpi_call_all(mpi_set_max_oif_objects_local, max_oif_objects);
}
