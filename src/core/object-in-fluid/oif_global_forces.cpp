/*
 * Copyright (C) 2012-2022 The ESPResSo project
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
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "oif_global_forces_params.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/math/triangle_functions.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/serialization/utility.hpp>

#include <cassert>
#include <cmath>
#include <functional>
#include <span>

/** Calculate the mesh volume and area. */
static auto calc_oif_mesh(int molType, BoxGeometry const &box_geo,
                          CellStructure &cs,
                          BondedInteractionsMap const &bonded_ias) {

  double area = 0.;
  double volume = 0.;

  cs.bond_loop([&area, &volume, &box_geo, &bonded_ias, molType](
                   Particle &p1, int bond_id, std::span<Particle *> partners) {
    if (p1.mol_id() != molType)
      return false;

    if (boost::get<OifGlobalForcesBond>(bonded_ias.at(bond_id).get())) {
      auto const p11 = box_geo.unfolded_position(p1.pos(), p1.image_box());
      auto const p22 = p11 + box_geo.get_mi_vector(partners[0]->pos(), p11);
      auto const p33 = p11 + box_geo.get_mi_vector(partners[1]->pos(), p11);

      auto const VOL_A = Utils::area_triangle(p11, p22, p33);
      area += VOL_A;

      auto const VOL_norm = Utils::get_n_triangle(p11, p22, p33);
      auto const VOL_dn = VOL_norm.norm();
      auto const VOL_hz = (1.0 / 3.0) * (p11[2] + p22[2] + p33[2]);
      volume -= VOL_A * VOL_norm[2] / VOL_dn * VOL_hz;
    }

    return false;
  });

  return Utils::Vector2d{{area, volume}};
}

/** Distribute the OIF global forces to all particles in the mesh. */
static void add_oif_global_forces(double area, double volume, int molType,
                                  BoxGeometry const &box_geo, CellStructure &cs,
                                  BondedInteractionsMap const &bonded_ias) {

  cs.bond_loop([&box_geo, &bonded_ias, area, volume, molType](
                   Particle &p1, int bond_id, std::span<Particle *> partners) {
    if (p1.mol_id() != molType)
      return false;

    auto const *bond_ptr = bonded_ias.at(bond_id).get();
    if (auto const *bond = boost::get<OifGlobalForcesBond>(bond_ptr)) {
      auto const p11 = box_geo.unfolded_position(p1.pos(), p1.image_box());
      auto const p22 = p11 + box_geo.get_mi_vector(partners[0]->pos(), p11);
      auto const p33 = p11 + box_geo.get_mi_vector(partners[1]->pos(), p11);

      // unfolded positions correct
      // starting code from volume force
      auto const VOL_norm = Utils::get_n_triangle(p11, p22, p33).normalize();
      auto const VOL_A = Utils::area_triangle(p11, p22, p33);
      auto const VOL_vv = (volume - bond->V0) / bond->V0;

      auto const VOL_force = (1. / 3.) * bond->kv * VOL_vv * VOL_A * VOL_norm;

      auto const h = (1. / 3.) * (p11 + p22 + p33);

      auto const deltaA = (area - bond->A0_g) / bond->A0_g;

      auto const m1 = h - p11;
      auto const m2 = h - p22;
      auto const m3 = h - p33;

      auto const m1_length = m1.norm();
      auto const m2_length = m2.norm();
      auto const m3_length = m3.norm();

      auto const fac = bond->ka_g * VOL_A * deltaA /
                       (m1_length * m1_length + m2_length * m2_length +
                        m3_length * m3_length);

      p1.force() += fac * m1 + VOL_force;
      partners[0]->force() += fac * m2 + VOL_force;
      partners[1]->force() += fac * m3 + VOL_force;
    }

    return false;
  });
}

void OifGlobal::run_force_kernel() {
  auto &system = get_system();
  auto &box_geo = *system.box_geo;
  auto &bonded_ias = *system.bonded_ias;
  auto &cell_structure = *system.cell_structure;
  for (int i = 0; i < max_oif_objects; ++i) {
    // There are two global quantities that need to be evaluated:
    // object's surface and object's volume.
    auto const local = calc_oif_mesh(i, box_geo, cell_structure, bonded_ias);
    auto const global = boost::mpi::all_reduce(comm_cart, local, std::plus());
    auto const area = std::abs(global[0]);
    auto const volume = std::abs(global[1]);
    if (area < 1e-100 and volume < 1e-100) {
      break;
    }
    add_oif_global_forces(area, volume, i, box_geo, cell_structure, bonded_ias);
  }
}
