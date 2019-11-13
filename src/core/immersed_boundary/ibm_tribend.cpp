/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "immersed_boundary/ibm_tribend.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>

std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
IBM_Tribend_CalcForce(Particle const &p1, Particle const &p2,
                      Particle const &p3, Particle const &p4,
                      Bonded_ia_parameters const &iaparams) {

  // Get vectors making up the two triangles
  auto const dx1 = get_mi_vector(p1.r.p, p3.r.p, box_geo);
  auto const dx2 = get_mi_vector(p2.r.p, p3.r.p, box_geo);
  auto const dx3 = get_mi_vector(p4.r.p, p3.r.p, box_geo);

  // Get normals on triangle; pointing outwards by definition of indices
  // sequence
  auto n1 = vector_product(dx1, dx2);
  auto n2 = vector_product(dx3, dx1);

  // Get 2*area of triangles out of the magnitude of the resulting normals and
  // make the latter unity
  auto const Ai = n1.norm();
  n1 /= Ai;

  auto const Aj = n2.norm();
  n2 /= Aj;

  // Get the prefactor for the force term
  auto sc = n1 * n2;
  if (sc > 1.0)
    sc = 1.0;

  // Get theta as angle between normals
  auto theta = acos(sc);

  auto const direc = vector_product(n1, n2);
  auto const desc = (dx1 * direc);

  if (desc < 0)
    theta *= -1;

  auto const DTh = theta - iaparams.p.ibm_tribend.theta0;

  auto Pre = iaparams.p.ibm_tribend.kb * DTh;
  // Correct version with linearized sin
  if (theta < 0)
    Pre *= -1;

  auto const v1 = (n2 - sc * n1).normalize();
  auto const v2 = (n1 - sc * n2).normalize();

  // Force on particles
  auto const force1 =
      Pre * (vector_product(get_mi_vector(p2.r.p, p3.r.p, box_geo), v1) / Ai +
             vector_product(get_mi_vector(p3.r.p, p4.r.p, box_geo), v2) / Aj);
  auto const force2 =
      Pre * (vector_product(get_mi_vector(p3.r.p, p1.r.p, box_geo), v1) / Ai);
  auto const force3 =
      Pre * (vector_product(get_mi_vector(p1.r.p, p2.r.p, box_geo), v1) / Ai +
             vector_product(get_mi_vector(p4.r.p, p1.r.p, box_geo), v2) / Aj);
  auto const force4 =
      Pre * (vector_product(get_mi_vector(p1.r.p, p3.r.p, box_geo), v2) / Aj);
  return std::make_tuple(force1, force2, force3, force4);
}

int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2,
                          const int ind3, const int ind4, const double kb,
                          const bool flat) {
  // Create bond
  make_bond_type_exist(bond_type);

  // General parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_TRIBEND;

  // Specific parameters
  //  bonded_ia_params[bond_type].p.ibm_tribend.method = method;

  // Distinguish bending methods
  //  if ( method == TriangleNormals )
  {
    double theta0;

    if (!flat) {
      // Compute theta0
      auto p1 = get_particle_data(ind1);
      auto p2 = get_particle_data(ind2);
      auto p3 = get_particle_data(ind3);
      auto p4 = get_particle_data(ind4);

      // Get vectors of triangles
      auto const dx1 = get_mi_vector(p1.r.p, p3.r.p, box_geo);
      auto const dx2 = get_mi_vector(p2.r.p, p3.r.p, box_geo);
      auto const dx3 = get_mi_vector(p4.r.p, p3.r.p, box_geo);

      // Get normals on triangle; pointing outwards by definition of indices
      // sequence
      auto const n1l = vector_product(dx1, dx2);
      auto const n2l = -vector_product(dx1, dx3);

      auto const n1 = n1l / n1l.norm();
      auto const n2 = n2l / n2l.norm();

      // calculate theta by taking the acos of the scalar n1*n2
      auto sc = n1 * n2;
      if (sc > 1.0)
        sc = 1.0;

      theta0 = acos(sc);

      auto const desc = dx1 * vector_product(n1, n2);
      if (desc < 0)
        theta0 = 2.0 * Utils::pi() - theta0;

    } else
      theta0 = 0; // Flat

    // Krüger always has three partners
    bonded_ia_params[bond_type].num = 3;
    bonded_ia_params[bond_type].p.ibm_tribend.theta0 = theta0;
    // NOTE: This is the bare bending modulus used by the program.
    // If triangle pairs appear only once, the total bending force should get a
    // factor 2. For the numerical model, a factor sqrt(3) should be added, see
    // Gompper&Kroll J. Phys. 1996 and Krüger thesis. This is an approximation,
    // it holds strictly only for a sphere
    bonded_ia_params[bond_type].p.ibm_tribend.kb = kb;
  }

  // Broadcast and return
  mpi_bcast_ia_params(bond_type, -1);
  return ES_OK;
}
