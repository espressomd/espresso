/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "grid.hpp"
#include "ibm_common.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <cmath>
#include <tuple>

std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
IBMTribend::calc_forces(Particle const &p1, Particle const &p2,
                        Particle const &p3, Particle const &p4) const {

  // Get vectors making up the two triangles
  auto const dx1 = box_geo.get_mi_vector(p1.r.p, p3.r.p);
  auto const dx2 = box_geo.get_mi_vector(p2.r.p, p3.r.p);
  auto const dx3 = box_geo.get_mi_vector(p4.r.p, p3.r.p);

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
  auto const sc = std::min(1.0, n1 * n2);

  // Get theta as angle between normals
  auto theta = acos(sc);

  auto const direc = vector_product(n1, n2);
  auto const desc = (dx1 * direc);

  if (desc < 0)
    theta *= -1;

  auto const DTh = theta - theta0;

  auto Pre = kb * DTh;
  // Correct version with linearized sin
  if (theta < 0)
    Pre *= -1;

  auto const v1 = (n2 - sc * n1).normalize();
  auto const v2 = (n1 - sc * n2).normalize();

  // Force on particles: eq. (C.28-C.31)
  auto const force1 =
      Pre * (vector_product(box_geo.get_mi_vector(p2.r.p, p3.r.p), v1) / Ai +
             vector_product(box_geo.get_mi_vector(p3.r.p, p4.r.p), v2) / Aj);
  auto const force2 =
      Pre * (vector_product(box_geo.get_mi_vector(p3.r.p, p1.r.p), v1) / Ai);
  auto const force3 =
      Pre * (vector_product(box_geo.get_mi_vector(p1.r.p, p2.r.p), v1) / Ai +
             vector_product(box_geo.get_mi_vector(p4.r.p, p1.r.p), v2) / Aj);
  auto const force4 =
      Pre * (vector_product(box_geo.get_mi_vector(p1.r.p, p3.r.p), v2) / Aj);
  return std::make_tuple(force1, force2, force3, force4);
}

IBMTribend::IBMTribend(const int ind1, const int ind2, const int ind3,
                       const int ind4, const double kb, const bool flat) {

  // Compute theta0
  if (flat) {
    theta0 = 0;
  } else {
    // Get particles
    auto const pos1 = get_ibm_particle_position(ind1);
    auto const pos2 = get_ibm_particle_position(ind2);
    auto const pos3 = get_ibm_particle_position(ind3);
    auto const pos4 = get_ibm_particle_position(ind4);

    // Get vectors of triangles
    auto const dx1 = box_geo.get_mi_vector(pos1, pos3);
    auto const dx2 = box_geo.get_mi_vector(pos2, pos3);
    auto const dx3 = box_geo.get_mi_vector(pos4, pos3);

    // Get normals on triangle; pointing outwards by definition of indices
    // sequence
    auto const n1l = vector_product(dx1, dx2);
    auto const n2l = -vector_product(dx1, dx3);

    auto const n1 = n1l / n1l.norm();
    auto const n2 = n2l / n2l.norm();

    // calculate theta0 by taking the acos of the scalar n1*n2
    auto const sc = std::min(1.0, n1 * n2);

    theta0 = acos(sc);

    auto const desc = dx1 * vector_product(n1, n2);
    if (desc < 0)
      theta0 = 2.0 * Utils::pi() - theta0;
  }

  // NOTE: This is the bare bending modulus used by the program.
  // If triangle pairs appear only once, the total bending force should get a
  // factor 2. For the numerical model, a factor sqrt(3) should be added, see
  // @cite gompper96a and @cite kruger12a. This is an approximation,
  // it holds strictly only for a sphere
  this->kb = kb;
}
