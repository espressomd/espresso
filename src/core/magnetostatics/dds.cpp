/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "config/config.hpp"

#ifdef DIPOLES

#include "magnetostatics/dds.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

/**
 * Calculate dipolar energy and optionally force between two particles.
 * @param[in,out] p1          First particle
 * @param[in]     dip1        Cached dipole moment of the first particle
 * @param[in,out] p2          Second particle
 * @param[in]     force_flag  If true, update the particle forces and torques
 */
double DipolarDirectSum::calc_dipole_dipole_ia(Particle &p1,
                                               Utils::Vector3d const &dip1,
                                               Particle &p2,
                                               bool force_flag) const {

  // Cache dipole moment
  auto const dip2 = p2.calc_dip();

  // Distance between particles
  auto const dr = box_geo.get_mi_vector(p1.pos(), p2.pos());

  // Powers of distance
  auto const r2 = dr.norm2();
  auto const r = sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;
  auto const r7 = r5 * r2;

  // Dot products
  auto const pe1 = dip1 * dip2;
  auto const pe2 = dip1 * dr;
  auto const pe3 = dip2 * dr;
  auto const pe4 = 3.0 / r5;

  // Energy
  auto const energy = prefactor * (pe1 / r3 - pe4 * pe2 * pe3);

  // Forces, if requested
  if (force_flag) {
    auto const a = pe4 * pe1;
    auto const b = -15.0 * pe2 * pe3 / r7;
    auto const ab = a + b;
    auto const cc = pe4 * pe3;
    auto const dd = pe4 * pe2;

    // Forces
    auto const ff = ab * dr + cc * dip1 + dd * dip2;
    p1.force() += prefactor * ff;
    p2.force() -= prefactor * ff;

    // Torques
    auto const aa = vector_product(dip1, dip2);
    auto const b1 = vector_product(dip1, dr);
    auto const b2 = vector_product(dip2, dr);
    p1.torque() += prefactor * (-aa / r3 + b1 * cc);
    p2.torque() += prefactor * (aa / r3 + b2 * dd);
  }

  return energy;
}

double DipolarDirectSum::kernel(bool force_flag, bool energy_flag,
                                ParticleRange const &particles) const {

  assert(n_nodes == 1);
  assert(force_flag || energy_flag);

  double energy = 0.0;
  // Iterate over all particles
  for (auto it = particles.begin(), end = particles.end(); it != end; ++it) {
    // If the particle has no dipole moment, ignore it
    if (it->dipm() == 0.0)
      continue;

    auto const dip1 = it->calc_dip();
    auto jt = it;
    /* Skip diagonal */
    ++jt;
    for (; jt != end; ++jt) {
      // If the particle has no dipole moment, ignore it
      if (jt->dipm() == 0.0)
        continue;
      // Calculate energy and/or force between the particles
      energy += calc_dipole_dipole_ia(*it, dip1, *jt, force_flag);
    }
  }

  return energy;
}

DipolarDirectSum::DipolarDirectSum(double prefactor) : prefactor{prefactor} {
  if (n_nodes > 1) {
    throw std::runtime_error(
        "MPI parallelization not supported by DipolarDirectSumCpu.");
  }
  if (prefactor <= 0.) {
    throw std::domain_error("Parameter 'prefactor' must be > 0");
  }
}

#endif // DIPOLES
