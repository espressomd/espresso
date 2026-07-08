/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef GAUSSIAN_ANISO_H
#define GAUSSIAN_ANISO_H

/** \file
 *  Routines to calculate an anisotropic 3D Gaussian potential between
 *  particle pairs with independent widths sigma_x, sigma_y, sigma_z.
 *
 *  Potential:
 *    U(dx, dy, dz) = eps * exp( -0.5 * [ (dx/sig_x)^2
 *                                        + (dy/sig_y)^2
 *                                        + (dz/sig_z)^2 ] )
 *
 *  If sig_x == sig_y == sig_z, this reduces to the ordinary isotropic
 *  3D Gaussian (up to using a radial cutoff on |r|, same semantics
 *  as the built-in GAUSSIAN potential).
 *
 *  Implementation in \ref gaussian_aniso.cpp.
 */

#include "config/config.hpp"

#include "nonbonded_interaction_data.hpp"

#include <utils/math/sqr.hpp>
#include <utils/Vector.hpp>
#include <cmath>

#ifdef ESPRESSO_GAUSSIAN_ANISO

/** Calculate anisotropic Gaussian energy.
 *
 *  Parameters from ia_params.gaussian_aniso:
 *    eps    : amplitude
 *    sig_x  : width along x
 *    sig_y  : width along y
 *    sig_z  : width along z
 *    cut    : radial cutoff on |r| (same semantics as GAUSSIAN)
 *
 *  Input:
 *    dx, dy, dz : displacement components (r_j - r_i)
 *
 *  Returns:
 *    U_ij
 */
inline double gaussian_aniso_pair_energy(IA_parameters const &ia_params,
                                         Utils::Vector3d const &d)
{
    double dx = d[0];
    double dy = d[1];
    double dz = d[2];

    auto const &p = ia_params.gaussian_aniso;

    const double r2 = Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz);

    if (r2 >= Utils::sqr(p.cut)) {
        return 0.0;
    }

    const double sig_x2 = Utils::sqr(p.sig_x);
    const double sig_y2 = Utils::sqr(p.sig_y);
    const double sig_z2 = Utils::sqr(p.sig_z);

    const double A = 0.5 * (Utils::sqr(dx) / sig_x2
                          + Utils::sqr(dy) / sig_y2
                          + Utils::sqr(dz) / sig_z2);

    return p.eps * std::exp(-A);
}


/** Calculate anisotropic Gaussian force.
 *
 *  This is the analogue of gaussian_pair_force_factor, but because
 *  the potential is anisotropic, the force is not simply parallel to r,
 *  so we return the components directly.
 *
 *  Input:
 *    dx, dy, dz : displacement components (r_j - r_i)
 *
 *  Output:
 *    Force vector on particle i due to particle j
 */
inline Utils::Vector3d gaussian_aniso_pair_force(
    IA_parameters const &ia_params, Utils::Vector3d const &d) {
  auto const &p = ia_params.gaussian_aniso;

  const double dx = d[0];
  const double dy = d[1];
  const double dz = d[2];

  const double r2 = Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz);

  if (r2 >= Utils::sqr(p.cut)) {
    return {};
  }

  const double sig_x2 = Utils::sqr(p.sig_x);
  const double sig_y2 = Utils::sqr(p.sig_y);
  const double sig_z2 = Utils::sqr(p.sig_z);

  const double A = 0.5 * (Utils::sqr(dx) / sig_x2
                        + Utils::sqr(dy) / sig_y2
                        + Utils::sqr(dz) / sig_z2);
                      
  const double U = p.eps * std::exp(-A);

  // F = -∇U
  // With d = r_j - r_i, check sign against ESPResSo's force convention.
  return {U * dx / sig_x2, U * dy / sig_y2, U * dz / sig_z2};
}

#endif /* ifdef ESPRESSO_GAUSSIAN_ANISO */
#endif /* GAUSSIAN_ANISO_H */

