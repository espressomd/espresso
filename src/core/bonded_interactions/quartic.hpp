/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef _QUARTIC_HPP
#define _QUARTIC_HPP
/** \file
 *  Routines to calculate the quartic potential between particle pairs.
 */

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>

/** Parameters for quartic bond Potential */
struct QuarticBond {
  double k0, k1;
  double r;
  double r_cut;

  double cutoff() const { return r_cut; }

  static constexpr int num = 1;

  QuarticBond(double k0, double k1, double r, double r_cut) {
    this->k0 = k0;
    this->k1 = k1;
    this->r = r;
    this->r_cut = r_cut;
  }

  boost::optional<Utils::Vector3d> force(Utils::Vector3d const &dx) const;
  boost::optional<double> energy(Utils::Vector3d const &dx) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &k0;
    ar &k1;
    ar &r;
    ar &r_cut;
  }
};

/** Compute the quartic bond force.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<Utils::Vector3d>
QuarticBond::force(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if ((r_cut > 0.0) && (dist > r_cut)) {
    return {};
  }

  auto const dr = dist - r;
  auto const fac = (k0 * dr + k1 * Utils::int_pow<3>(dr)) / dist;
  return -fac * dx;
}

/** Compute the quartic bond energy.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<double>
QuarticBond::energy(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if ((r_cut > 0.0) && (dist > r_cut)) {
    return {};
  }

  double dr2 = Utils::sqr(dist - r);

  return 0.5 * k0 * dr2 + 0.25 * k1 * Utils::sqr(dr2);
}

#endif
