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
#ifndef _HARMONIC_HPP
#define _HARMONIC_HPP
/** \file
 *  Routines to calculate the harmonic bond potential between particle pairs.
 */

#include "config.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>

/** Parameters for harmonic bond Potential */
struct HarmonicBond {
  /** spring constant */
  double k;
  /** equilibrium bond length */
  double r;
  /** cutoff length */
  double r_cut;

  double cutoff() const { return r_cut; }

  static constexpr int num = 1;

  HarmonicBond(double k, double r, double r_cut) {
    this->k = k;
    this->r = r;
    this->r_cut = r_cut;
  }

  boost::optional<Utils::Vector3d> force(Utils::Vector3d const &dx) const;
  boost::optional<double> energy(Utils::Vector3d const &dx) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &k;
    ar &r;
    ar &r_cut;
  }
};

/** Compute the harmonic bond force.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<Utils::Vector3d>
HarmonicBond::force(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if ((r_cut > 0.0) && (dist > r_cut)) {
    return {};
  }

  auto const dr = dist - r;
  auto fac = -k * dr;
  if (dist > ROUND_ERROR_PREC) { /* Regular case */
    fac /= dist;
  } else {
    fac = 0;
  }
  return fac * dx;
}

/** Compute the harmonic bond energy.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<double>
HarmonicBond::energy(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if ((r_cut > 0.0) && (dist > r_cut)) {
    return {};
  }

  return 0.5 * k * Utils::sqr(dist - r);
}

#endif
