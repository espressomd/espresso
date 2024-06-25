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

#pragma once

/** \file
 *  Routines to calculate the FENE potential between particle pairs.
 */

#include "config/config.hpp"
#include "errorhandling.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <optional>

/** Parameters for FENE bond Potential. */
struct FeneBond {
  /** spring constant */
  double k;
  /** maximal bond stretching */
  double drmax;
  /** equilibrium bond length */
  double r0;
  /** square of @p drmax (internal parameter) */
  double drmax2;
  /** inverse square of @p drmax (internal parameter) */
  double drmax2i;

  double cutoff() const { return r0 + drmax; }

  static constexpr int num = 1;

  FeneBond(double k, double drmax, double r0) {
    this->k = k;
    this->drmax = drmax;
    this->r0 = r0;

    this->drmax2 = Utils::sqr(this->drmax);
    this->drmax2i = 1. / this->drmax2;
  }

  std::optional<Utils::Vector3d> force(Utils::Vector3d const &dx) const;
  std::optional<double> energy(Utils::Vector3d const &dx) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & k;
    ar & drmax;
    ar & r0;
    ar & drmax2;
    ar & drmax2i;
  }
};

/** Compute the FENE bond force.
 *  @param[in]  dx        Distance between the particles.
 */
inline std::optional<Utils::Vector3d>
FeneBond::force(Utils::Vector3d const &dx) const {
  auto const len = dx.norm();
  auto const dr = len - r0;

  if (dr >= drmax) {
    return {};
  }

  auto fac = -k * dr / (1.0 - dr * dr * drmax2i);
  if (len > ROUND_ERROR_PREC) {
    fac /= len;
  } else {
    if (r0 > 0.) {
      runtimeErrorMsg() << "FENE bond: Particles have zero distance. "
                           "This is most likely an error in the system setup.";
    }
  }

  return fac * dx;
}

/** Compute the FENE bond energy.
 *  @param[in]  dx        Distance between the particles.
 */
inline std::optional<double> FeneBond::energy(Utils::Vector3d const &dx) const {
  /* compute bond stretching (r-r0) */
  double const dr = dx.norm() - r0;

  /* check bond stretching */
  if (dr >= drmax) {
    return {};
  }

  return -0.5 * k * drmax2 * log(1.0 - dr * dr * drmax2i);
}
