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
 *  Routines to calculate the short-range part of the bonded Coulomb potential
 *  between particle pairs. Can be used to subtract certain intramolecular
 *  interactions in combination with Thole damping.
 */

#include "config/config.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cmath>
#include <functional>

/** Parameters for Coulomb bond short-range Potential */
struct BondedCoulombSR {
  /** charge factor */
  double q1q2;

  double cutoff() const { return 0.; }

  static constexpr int num = 1;

  BondedCoulombSR(double q1q2) { this->q1q2 = q1q2; }

  boost::optional<Utils::Vector3d>
  force(Utils::Vector3d const &dx,
        std::function<Utils::Vector3d(double, Utils::Vector3d const &,
                                      double)> const &kernel) const;
  boost::optional<double>
  energy(Particle const &p1, Particle const &p2, Utils::Vector3d const &dx,
         std::function<double(Particle const &, Particle const &, double,
                              Utils::Vector3d const &, double)> const &kernel)
      const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & q1q2;
  }
};

/** Compute the short-range bonded Coulomb pair force.
 *  @param[in]  dx        Distance between the particles.
 *  @param[in]  kernel    Coulomb force kernel.
 */
inline boost::optional<Utils::Vector3d> BondedCoulombSR::force(
    Utils::Vector3d const &dx,
    std::function<Utils::Vector3d(double, Utils::Vector3d const &,
                                  double)> const &kernel) const {
#ifdef ELECTROSTATICS
  return kernel(q1q2, dx, dx.norm());
#else
  return Utils::Vector3d{};
#endif
}

/** Compute the short-range bonded Coulomb pair energy.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  dx        Distance between the particles.
 *  @param[in]  kernel    Coulomb energy kernel.
 */
inline boost::optional<double> BondedCoulombSR::energy(
    Particle const &p1, Particle const &p2, Utils::Vector3d const &dx,
    std::function<double(Particle const &, Particle const &, double,
                         Utils::Vector3d const &, double)> const &kernel)
    const {
#ifdef ELECTROSTATICS
  return kernel(p1, p2, q1q2, dx, dx.norm());
#else
  return 0.;
#endif
}
