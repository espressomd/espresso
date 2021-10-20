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

#ifndef IBM_TRIEL_H
#define IBM_TRIEL_H

#include "Particle.hpp"
#include "config.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <tuple>

enum class tElasticLaw { NeoHookean, Skalak };

/** Parameters for IBM elastic triangle (triel) */
struct IBMTriel {
  // These values encode the reference state
  double l0;
  double lp0;
  double sinPhi0;
  double cosPhi0;
  double area0;

  // These values are cache values to speed up computation
  double a1;
  double a2;
  double b1;
  double b2;

  // These are interaction parameters
  // k1 is used for Neo-Hookean
  // k1 and k2 are used Skalak
  double maxDist;
  tElasticLaw elasticLaw;
  double k1;
  double k2;

  double cutoff() const { return maxDist; }

  static constexpr int num = 2;

  /** Set the IBM Triel parameters.
   *  Also calculate and store the reference state.
   */
  IBMTriel(int ind1, int ind2, int ind3, double maxDist, tElasticLaw elasticLaw,
           double k1, double k2);

  /** Calculate the forces.
   *  The equations can be found in Appendix C of @cite kruger12a.
   *  @return the forces on @p p1, @p p2, @p p3
   */
  boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
  calc_forces(Particle const &p1, Particle const &p2, Particle const &p3) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &l0;
    ar &lp0;
    ar &sinPhi0;
    ar &cosPhi0;
    ar &area0;
    ar &a1;
    ar &a2;
    ar &b1;
    ar &b2;
    ar &maxDist;
    ar &elasticLaw;
    ar &k1;
    ar &k2;
  }
};

#endif
