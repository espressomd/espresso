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

#ifndef IBM_TRIBEND_H
#define IBM_TRIBEND_H

#include "config.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <tuple>

/** Parameters for IBM tribend */
struct IBMTribend {
  /** Interaction data */
  double kb;

  /** Reference angle */
  double theta0;

  double cutoff() const { return 0.; }

  // Krüger always has three partners
  static constexpr int num = 3;

  /** Set the IBM Tribend parameters.
   *  Also calculate and store the reference state.
   *  See details in @cite gompper96a and @cite kruger12a.
   */
  IBMTribend(int ind1, int ind2, int ind3, int ind4, double kb, bool flat);

  /** Calculate the forces
   *  The equations can be found in Appendix C of @cite kruger12a.
   *  @return forces on @p p1, @p p2, @p p3, @p p4
   */
  std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
  calc_forces(Particle const &p1, Particle const &p2, Particle const &p3,
              Particle const &p4) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &kb;
    ar &theta0;
  }
};

#endif
