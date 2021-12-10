
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

#ifndef THERMALIZED_DIST_H
#define THERMALIZED_DIST_H
/** \file
 *  Routines to thermalize the center of mass and distance of a particle pair.
 *
 *  Implementation in \ref thermalized_bond.cpp.
 */

/** number of thermalized bonds */
extern int n_thermalized_bonds;

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <tuple>

/** Parameters for Thermalized bond */
struct ThermalizedBond {
  double temp_com;
  double gamma_com;
  double temp_distance;
  double gamma_distance;
  double r_cut;
  double pref1_com;
  double pref2_com;
  double pref1_dist;
  double pref2_dist;

  double cutoff() const { return r_cut; }

  static constexpr int num = 1;

  ThermalizedBond(double temp_com, double gamma_com, double temp_distance,
                  double gamma_distance, double r_cut);
  boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d>>
  forces(Particle const &p1, Particle const &p2,
         Utils::Vector3d const &dx) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &temp_com;
    ar &gamma_com;
    ar &temp_distance;
    ar &gamma_distance;
    ar &r_cut;
    ar &pref1_com;
    ar &pref2_com;
    ar &pref1_dist;
    ar &pref2_dist;
  }
};

#endif
