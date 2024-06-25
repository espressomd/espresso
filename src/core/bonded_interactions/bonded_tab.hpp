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
 *  Routines to calculate the energy and/or force for particle bonds, angles
 *  and dihedrals via interpolation of lookup tables.
 */

#include "config/config.hpp"

#include "TabulatedPotential.hpp"
#include "angle_common.hpp"
#include "bonded_interactions/dihedral.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/serialization/shared_ptr.hpp>

#include <cassert>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <tuple>
#include <vector>

/** Base class for n-body tabulated potential (n=2,3,4). */
struct TabulatedBond {
  std::shared_ptr<TabulatedPotential> pot;

  /** Set the parameters of a bonded tabulated potential.
   *  ia_params and force/energy tables are communicated to each node.
   *
   *  @param min          @copybrief TabulatedPotential::minval
   *  @param max          @copybrief TabulatedPotential::maxval
   *  @param energy       @copybrief TabulatedPotential::energy_tab
   *  @param force        @copybrief TabulatedPotential::force_tab
   */
  TabulatedBond(double min, double max, std::vector<double> const &energy,
                std::vector<double> const &force) {
    pot = std::make_shared<TabulatedPotential>(min, max, force, energy);
  }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & pot;
  }
};

/** Parameters for 2-body tabulated potential. */
struct TabulatedDistanceBond : public TabulatedBond {
  double cutoff() const { return assert(pot), pot->cutoff(); }

  static constexpr int num = 1;

  TabulatedDistanceBond(double min, double max,
                        std::vector<double> const &energy,
                        std::vector<double> const &force)
      : TabulatedBond(min, max, energy, force) {
    this->pot->minval = min;
    this->pot->maxval = max;
  }

  std::optional<Utils::Vector3d> force(Utils::Vector3d const &dx) const;
  std::optional<double> energy(Utils::Vector3d const &dx) const;
};

/** Parameters for 3-body tabulated potential. */
struct TabulatedAngleBond : public TabulatedBond {
  double cutoff() const { return 0.; }

  static constexpr int num = 2;

  TabulatedAngleBond(double min, double max, std::vector<double> const &energy,
                     std::vector<double> const &force)
      : TabulatedBond(min, max, energy, force) {
    this->pot->minval = 0.;
    this->pot->maxval = std::numbers::pi + ROUND_ERROR_PREC;
  }

  std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
  forces(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2) const;
  double energy(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2) const;
};

/** Parameters for 4-body tabulated potential. */
struct TabulatedDihedralBond : public TabulatedBond {
  double cutoff() const { return 0.; }

  static constexpr int num = 3;

  TabulatedDihedralBond(double min, double max,
                        std::vector<double> const &energy,
                        std::vector<double> const &force)
      : TabulatedBond(min, max, energy, force) {
    this->pot->minval = 0.;
    this->pot->maxval = 2. * std::numbers::pi + ROUND_ERROR_PREC;
  }

  std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d,
                           Utils::Vector3d>>
  forces(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
         Utils::Vector3d const &v34) const;
  std::optional<double> energy(Utils::Vector3d const &v12,
                               Utils::Vector3d const &v23,
                               Utils::Vector3d const &v34) const;
};

/** Compute a tabulated bond length force.
 *
 *  The force acts in the direction of the connecting vector between the
 *  particles. For distances smaller than the tabulated range it uses a linear
 *  extrapolation based on the first two tabulated force values.
 *
 *  @param[in]  dx        Distance between the particles.
 */
inline std::optional<Utils::Vector3d>
TabulatedDistanceBond::force(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if (dist < pot->cutoff()) {
    auto const fac = pot->force(dist) / dist;
    return fac * dx;
  }
  return {};
}

/** Compute a tabulated bond length energy.
 *
 *  For distances smaller than the tabulated range it uses a quadratic
 *  extrapolation based on the first two tabulated force values and the first
 *  tabulated energy value.
 *
 *  @param[in]  dx        Distance between the particles.
 */
inline std::optional<double>
TabulatedDistanceBond::energy(Utils::Vector3d const &dx) const {
  auto const dist = dx.norm();

  if (dist < pot->cutoff()) {
    return pot->energy(dist);
  }
  return {};
}

/** Compute the three-body angle interaction force.
 *  @param[in]  vec1  Vector from central particle to left particle.
 *  @param[in]  vec2  Vector from central particle to right particle.
 *  @return Forces on the second, first and third particles, in that order.
 */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
TabulatedAngleBond::forces(Utils::Vector3d const &vec1,
                           Utils::Vector3d const &vec2) const {

  auto forceFactor = [this](double const cos_phi) {
    auto const sin_phi = sqrt(1 - Utils::sqr(cos_phi));
#ifdef TABANGLEMINUS
    auto const phi = acos(-cos_phi);
#else
    auto const phi = acos(cos_phi);
#endif
    auto const tab_pot = pot;
    auto const gradient = tab_pot->force(phi);
    return -gradient / sin_phi;
  };

  return angle_generic_force(vec1, vec2, forceFactor, true);
}

/** Compute the three-body angle interaction energy.
 *  It is assumed that the potential is tabulated
 *  for all angles between 0 and Pi.
 *
 *  @param[in]  vec1  Vector from central particle to left particle.
 *  @param[in]  vec2  Vector from central particle to right particle.
 */
inline double TabulatedAngleBond::energy(Utils::Vector3d const &vec1,
                                         Utils::Vector3d const &vec2) const {
  auto const cos_phi = calc_cosine(vec1, vec2, true);
  /* calculate phi */
#ifdef TABANGLEMINUS
  auto const phi = acos(-cos_phi);
#else
  auto const phi = acos(cos_phi);
#endif
  return pot->energy(phi);
}

/** Compute the four-body dihedral interaction force.
 *  The forces have a singularity at @f$ \phi = 0 @f$ and @f$ \phi = \pi @f$
 *  (see @cite swope92a page 592).
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 *  @return the forces on @p p2, @p p1, @p p3
 */
inline std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                Utils::Vector3d, Utils::Vector3d>>
TabulatedDihedralBond::forces(Utils::Vector3d const &v12,
                              Utils::Vector3d const &v23,
                              Utils::Vector3d const &v34) const {
  /* vectors for dihedral angle calculation */
  Utils::Vector3d v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cos_phi;

  /* dihedral angle */
  auto const angle_is_undefined = calc_dihedral_angle(
      v12, v23, v34, v12Xv23, l_v12Xv23, v23Xv34, l_v23Xv34, cos_phi, phi);
  /* dihedral angle not defined - force zero */
  if (angle_is_undefined) {
    return {};
  }

  auto const f1 = (v23Xv34 - cos_phi * v12Xv23) / l_v12Xv23;
  auto const f4 = (v12Xv23 - cos_phi * v23Xv34) / l_v23Xv34;

  auto const v23Xf1 = vector_product(v23, f1);
  auto const v23Xf4 = vector_product(v23, f4);
  auto const v34Xf4 = vector_product(v34, f4);
  auto const v12Xf1 = vector_product(v12, f1);

  /* table lookup */
  auto const fac = pot->force(phi);

  /* store dihedral forces */
  auto const force1 = fac * v23Xf1;
  auto const force2 = fac * (v34Xf4 - v12Xf1 - v23Xf1);
  auto const force3 = fac * (v12Xf1 - v23Xf4 - v34Xf4);

  return std::make_tuple(force2, force1, force3, -(force2 + force1 + force3));
}

/** Compute the four-body dihedral interaction energy.
 *  The energy doesn't have any singularity if the angle phi is well-defined.
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 */
inline std::optional<double>
TabulatedDihedralBond::energy(Utils::Vector3d const &v12,
                              Utils::Vector3d const &v23,
                              Utils::Vector3d const &v34) const {
  /* vectors for dihedral calculations. */
  Utils::Vector3d v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cos_phi;
  auto const angle_is_undefined = calc_dihedral_angle(
      v12, v23, v34, v12Xv23, l_v12Xv23, v23Xv34, l_v23Xv34, cos_phi, phi);
  /* dihedral angle not defined - energy zero */
  if (angle_is_undefined) {
    return {};
  }

  return pot->energy(phi);
}
