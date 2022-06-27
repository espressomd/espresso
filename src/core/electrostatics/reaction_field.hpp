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

/**
 * @file
 * Calculate the Reaction Field energy and force
 * for a particle pair @cite neumann85b, @cite tironi95a.
 */

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_REACTION_FIELD_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_REACTION_FIELD_HPP

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/actor.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <stdexcept>

/** @brief Reaction Field parameters. */
struct ReactionField : public Coulomb::Actor<ReactionField> {
  /** @brief Ionic strength. */
  double kappa;
  /** @brief Continuum dielectric constant inside the cavity. */
  double epsilon1;
  /** @brief Continuum dielectric constant outside the cavity. */
  double epsilon2;
  /** @brief Interaction cutoff. */
  double r_cut;
  /** @brief Interaction prefactor. Corresponds to the quantity
   *  @f$ 1 + B_1 @f$ from eq. 22 in @cite tironi95a.
   */
  double B;
  ReactionField(double prefactor, double kappa, double epsilon1,
                double epsilon2, double r_cut) {
    if (kappa < 0.0) {
      throw std::domain_error("Parameter 'kappa' must be >= 0");
    }
    if (epsilon1 < 0.0) {
      throw std::domain_error("Parameter 'epsilon1' must be >= 0");
    }
    if (epsilon2 < 0.0) {
      throw std::domain_error("Parameter 'epsilon2' must be >= 0");
    }
    if (r_cut < 0.0) {
      throw std::domain_error("Parameter 'r_cut' must be >= 0");
    }

    set_prefactor(prefactor);
    this->kappa = kappa;
    this->epsilon1 = epsilon1;
    this->epsilon2 = epsilon2;
    this->r_cut = r_cut;
    B = (2. * (epsilon1 - epsilon2) * (1. + kappa * r_cut) -
         epsilon2 * kappa * kappa * r_cut * r_cut) /
        ((epsilon1 + 2. * epsilon2) * (1. + kappa * r_cut) +
         epsilon2 * kappa * kappa * r_cut * r_cut);
  }

  void on_activation() const { sanity_checks(); }
  void on_boxl_change() const {}
  void on_node_grid_change() const {}
  void on_periodicity_change() const {}
  void on_cell_structure_change() const {}
  void init() const {}

  void sanity_checks() const { sanity_checks_charge_neutrality(); }

  /** @brief Compute the Reaction Field pair force.
   *  @param[in]  q1q2      Product of the charges on p1 and p2.
   *  @param[in]  d         Vector pointing from p1 to p2.
   *  @param[in]  dist      Distance between p1 and p2.
   */
  Utils::Vector3d pair_force(double const q1q2, Utils::Vector3d const &d,
                             double const dist) const {
    if (dist >= r_cut) {
      return {};
    }
    auto fac = 1. / Utils::int_pow<3>(dist) + B / Utils::int_pow<3>(r_cut);
    return (prefactor * q1q2 * fac) * d;
  }

  /** @brief Compute the Reaction Field pair energy.
   *  @param q1q2      Product of the charges on p1 and p2.
   *  @param dist      Distance between p1 and p2.
   */
  double pair_energy(double const q1q2, double const dist) const {
    if (dist >= r_cut) {
      return 0.;
    }
    auto fac = 1. / dist - (B * dist * dist) / (2. * Utils::int_pow<3>(r_cut));
    // remove discontinuity at dist = r_cut
    fac -= (1. - B / 2.) / r_cut;
    return prefactor * q1q2 * fac;
  }
};

#endif // ELECTROSTATICS

#endif
