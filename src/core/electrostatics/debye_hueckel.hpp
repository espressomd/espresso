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
 * Calculate the Debye-Hückel energy and force for a particle pair.
 */

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_DEBYE_HUECKEL_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_DEBYE_HUECKEL_HPP

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/actor.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <cmath>
#include <stdexcept>

/** @brief Debye-Hückel parameters. */
struct DebyeHueckel : public Coulomb::Actor<DebyeHueckel> {
  /** @brief Ionic strength. */
  double kappa;
  /** @brief Interaction cutoff. */
  double r_cut;

  DebyeHueckel(double prefactor, double kappa, double r_cut) {
    if (kappa < 0.0) {
      throw std::domain_error("Parameter 'kappa' must be >= 0");
    }
    if (r_cut < 0.0) {
      throw std::domain_error("Parameter 'r_cut' must be >= 0");
    }

    set_prefactor(prefactor);
    this->kappa = kappa;
    this->r_cut = r_cut;
  }

  void on_activation() const { sanity_checks(); }
  void on_boxl_change() const {}
  void on_node_grid_change() const {}
  void on_periodicity_change() const {}
  void on_cell_structure_change() const {}
  void init() const {}

  void sanity_checks() const { sanity_checks_charge_neutrality(); }

  /** @brief Compute the pair force.
   *  @param[in]  q1q2      Product of the charges on p1 and p2.
   *  @param[in]  d         Vector pointing from p1 to p2.
   *  @param[in]  dist      Distance between p1 and p2.
   */
  Utils::Vector3d pair_force(double const q1q2, Utils::Vector3d const &d,
                             double const dist) const {
    if (dist >= r_cut) {
      return {};
    }
    // pure Coulomb case
    auto fac = prefactor * q1q2 / Utils::int_pow<3>(dist);
    if (kappa > 0.) {
      // Debye-Hueckel case
      fac *= std::exp(-kappa * dist) * (1. + kappa * dist);
    }
    return fac * d;
  }

  /** @brief Compute the pair energy.
   *  @param q1q2      Product of the charges on p1 and p2.
   *  @param dist      Distance between p1 and p2.
   */
  double pair_energy(double const q1q2, double const dist) const {
    if (dist >= r_cut) {
      return 0.;
    }
    // pure Coulomb case
    auto energy = prefactor * q1q2 / dist;
    if (kappa > 0.) {
      // Debye-Hueckel case
      energy *= std::exp(-kappa * dist);
    }
    return energy;
  }
};

#endif // ELECTROSTATICS

#endif
