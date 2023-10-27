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
 * MMM1D algorithm for long-range %Coulomb interactions on the CPU.
 * Implementation of the MMM1D method for the calculation of the electrostatic
 * interaction in one-dimensionally periodic systems. For details on the
 * method see MMM in general. The MMM1D method works only with the N-squared
 * cell system since neither the near nor far formula can be decomposed.
 */

#pragma once

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/actor.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <array>

/** @brief Parameters for the MMM1D electrostatic interaction */
struct CoulombMMM1D : public Coulomb::Actor<CoulombMMM1D> {
  /**
   * @brief Maximal allowed pairwise error for the potential and force.
   * This error ignores prefactors, i.e. this is for a pure lattice 1/r-sum.
   */
  double maxPWerror;
  /**
   * @brief Far switch radius. Represents the xy-distance at which
   * the calculation switches from the far to the near formula.
   */
  double far_switch_radius;
  int tune_timings;
  bool tune_verbose;

  CoulombMMM1D(double prefactor, double maxPWerror, double switch_rad,
               int tune_timings, bool tune_verbose);

  /** Compute the pair force.
   *  @param[in]  q1q2      Product of the charges on p1 and p2.
   *  @param[in]  d         Vector pointing from p1 to p2.
   *  @param[in]  dist      Distance between p1 and p2.
   */
  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const;

  /** Compute the pair energy.
   *  @param[in] q1q2      Product of the charges on p1 and p2.
   *  @param[in] d         Vector pointing from p1 to p2.
   *  @param[in] dist      Distance between p1 and p2.
   */
  double pair_energy(double q1q2, Utils::Vector3d const &d, double dist) const;

  void tune();
  bool is_tuned() const { return m_is_tuned; }

  void on_activation() {
    sanity_checks();
    tune();
  }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { recalc_boxl_parameters(); }
  void on_node_grid_change() const {}
  void on_periodicity_change() const { sanity_checks_periodicity(); }
  void on_cell_structure_change() { sanity_checks_cell_structure(); }
  /** @brief Recalculate all derived parameters. */
  void init() { recalc_boxl_parameters(); }

  void sanity_checks() const {
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
  }

private:
  bool m_is_tuned;
  /** @brief Square of the far switch radius. */
  double far_switch_radius_sq;
  /** @brief Squared inverse box length in z-direction. */
  double uz2;
  /** @brief Squared inverse box length in z-direction times prefactor. */
  double prefuz2;
  /** @brief Cubed inverse box length in z-direction times prefactor. */
  double prefL3_i;
  /**
   * @brief Largest numerically stable cutoff for Bessel function.
   * Don't change without improving the formulas.
   */
  static constexpr auto MAXIMAL_B_CUT = 30;
  /** @brief From which distance a certain Bessel cutoff is valid. */
  std::array<double, MAXIMAL_B_CUT> bessel_radii;

  void determine_bessel_radii();
  void prepare_polygamma_series();
  void recalc_boxl_parameters();
  void sanity_checks_periodicity() const;
  void sanity_checks_cell_structure() const;
};

#endif // ELECTROSTATICS
