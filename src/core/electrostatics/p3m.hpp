/*
 * Copyright (C) 2010-2024 The ESPResSo project
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

/** @file
 *  P3M algorithm for long-range Coulomb interaction.
 *
 *  We use a P3M (Particle-Particle Particle-Mesh) method based on the
 *  Ewald summation. Details of the used method can be found in
 *  @cite hockney88a and @cite deserno98a @cite deserno98b.
 *
 *  Further reading: @cite ewald21a, @cite hockney88a, @cite deserno98a,
 *  @cite deserno98b, @cite deserno00e, @cite deserno00b, @cite cerda08d.
 *
 *  Implementation in p3m.cpp.
 */

#pragma once

#include "config/config.hpp"

#ifdef P3M

#include "electrostatics/actor.hpp"

#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/math/AS_erfc_part.hpp>

#include <cmath>
#include <numbers>

/** @brief P3M solver. */
struct CoulombP3M : public Coulomb::Actor<CoulombP3M> {
  P3MParameters const &p3m_params;

public:
  CoulombP3M(P3MParameters const &p3m_params) : p3m_params{p3m_params} {}

  virtual ~CoulombP3M() = default;

  [[nodiscard]] virtual bool is_tuned() const noexcept = 0;
  [[nodiscard]] virtual bool is_gpu() const noexcept = 0;
  [[nodiscard]] virtual bool is_double_precision() const noexcept = 0;

  /** @brief Recalculate all derived parameters. */
  virtual void init() = 0;
  virtual void on_activation() = 0;

  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { scaleby_box_l(); }
  void on_node_grid_change() const { sanity_checks_node_grid(); }
  void on_periodicity_change() const { sanity_checks_periodicity(); }
  void on_cell_structure_change() {
    sanity_checks_cell_structure();
    init();
  }
  void sanity_checks() const {
    sanity_checks_boxl();
    sanity_checks_node_grid();
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
  }

  /**
   * Count the number of charged particles and calculate
   * the sum of the squared charges.
   */
  virtual void count_charged_particles() = 0;
  virtual void count_charged_particles_elc(int, double, double) = 0;
  virtual void adapt_epsilon_elc() = 0;

  /**
   * @brief Tune P3M parameters to desired accuracy.
   *
   * The parameters
   * @ref P3MParameters::mesh "mesh",
   * @ref P3MParameters::cao "cao",
   * @ref P3MParameters::r_cut_iL "r_cut_iL" and
   * @ref P3MParameters::alpha_L "alpha_L" are tuned to obtain the target
   * @ref P3MParameters::accuracy "accuracy" in optimal time.
   * These parameters are stored in the @ref p3m_params object.
   *
   * The function utilizes the analytic expression of the error estimate
   * for the P3M method in @cite hockney88a (eq. (8.23)) in
   * order to obtain the rms error in the force for a system of N randomly
   * distributed particles in a cubic box.
   * For the real space error the estimate of Kolafa/Perram is used.
   *
   * Parameter ranges if not given explicitly in the constructor:
   * - @p mesh explores the range 0-512 (biased toward values less than 128)
   * - @p cao explores all possible values
   * - @p alpha_L is tuned for each tuple (@p r_cut_iL, @p mesh, @p cao) and
   *   calculated assuming that the error contributions of real and reciprocal
   *   space should be equal
   *
   * After checking if the total error lies below the target accuracy,
   * the time needed for one force calculation is measured. Parameters
   * that minimize the runtime are kept.
   *
   * The function is based on routines of the program HE_Q.cpp written by M.
   * Deserno.
   */
  virtual void tune() = 0;

  /** Assign the physical charges using the tabulated charge assignment
   * function.
   */
  virtual void charge_assign(ParticleRange const &particles) = 0;

  /**
   * @brief Assign a single charge into the current charge grid.
   *
   * @param[in] q          Particle charge
   * @param[in] real_pos   Particle position in real space
   * @param[in] skip_cache Skip interpolation weights caching.
   */
  virtual void assign_charge(double q, Utils::Vector3d const &real_pos,
                             bool skip_cache) = 0;

  virtual void prepare_fft_mesh(bool reset_weights) = 0;

  /** Calculate real-space contribution of p3m Coulomb pair forces. */
  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const {
    if ((q1q2 == 0.) || dist >= p3m_params.r_cut || dist <= 0.) {
      return {};
    }
    auto const alpha = p3m_params.alpha;
    auto const adist = alpha * dist;
    auto const exp_adist_sq = exp(-adist * adist);
    auto const dist_sq = dist * dist;
    auto const two_a_sqrt_pi_i = 2. * alpha * std::numbers::inv_sqrtpi;
#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    auto const fac = exp_adist_sq * (erfc_part_ri + two_a_sqrt_pi_i) / dist_sq;
#else
    auto const erfc_part_ri = erfc(adist) / dist;
    auto const fac = (erfc_part_ri + two_a_sqrt_pi_i * exp_adist_sq) / dist_sq;
#endif
    return (fac * prefactor * q1q2) * d;
  }

  /** Calculate real-space contribution of Coulomb pair energy. */
  // Eq. (3.6) @cite deserno00b
  double pair_energy(double q1q2, double dist) const {
    if ((q1q2 == 0.) || dist >= p3m_params.r_cut || dist <= 0.) {
      return {};
    }
    auto const adist = p3m_params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    return prefactor * q1q2 * erfc_part_ri * exp(-adist * adist);
#else
    auto const erfc_part_ri = erfc(adist) / dist;
    return prefactor * q1q2 * erfc_part_ri;
#endif
  }

  /** Compute the k-space part of the pressure tensor */
  virtual Utils::Vector9d long_range_pressure(ParticleRange const &) = 0;

  /** Compute the k-space part of energies. */
  virtual double long_range_energy(ParticleRange const &) = 0;

  /** Compute the k-space part of forces. */
  virtual void add_long_range_forces(ParticleRange const &) = 0;

protected:
  virtual void calc_influence_function_force() = 0;
  virtual void calc_influence_function_energy() = 0;

  /** Checks for correctness of the k-space cutoff. */
  void sanity_checks_boxl() const;
  void sanity_checks_node_grid() const;
  void sanity_checks_periodicity() const;
  void sanity_checks_cell_structure() const;

  virtual void scaleby_box_l() = 0;
};

#endif // P3M
