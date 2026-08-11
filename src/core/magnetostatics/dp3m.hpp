/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
 *  P3M algorithm for long-range magnetic dipole-dipole interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be
 *  found in @cite hockney88a and @cite deserno98a @cite deserno98b.
 *
 *  Further reading: @cite cerda08d
 */

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_DP3M

#include "magnetostatics/actor.hpp"

#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/math.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/math/AS_erfc_part.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <numbers>

/** @brief Dipolar P3M solver. */
struct DipolarP3M : public Dipoles::Actor<DipolarP3M> {
  P3MParameters const &dp3m_params;

public:
  DipolarP3M(P3MParameters const &dp3m_params) : dp3m_params{dp3m_params} {}

  virtual ~DipolarP3M() = default;

  [[nodiscard]] virtual bool is_tuned() const noexcept = 0;
  [[nodiscard]] virtual bool is_gpu() const noexcept = 0;
  [[nodiscard]] virtual bool is_double_precision() const noexcept = 0;

  virtual void on_activation() = 0;
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { scaleby_box_l(); }
  void on_node_grid_change() const { sanity_checks_node_grid(); }
  void on_periodicity_change() const { sanity_checks_periodicity(); }
  void on_cell_structure_change() {
    sanity_checks_cell_structure();
    init();
  }
  /** @brief Recalculate all derived parameters. */
  virtual void init() = 0;

  void sanity_checks() const {
    sanity_checks_boxl();
    sanity_checks_node_grid();
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
  }

  /**
   * @brief Count the number of magnetic particles and calculate
   * the sum of the squared dipole moments.
   */
  virtual void count_magnetic_particles() = 0;

  /** Assign the physical dipoles using the tabulated assignment function. */
  virtual void dipole_assign() = 0;

  /**
   * @brief Tune dipolar P3M parameters to desired accuracy.
   *
   * The parameters
   * @ref P3MParameters::mesh "mesh",
   * @ref P3MParameters::cao "cao",
   * @ref P3MParameters::r_cut_iL "r_cut_iL" and
   * @ref P3MParameters::alpha_L "alpha_L" are tuned to obtain the target
   * @ref P3MParameters::accuracy "accuracy" in optimal time.
   * These parameters are stored in the @ref dp3m_params object.
   *
   * The function utilizes the analytic expression of the error estimate
   * for the dipolar P3M method in the paper of @cite cerda08d in
   * order to obtain the rms error in the force for a system of N randomly
   * distributed particles in a cubic box. For the real space error, the
   * estimate in @cite kolafa92a is used.
   *
   * Parameter ranges if not given explicitly in the constructor:
   * - @p mesh is set up such that the number of mesh points is equal to the
   *   number of magnetic dipolar particles
   * - @p cao explores all possible values
   * - @p alpha_L is tuned for each tuple (@p r_cut_iL, @p mesh, @p cao) and
   *   calculated assuming that the error contributions of real and reciprocal
   *   space should be equal
   *
   * After checking if the total error lies below the target accuracy,
   * the time needed for one force calculation is measured. Parameters
   * that minimize the runtime are kept.
   *
   * The function is based on routines of the program HE_Q.cpp for charges
   * written by M. Deserno.
   */
  virtual void tune() = 0;

  /** Compute the k-space part of the pressure tensor */
  virtual Utils::Vector9d long_range_pressure() = 0;

  /** Compute the k-space part of energies. */
  virtual double long_range_energy() = 0;

  /** Compute the k-space part of forces. */
  virtual void add_long_range_forces() = 0;

  /** Calculate real-space contribution of p3m dipolar pair forces and torques.
   *  If NPT is compiled in, update the NpT virial.
   */
  inline ParticleForce pair_force(double d1d2, Utils::Vector3d const &dip1,
                                  Utils::Vector3d const &dip2,
                                  Utils::Vector3d const &d, double dist,
                                  double dist2) const {
    if (d1d2 == 0. or dist >= dp3m_params.r_cut or dist <= 0.)
      return {};

    auto const alpsq = dp3m_params.alpha * dp3m_params.alpha;
    auto const adist = dp3m_params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
#else
    auto const erfc_part_ri = erfc(adist) / dist;
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    auto const mimj = dip1 * dip2;

    auto const mir = dip1 * d;
    auto const mjr = dip2 * d;

    auto const coeff = 2. * dp3m_params.alpha * std::numbers::inv_sqrtpi;
    auto const dist2i = 1. / dist2;
    auto const exp_adist2 = exp(-Utils::sqr(adist));

    auto const B_r = (dp3m_params.accuracy > 5e-06)
                         ? (erfc_part_ri + coeff) * exp_adist2 * dist2i
                         : (erfc(adist) / dist + coeff * exp_adist2) * dist2i;

    auto const common_term = alpsq * coeff * exp_adist2;
    auto const C_r = dist2i * (3. * B_r + 2. * common_term);
    auto const D_r = dist2i * (5. * C_r + 4. * common_term * alpsq);

    // Calculate real-space forces
    auto const force = prefactor * ((mimj * d + dip1 * mjr + dip2 * mir) * C_r -
                                    mir * mjr * D_r * d);

    // Calculate vector multiplications for vectors mi, mj, rij
    auto const mixmj = vector_product(dip1, dip2);
    auto const mixr = vector_product(dip1, d);

    // Calculate real-space torques
    auto const torque = prefactor * (-mixmj * B_r + mixr * (mjr * C_r));
    // NOTE: the NpT virial contribution of this pair force (d * force,
    // the trace of the pairwise virial tensor d (x) force, computed
    // explicitly rather than approximated from the pair energy -- see
    // forces_cabana.hpp for why) is accumulated by the caller (see
    // ForcesKernel in forces_cabana.hpp), not here:
    // pair_force() is invoked concurrently by many threads over all pairs
    // (both for force calculation and for the pressure observable via
    // pair_pressure_kernel() in dipoles_inline.hpp), so calling
    // npt_add_virial_contribution() from here would race on the shared
    // NpT virial accumulator.
    return ParticleForce{force, torque};
  }

  /** Calculate real-space contribution of dipolar pair energy. */
  inline double pair_energy(Utils::Vector3d const &dip1,
                            Utils::Vector3d const &dip2,
                            Utils::Vector3d const &d, double dist,
                            double dist2) const {
    if (dist >= dp3m_params.r_cut or dist <= 0.)
      return {};

    auto const alpsq = dp3m_params.alpha * dp3m_params.alpha;
    auto const adist = dp3m_params.alpha * dist;

#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
#else
    auto const erfc_part_ri = erfc(adist) / dist;
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    auto const mimj = dip1 * dip2;
    auto const mir = dip1 * d;
    auto const mjr = dip2 * d;

    auto const coeff = 2. * dp3m_params.alpha * std::numbers::inv_sqrtpi;
    auto const dist2i = 1. / dist2;
    auto const exp_adist2 = exp(-Utils::sqr(adist));

    auto const B_r = (dp3m_params.accuracy > 5e-06)
                         ? dist2i * (erfc_part_ri + coeff) * exp_adist2
                         : dist2i * (erfc(adist) / dist + coeff * exp_adist2);
    auto const C_r = (3. * B_r + 2. * alpsq * coeff * exp_adist2) * dist2i;

    return prefactor * (mimj * B_r - mir * mjr * C_r);
  }

protected:
  /** Calculate self-energy in k-space. */
  virtual double calc_average_self_energy_k_space() const = 0;

  /** Calculate energy correction that minimizes the error.
   *  This quantity is only calculated once and is cached.
   */
  virtual void calc_energy_correction() = 0;

  /** Calculate the influence function for the dipolar forces and torques. */
  virtual void calc_influence_function_force() = 0;

  /** Calculate the influence function for the dipolar energy. */
  virtual void calc_influence_function_energy() = 0;

  /** Compute the dipolar surface terms */
  virtual double calc_surface_term(bool force_flag, bool energy_flag) = 0;

  /** Checks for correctness of the k-space cutoff. */
  void sanity_checks_boxl() const;
  void sanity_checks_node_grid() const;
  void sanity_checks_periodicity() const;
  void sanity_checks_cell_structure() const;

  virtual void scaleby_box_l() = 0;

#ifdef ESPRESSO_NPT
  /** Update the NpT virial */
  virtual void npt_add_virial_contribution(double virial) const = 0;
#endif
};

std::shared_ptr<DipolarP3M>
new_dipolar_p3m_heffte(P3MParameters &&p3m_params,
                       TuningParameters const &tuning_params, double prefactor,
                       bool single_precision, Arch arch);

#endif // ESPRESSO_DP3M
