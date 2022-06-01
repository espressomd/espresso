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

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_P3M_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_P3M_HPP

#include "config.hpp"

#ifdef P3M

#include "electrostatics/actor.hpp"

#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/fft.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/send_mesh.hpp"

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/AS_erfc_part.hpp>

#include <array>
#include <cmath>

struct p3m_data_struct : public p3m_data_struct_base {
  explicit p3m_data_struct(P3MParameters &&parameters)
      : p3m_data_struct_base{std::move(parameters)} {}

  /** local mesh. */
  P3MLocalMesh local_mesh;
  /** real space mesh (local) for CA/FFT. */
  fft_vector<double> rs_mesh;
  /** mesh (local) for the electric field. */
  std::array<fft_vector<double>, 3> E_mesh;

  /** number of charged particles (only on head node). */
  int sum_qpart = 0;
  /** Sum of square of charges (only on head node). */
  double sum_q2 = 0.;
  /** square of sum of charges (only on head node). */
  double square_sum_q = 0.;

  p3m_interpolation_cache inter_weights;

  /** send/recv mesh sizes */
  p3m_send_mesh sm;

  fft_data_struct fft;
};

/** @brief P3M solver. */
struct CoulombP3M : public Coulomb::Actor<CoulombP3M> {
  /** P3M parameters. */
  p3m_data_struct p3m;

  int tune_timings;
  bool tune_verbose;

private:
  bool m_is_tuned;

public:
  CoulombP3M(P3MParameters &&parameters, double prefactor, int tune_timings,
             bool tune_verbose);

  bool is_tuned() const { return m_is_tuned; }

  /** Compute the k-space part of forces and energies. */
  double kernel(bool force_flag, bool energy_flag,
                ParticleRange const &particles);

  /** @brief Recalculate all derived parameters. */
  void init();
  void on_activation() {
    sanity_checks();
    tune();
  }
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
  void count_charged_particles();

  /**
   * @brief Tune P3M parameters to desired accuracy.
   *
   * The parameters
   * @ref P3MParameters::mesh "mesh",
   * @ref P3MParameters::cao "cao",
   * @ref P3MParameters::r_cut_iL "r_cut_iL" and
   * @ref P3MParameters::alpha_L "alpha_L" are tuned to obtain the target
   * @ref P3MParameters::accuracy "accuracy" in optimal time.
   * These parameters are stored in the @ref p3m object.
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
  void tune();

  /** Assign the physical charges using the tabulated charge assignment
   * function.
   */
  void charge_assign(ParticleRange const &particles);

  /**
   * @brief Assign a single charge into the current charge grid.
   *
   * @param[in] q          %Particle charge
   * @param[in] real_pos   %Particle position in real space
   * @param[out] inter_weights Cached interpolation weights to be used.
   */
  void assign_charge(double q, Utils::Vector3d const &real_pos,
                     p3m_interpolation_cache &inter_weights);
  /** @overload */
  void assign_charge(double q, Utils::Vector3d const &real_pos);

  /** Calculate real-space contribution of p3m Coulomb pair forces. */
  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const {
    if ((q1q2 == 0.) || dist >= p3m.params.r_cut || dist <= 0.) {
      return {};
    }
    auto const adist = p3m.params.alpha * dist;
    auto const exp_adist_sq = exp(-adist * adist);
    auto const dist_sq = dist * dist;
    auto const two_a_sqrt_pi_i = 2.0 * p3m.params.alpha * Utils::sqrt_pi_i();
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
  double pair_energy(double q1q2, double dist) const {
    if ((q1q2 == 0.) || dist >= p3m.params.r_cut || dist <= 0.) {
      return {};
    }
    auto const adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    return prefactor * q1q2 * erfc_part_ri * exp(-adist * adist);
#else
    auto const erfc_part_ri = erfc(adist) / dist;
    return prefactor * q1q2 * erfc_part_ri;
#endif
  }

  /** Compute the k-space part of the pressure tensor */
  Utils::Vector9d p3m_calc_kspace_pressure_tensor();

  /** Compute the k-space part of energies. */
  double long_range_energy(ParticleRange const &particles) {
    return long_range_kernel(false, true, particles);
  }

  /** Compute the k-space part of forces. */
  void add_long_range_forces(ParticleRange const &particles) {
    long_range_kernel(true, false, particles);
  }

  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag,
                           ParticleRange const &particles);

private:
  void calc_influence_function_force();
  void calc_influence_function_energy();

  /** Checks for correctness of the k-space cutoff. */
  void sanity_checks_boxl() const;
  void sanity_checks_node_grid() const;
  void sanity_checks_periodicity() const;
  void sanity_checks_cell_structure() const;

  void scaleby_box_l();
};

#endif // P3M

#endif
