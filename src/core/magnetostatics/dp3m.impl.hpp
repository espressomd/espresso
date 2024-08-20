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
 *  P3M algorithm for long-range magnetic dipole-dipole interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be
 *  found in @cite hockney88a and @cite deserno98a @cite deserno98b.
 *
 *  Further reading: @cite cerda08d
 */

#pragma once

#include "config/config.hpp"

#ifdef DP3M

#include "magnetostatics/dp3m.hpp"

#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/interpolation.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/math/AS_erfc_part.hpp>

#include <memory>
#include <utility>
#include <vector>

template <typename FloatType>
struct p3m_data_struct_dipoles : public p3m_data_struct_fft<FloatType> {
  using p3m_data_struct_fft<FloatType>::p3m_data_struct_fft;

  /** number of dipolar particles (only on head node). */
  int sum_dip_part = 0;
  /** Sum of square of magnetic dipoles (only on head node). */
  double sum_mu2 = 0.;

  /** position shift for calculation of first assignment mesh point. */
  double pos_shift = 0.;

  /** cached k-space self-energy correction */
  double energy_correction = 0.;
  /** k-space scalar mesh for k-space calculations. */
  std::vector<FloatType> ks_scalar;

  p3m_interpolation_cache inter_weights;
};

template <typename FloatType, Arch Architecture>
struct DipolarP3MImpl : public DipolarP3M {
  ~DipolarP3MImpl() override = default;

  /** @brief Coulomb P3M meshes and FFT algorithm. */
  std::unique_ptr<p3m_data_struct_dipoles<FloatType>> dp3m_impl;
  /** @brief Dipolar P3M parameters. */
  p3m_data_struct_dipoles<FloatType> &dp3m;

  template <typename... Args>
  DipolarP3MImpl(
      std::unique_ptr<p3m_data_struct_dipoles<FloatType>> &&dp3m_handle,
      Args &&...args)
      : DipolarP3M(*dp3m_handle, std::forward<Args>(args)...),
        dp3m_impl{std::move(dp3m_handle)}, dp3m{*dp3m_impl} {}

  void init() override {
    if constexpr (Architecture == Arch::CPU) {
      init_cpu_kernels();
    }
  }
  void tune() override;
  void count_magnetic_particles() override;

  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == Arch::GPU;
  }

  void on_activation() override {
    sanity_checks();
    tune();
  }

  double long_range_energy(ParticleRange const &particles) override {
    return long_range_kernel(false, true, particles);
  }

  void add_long_range_forces(ParticleRange const &particles) override {
    if constexpr (Architecture == Arch::CPU) {
      long_range_kernel(true, false, particles);
    }
  }

  void dipole_assign(ParticleRange const &particles) override;

protected:
  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag,
                           ParticleRange const &particles);
  double calc_average_self_energy_k_space() const override;
  void calc_energy_correction() override;
  void calc_influence_function_force() override;
  void calc_influence_function_energy() override;
  double calc_surface_term(bool force_flag, bool energy_flag,
                           ParticleRange const &particles) override;
  void init_cpu_kernels();
  void scaleby_box_l() override;
};

template <typename FloatType, Arch Architecture,
          template <typename> class FFTBackendImpl, class... Args>
std::shared_ptr<DipolarP3M> new_dp3m_handle(P3MParameters &&p3m,
                                            Args &&...args) {
  auto obj = std::make_shared<DipolarP3MImpl<FloatType, Architecture>>(
      std::make_unique<p3m_data_struct_dipoles<FloatType>>(std::move(p3m)),
      std::forward<Args>(args)...);
  obj->dp3m.template make_fft_instance<FFTBackendImpl<FloatType>>();
  return obj;
}

#endif // DP3M
