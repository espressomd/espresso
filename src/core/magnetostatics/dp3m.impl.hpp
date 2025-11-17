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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_DP3M

#include "magnetostatics/dp3m.hpp"

#include "ParticleRange.hpp"
#include "communication.hpp"
#include "p3m/FFTBackendLegacy.hpp"
#include "p3m/FFTBuffersLegacy.hpp"
#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/interpolation.hpp"

#include <utils/Vector.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#include <omp.h>
#endif

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <typename FloatType>
struct p3m_data_struct_dipoles : public p3m_data_struct<FloatType> {
  using p3m_data_struct<FloatType>::p3m_data_struct;

  /** number of dipolar particles. */
  std::size_t sum_dip_part = 0;
  /** Sum of square of magnetic dipoles. */
  double sum_mu2 = 0.;

  /** position shift for calculation of first assignment mesh point. */
  double pos_shift = 0.;

  /** cached k-space self-energy correction */
  double energy_correction = 0.;
  /** k-space scalar mesh for k-space calculations. */
  std::vector<FloatType> ks_scalar;

  p3m_interpolation_cache inter_weights;

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace>
      rs_fields_kokkos;

  void init_labels() {
    assert(ks_scalar.empty());
    rs_fields_kokkos = decltype(rs_fields_kokkos)(
        "p3m_data_struct_dipoles::rs_fields_kokkos", 0, 0, 0);
  }
#endif
};

template <typename FloatType, Arch Architecture>
struct DipolarP3MImpl : public DipolarP3M {
  ~DipolarP3MImpl() override = default;

  /** @brief Dipolar P3M parameters. */
  p3m_data_struct_dipoles<FloatType> &dp3m;

private:
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  // kokkos handle must outlive kokkos data structures from other class members
  std::shared_ptr<KokkosHandle> m_kokkos_handle;
#endif
  /** @brief Dipolar P3M meshes and FFT algorithm. */
  std::unique_ptr<p3m_data_struct_dipoles<FloatType>> dp3m_impl;
  TuningParameters tuning;
  bool m_is_tuned;

public:
  DipolarP3MImpl(
      std::unique_ptr<p3m_data_struct_dipoles<FloatType>> &&dp3m_handle,
      TuningParameters tuning_params, double prefactor)
      : DipolarP3M(dp3m_handle->params), dp3m{*dp3m_handle},
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
        m_kokkos_handle{::kokkos_handle},
#endif
        dp3m_impl{std::move(dp3m_handle)}, tuning{std::move(tuning_params)} {

    if (tuning.timings <= 0) {
      throw std::domain_error("Parameter 'timings' must be > 0");
    }
    if (dp3m.params.mesh != Utils::Vector3i::broadcast(dp3m.params.mesh[0])) {
      throw std::domain_error("DipolarP3M requires a cubic mesh");
    }
    m_is_tuned = not dp3m.params.tuning;
    dp3m.params.tuning = false;
    set_prefactor(prefactor);
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    dp3m.init_labels();
#endif
  }

  void init() override {
    if constexpr (Architecture == Arch::CPU) {
      init_cpu_kernels();
    }
  }
  void tune() override;
  void count_magnetic_particles() override;

  [[nodiscard]] bool is_tuned() const noexcept override { return m_is_tuned; }
  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture != Arch::CPU;
  }
  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
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

private:
  void prepare_fft_mesh() {
    dp3m.inter_weights.reset(dp3m.params.cao);
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    using execution_space = Kokkos::DefaultExecutionSpace;
    auto const num_threads = execution_space().concurrency();
    Kokkos::realloc(Kokkos::WithoutInitializing, dp3m.rs_fields_kokkos,
                    num_threads, 3, dp3m.local_mesh.size);
    Kokkos::deep_copy(dp3m.rs_fields_kokkos, FloatType{0});
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
    for (auto &rs_mesh_field : dp3m.mesh.rs_fields) {
      std::ranges::fill(rs_mesh_field, FloatType{0});
    }
  }

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
#ifdef ESPRESSO_NPT
  void npt_add_virial_contribution(double energy) const override;
#endif
};

#endif // ESPRESSO_DP3M
