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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "electrostatics/p3m.hpp"

#include "communication.hpp"
#include "kokkos_helpers.hpp"
#include "p3m/P3MFFTBackend.hpp"
#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/send_mesh.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <Kokkos_Core.hpp>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <typename FloatType, Arch Architecture, class FFTConfig> class P3MFFT;

/**
 * @brief Base class for the electrostatics P3M algorithm.
 * Contains a handle to the FFT backend, information about the local
 * mesh, the differential operator, and various buffers.
 */
template <typename FloatType, Arch Architecture, class FFTConfig>
struct CoulombP3MState : public P3MStateCommon<FloatType, Architecture> {
  using P3MStateCommon<FloatType, Architecture>::P3MStateCommon;
  using value_type = FloatType;
  using ComplexType = std::complex<value_type>;
  using memory_space = Kokkos::HostSpace;
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  using r_space_layout = FFTConfig::r_space_layout;
  using k_space_layout = FFTConfig::k_space_layout;

  /** number of charged particles. */
  std::size_t sum_qpart = 0;
  /** Sum of square of charges. */
  double sum_q2 = 0.;
  /** square of sum of charges. */
  double square_sum_q = 0.;

  p3m_interpolation_cache inter_weights;

  /** charge density in real-space with halo */
  std::vector<FloatType> rs_charge_density;
  /** charge density in k-space without halo */
  std::vector<ComplexType> ks_charge_density;
  /** electric fields in real-space with halo */
  std::array<std::vector<FloatType>, 3> rs_E_fields;
  /** electric fields in k-space without halo */
  std::array<std::vector<ComplexType>, 3> ks_E_fields;
  /** electric fields in real-space without halo */
  std::array<std::vector<FloatType>, 3> rs_E_fields_no_halo;
  p3m_send_mesh<FloatType> halo_comm;
  // FFT reciprocal-space transform, hidden behind the backend interface so the
  // solver is agnostic to heFFTe vs kokkos-fft (see init_cpu_kernels).
  std::shared_ptr<P3MFFTBackend<FloatType, FFTConfig>> fft;
  Kokkos::View<FloatType **, r_space_layout, Kokkos::HostSpace>
      rs_charge_density_kokkos;

  void init_labels() {
    assert(rs_charge_density.empty());
    rs_charge_density_kokkos = decltype(rs_charge_density_kokkos)(
        "CoulombP3MState::rs_charge_density_kokkos", 0, 0);
  }
};

#ifdef ESPRESSO_CUDA
struct P3MGpuParams;
#endif

template <typename FloatType, Arch Architecture, class FFTConfig>
struct CoulombP3MImpl : public CoulombP3M {
  ~CoulombP3MImpl() override = default;

  using CoulombP3MStateClass =
      CoulombP3MState<FloatType, Architecture, FFTConfig>;
  /** @brief Coulomb P3M parameters. */
  CoulombP3MStateClass &p3m;

private:
  // kokkos handle must outlive kokkos data structures from other class members
  std::shared_ptr<KokkosHandle> m_kokkos_handle;
  std::unique_ptr<CoulombP3MStateClass> p3m_state_ptr;
  TuningParameters tuning;
  bool m_is_tuned;

  constexpr const Utils::Vector3i get_memory_layout() const {
    auto constexpr memory_order = CoulombP3MStateClass::memory_order;
    if constexpr (memory_order == Utils::MemoryOrder::COLUMN_MAJOR) {
      return {2, 1, 0};
    }
    return {0, 1, 2};
  }

public:
  CoulombP3MImpl(std::unique_ptr<CoulombP3MStateClass> &&p3m_state,
                 TuningParameters tuning_params, double prefactor)
      : CoulombP3M(p3m_state->params), p3m{*p3m_state},
        m_kokkos_handle{::kokkos_handle}, p3m_state_ptr{std::move(p3m_state)},
        tuning{std::move(tuning_params)} {

    if (tuning.timings <= 0) {
      throw std::domain_error("Parameter 'timings' must be > 0");
    }
    m_is_tuned = not p3m.params.tuning;
    p3m.params.tuning = false;
    set_prefactor(prefactor);
    p3m.init_labels();
  }

  void init() override {
    if constexpr (Architecture == Arch::CPU) {
      init_cpu_kernels();
    }
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::CUDA) {
      init_gpu_kernels();
    }
#endif
  }
  void tune() override;
  void count_charged_particles() override;
  void count_charged_particles_elc(std::size_t n, double sum_q2,
                                   double square_sum_q) override {
    p3m.sum_qpart = n;
    p3m.sum_q2 = sum_q2;
    p3m.square_sum_q = square_sum_q;
  }
  void adapt_epsilon_elc() override {
    p3m.params.epsilon = P3M_EPSILON_METALLIC;
  }

  [[nodiscard]] bool is_tuned() const noexcept override { return m_is_tuned; }
  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture != Arch::CPU;
  }
  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

  void on_activation() override {
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::CUDA) {
      request_gpu();
    }
#endif
    sanity_checks();
    tune();
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::CUDA) {
      if (is_tuned()) {
        init_cpu_kernels();
      }
    }
#endif
  }

  double long_range_energy() override { return long_range_kernel(false, true); }

  void add_long_range_forces() override {
    if constexpr (Architecture == Arch::CPU) {
      long_range_kernel(true, false);
    }
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::CUDA) {
      add_long_range_forces_gpu();
    }
#endif
  }

  Utils::Vector9d long_range_pressure() override;

  void charge_assign() override;
  void assign_charge(double q, Utils::Vector3d const &real_pos,
                     bool skip_cache) override;
  void prepare_fft_mesh(bool reset_weights) override {
    if (reset_weights) {
      p3m.inter_weights.reset(p3m.params.cao);
    }
    p3m.rs_charge_density.resize(p3m.local_mesh.size);
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    auto const num_threads = execution_space().concurrency();
    Kokkos::realloc(Kokkos::WithoutInitializing, p3m.rs_charge_density_kokkos,
                    num_threads, p3m.local_mesh.size);
    kokkos_deep_copy(execution_space{}, p3m.rs_charge_density_kokkos,
                     FloatType{0});
    std::ranges::fill(p3m.rs_charge_density, FloatType{0});
  }

protected:
  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag);
  void calc_influence_function_force() override;
  void calc_influence_function_energy() override;
  void scaleby_box_l() override;
  void init_cpu_kernels();
#ifdef ESPRESSO_CUDA
  void init_gpu_kernels();
  void add_long_range_forces_gpu();
  std::shared_ptr<P3MGpuParams> m_gpu_data = nullptr;
  void request_gpu() const;
#endif
private:
  void kernel_ks_charge_density();
  void kernel_rs_electric_field();
};

#endif // ESPRESSO_P3M
