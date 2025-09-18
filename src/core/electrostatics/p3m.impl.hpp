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

#include "config/config.hpp"

#ifdef ESPRESSO_P3M

#include "electrostatics/p3m.hpp"

#include "ParticleRange.hpp"
#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/send_mesh.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>

template <typename FloatType> class P3MFFT;

template <typename FloatType>
struct CoulombP3MState : public P3MStateCommon<FloatType> {
  using P3MStateCommon<FloatType>::P3MStateCommon;

  constexpr static auto memory_order = Utils::MemoryOrder::ROW_MAJOR;

  /** number of charged particles. */
  std::size_t sum_qpart = 0;
  /** Sum of square of charges. */
  double sum_q2 = 0.;
  /** square of sum of charges. */
  double square_sum_q = 0.;

  p3m_interpolation_cache inter_weights;

  /* fields */
  std::vector<FloatType> rs_charge_density;
  std::vector<std::complex<FloatType>> ks_charge_density;
  std::array<std::vector<FloatType>, 3> rs_E_fields;
  std::vector<std::complex<FloatType>> ks_E_fields_storage;
  std::vector<std::complex<FloatType>> rs_E_fields_no_halo;
  p3m_send_mesh<FloatType> halo_comm;
  std::shared_ptr<P3MFFT<FloatType>> fft;
};

#ifdef ESPRESSO_CUDA
struct P3MGpuParams;
#endif

template <typename FloatType, Arch Architecture>
struct CoulombP3MImpl : public CoulombP3M {
  ~CoulombP3MImpl() override = default;

  /** @brief Coulomb P3M parameters. */
  CoulombP3MState<FloatType> &p3m;

private:
  std::unique_ptr<CoulombP3MState<FloatType>> p3m_state_ptr;
  int tune_timings;
  std::pair<std::optional<int>, std::optional<int>> tune_limits;
  bool tune_verbose;
  bool check_complex_residuals;
  bool m_is_tuned;

  constexpr const Utils::Vector3i get_memory_layout() const {
    auto constexpr memory_order =
        std::remove_reference<decltype(p3m)>::type::memory_order;
    if constexpr (memory_order == Utils::MemoryOrder::COLUMN_MAJOR) {
      return {2, 1, 0};
    }
    return {0, 1, 2};
  }

public:
  CoulombP3MImpl(std::unique_ptr<CoulombP3MState<FloatType>> &&p3m_state,
                 double prefactor, int tune_timings, bool tune_verbose,
                 decltype(tune_limits) tune_limits,
                 bool check_complex_residuals)
      : CoulombP3M(p3m_state->params), p3m{*p3m_state},
        p3m_state_ptr{std::move(p3m_state)}, tune_timings{tune_timings},
        tune_limits{std::move(tune_limits)}, tune_verbose{tune_verbose},
        check_complex_residuals{check_complex_residuals} {

    if (tune_timings <= 0) {
      throw std::domain_error("Parameter 'timings' must be > 0");
    }
    m_is_tuned = not p3m.params.tuning;
    p3m.params.tuning = false;
    set_prefactor(prefactor);
  }

  void init() override {
    if constexpr (Architecture == Arch::CPU) {
      init_cpu_kernels();
    }
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::GPU) {
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
    return Architecture == Arch::GPU;
  }
  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

  void on_activation() override {
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::GPU) {
      request_gpu();
    }
#endif
    sanity_checks();
    tune();
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::GPU) {
      if (is_tuned()) {
        init_cpu_kernels();
      }
    }
#endif
  }

  double long_range_energy(ParticleRange const &particles) override {
    return long_range_kernel(false, true, particles);
  }

  void add_long_range_forces(ParticleRange const &particles) override {
    if constexpr (Architecture == Arch::CPU) {
      long_range_kernel(true, false, particles);
    }
#ifdef ESPRESSO_CUDA
    if constexpr (Architecture == Arch::GPU) {
      add_long_range_forces_gpu(particles);
    }
#endif
  }

  Utils::Vector9d long_range_pressure(ParticleRange const &particles) override;

  void charge_assign(ParticleRange const &particles) override;
  void assign_charge(double q, Utils::Vector3d const &real_pos,
                     bool skip_cache) override;
  void prepare_fft_mesh(bool reset_weights) override {
    if (reset_weights) {
      p3m.inter_weights.reset(p3m.params.cao);
    }
    p3m.rs_charge_density.resize(Utils::product(p3m.local_mesh.dim));
    std::ranges::fill(p3m.rs_charge_density, FloatType{});
  }

protected:
  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag,
                           ParticleRange const &particles);
  void calc_influence_function_force() override;
  void calc_influence_function_energy() override;
  void scaleby_box_l() override;
  void init_cpu_kernels();
#ifdef ESPRESSO_CUDA
  void init_gpu_kernels();
  void add_long_range_forces_gpu(ParticleRange const &particles);
  std::shared_ptr<P3MGpuParams> m_gpu_data = nullptr;
  void request_gpu() const;
#endif
private:
  void kernel_ks_charge_density();
  void kernel_rs_electric_field();
};

template <typename FloatType, Arch Architecture, typename... Args>
std::shared_ptr<CoulombP3M> new_coulomb_p3m(P3MParameters &&p3m_params,
                                            Args &&...args) {
  auto state_ptr =
      std::make_unique<CoulombP3MState<FloatType>>(std::move(p3m_params));
  auto obj = std::make_shared<CoulombP3MImpl<FloatType, Architecture>>(
      std::move(state_ptr), std::forward<Args>(args)...);
  return obj;
}

#endif // ESPRESSO_P3M
