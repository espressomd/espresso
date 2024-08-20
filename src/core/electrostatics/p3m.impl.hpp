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

#ifdef P3M

#include "electrostatics/p3m.hpp"

#include "p3m/common.hpp"
#include "p3m/data_struct.hpp"
#include "p3m/interpolation.hpp"

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

template <typename FloatType>
struct p3m_data_struct_coulomb : public p3m_data_struct<FloatType> {
  using p3m_data_struct<FloatType>::p3m_data_struct;

  /** number of charged particles (only on head node). */
  int sum_qpart = 0;
  /** Sum of square of charges (only on head node). */
  double sum_q2 = 0.;
  /** square of sum of charges (only on head node). */
  double square_sum_q = 0.;

  p3m_interpolation_cache inter_weights;
};

#ifdef CUDA
struct P3MGpuParams;
#endif

template <typename FloatType, Arch Architecture>
struct CoulombP3MImpl : public CoulombP3M {
  ~CoulombP3MImpl() override = default;

  /** @brief Coulomb P3M parameters. */
  p3m_data_struct_coulomb<FloatType> &p3m;

private:
  /** @brief Coulomb P3M meshes and FFT algorithm. */
  std::unique_ptr<p3m_data_struct_coulomb<FloatType>> p3m_impl;
  int tune_timings;
  bool tune_verbose;
  bool check_complex_residuals;
  bool m_is_tuned;

public:
  CoulombP3MImpl(
      std::unique_ptr<p3m_data_struct_coulomb<FloatType>> &&p3m_handle,
      double prefactor, int tune_timings, bool tune_verbose,
      bool check_complex_residuals)
      : CoulombP3M(p3m_handle->params), p3m{*p3m_handle},
        p3m_impl{std::move(p3m_handle)}, tune_timings{tune_timings},
        tune_verbose{tune_verbose},
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
#ifdef CUDA
    if constexpr (Architecture == Arch::GPU) {
      init_gpu_kernels();
    }
#endif
  }
  void tune() override;
  void count_charged_particles() override;
  void count_charged_particles_elc(int n, double sum_q2,
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
#ifdef CUDA
    if constexpr (Architecture == Arch::GPU) {
      request_gpu();
    }
#endif
    sanity_checks();
    tune();
#ifdef CUDA
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
#ifdef CUDA
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
    for (int i = 0; i < p3m.local_mesh.size; i++) {
      p3m.mesh.rs_scalar[i] = FloatType(0);
    }
  }

protected:
  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag,
                           ParticleRange const &particles);
  void calc_influence_function_force() override;
  void calc_influence_function_energy() override;
  void scaleby_box_l() override;
  void init_cpu_kernels();
#ifdef CUDA
  void init_gpu_kernels();
  void add_long_range_forces_gpu(ParticleRange const &particles);
  std::shared_ptr<P3MGpuParams> m_gpu_data = nullptr;
  void request_gpu() const;
#endif
};

template <typename FloatType, Arch Architecture,
          template <typename> class FFTBackendImpl,
          template <typename> class P3MFFTMeshImpl, class... Args>
std::shared_ptr<CoulombP3M> new_p3m_handle(P3MParameters &&p3m,
                                           Args &&...args) {
  auto obj = std::make_shared<CoulombP3MImpl<FloatType, Architecture>>(
      std::make_unique<p3m_data_struct_coulomb<FloatType>>(std::move(p3m)),
      std::forward<Args>(args)...);
  obj->p3m.template make_mesh_instance<P3MFFTMeshImpl<FloatType>>();
  obj->p3m.template make_fft_instance<FFTBackendImpl<FloatType>>();
  return obj;
}

#endif // P3M
