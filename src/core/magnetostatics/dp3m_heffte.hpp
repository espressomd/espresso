/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef ESPRESSO_ADDITIONAL_CHECKS
#define ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
#endif

template <typename FloatType, class FFTConfig> class P3MFFT;

/**
 * @brief Base class for the magnetostatics P3M algorithm.
 * Contains a handle to the FFT backend, information about the local
 * mesh, the differential operator, and various buffers.
 */
template <typename FloatType, class FFTConfig>
struct DipolarP3MState : public P3MStateCommon<FloatType> {
  using value_type = FloatType;
  using ComplexType = std::complex<FloatType>;
  using P3MStateCommon<FloatType>::P3MStateCommon;
  using P3MStateCommon<FloatType>::local_mesh;

  /** number of dipolar particles. */
  std::size_t sum_dip_part = 0;
  /** Sum of square of magnetic dipoles. */
  double sum_mu2 = 0.;

  /** position shift for calculation of first assignment mesh point. */
  double pos_shift = 0.;

  /** cached k-space self-energy correction */
  double energy_correction = 0.;
  /** k-space scalar mesh for k-space calculations. */
  std::vector<ComplexType> ks_scalar;

  p3m_interpolation_cache inter_weights;

  P3MFFTMesh<FloatType> mesh;

  /** @brief FFT algorithm. */
  std::unique_ptr<FFTBackend<FloatType>> fft;
  /** @brief FFT buffers. */
  std::unique_ptr<FFTBuffers<FloatType>> fft_buffers;

  void update_mesh_views() {
    auto const mesh_size_ptr = fft->get_mesh_size();
    auto const mesh_start_ptr = fft->get_mesh_start();
    for (auto i = 0u; i < 3u; ++i) {
      mesh.size[i] = mesh_size_ptr[i];
      mesh.start[i] = mesh_start_ptr[i];
    }
    mesh.stop = mesh.start + mesh.size;
    fft_buffers->update_mesh_views(mesh);
  }

  template <typename T, class... Args> void make_fft_instance(Args... args) {
    assert(fft == nullptr);
    fft = std::make_unique<T>(std::as_const(local_mesh), args...);
  }

  template <typename T, class... Args> void make_mesh_instance(Args... args) {
    assert(fft_buffers == nullptr);
    fft_buffers = std::make_unique<T>(std::as_const(local_mesh), args...);
  }

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  struct {
    /** dipole density in real-space with halo */
    std::array<std::vector<FloatType>, 3u> rs_dipole_density;
    /** dipole density in k-space without halo */
    std::array<std::vector<ComplexType>, 3u> ks_dipole_density;
    /** magnetic fields in real-space with halo */
    std::array<std::vector<FloatType>, 3> rs_B_fields;
    /** magnetic fields in k-space without halo */
    std::vector<ComplexType> ks_B_field_storage;
    /** magnetic fields in real-space without halo */
    std::array<std::vector<FloatType>, 3u> rs_B_fields_no_halo;
    /** k-space workspace without halo */
    std::vector<ComplexType> ks_scalar;
    /** @brief Force optimised influence function (k-space) */
    std::vector<FloatType> g_force;
    /** @brief Energy optimised influence function (k-space) */
    std::vector<FloatType> g_energy;
    p3m_send_mesh<FloatType> halo_comm;
    std::shared_ptr<P3MFFT<FloatType, FFTConfig>> fft;
    int world_size;
  } heffte;
  void resize_heffte_buffers();
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace>
      rs_fields_kokkos;

  void init_labels() {
    assert(ks_scalar.empty());
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
    assert(heffte.rs_dipole_density[0u].empty());
    assert(heffte.rs_dipole_density[1u].empty());
    assert(heffte.rs_dipole_density[2u].empty());
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
    rs_fields_kokkos = decltype(rs_fields_kokkos)(
        "DipolarP3MState::rs_fields_kokkos", 0, 0, 0);
  }
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
};

template <typename FloatType, Arch Architecture, class FFTConfig>
struct DipolarP3MHeffte : public DipolarP3M {
  ~DipolarP3MHeffte() override = default;

  // heFFTe checks are only implemented for row-major data
  static_assert(FFTConfig::r_space_order == Utils::MemoryOrder::ROW_MAJOR);
  static_assert(FFTConfig::k_space_order == Utils::MemoryOrder::ROW_MAJOR);

  using DipolarP3MStateClass = DipolarP3MState<FloatType, FFTConfig>;
  /** @brief Dipolar P3M parameters. */
  DipolarP3MStateClass &dp3m;

private:
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  // kokkos handle must outlive kokkos data structures from other class members
  std::shared_ptr<KokkosHandle> m_kokkos_handle;
#endif
  /** @brief Dipolar P3M meshes and FFT algorithm. */
  std::unique_ptr<DipolarP3MStateClass> dp3m_impl;
  TuningParameters tuning;
  bool m_is_tuned;

public:
  DipolarP3MHeffte(std::unique_ptr<DipolarP3MStateClass> &&dp3m_state,
                   TuningParameters tuning_params, double prefactor)
      : DipolarP3M(dp3m_state->params), dp3m{*dp3m_state},
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
        m_kokkos_handle{::kokkos_handle},
#endif
        dp3m_impl{std::move(dp3m_state)}, tuning{std::move(tuning_params)} {

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

  double long_range_energy() override { return long_range_kernel(false, true); }

  void add_long_range_forces() override {
    if constexpr (Architecture == Arch::CPU) {
      long_range_kernel(true, false);
    }
  }

  void dipole_assign() override;

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
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
    for (auto dir : {0u, 1u, 2u}) {
      dp3m.heffte.rs_dipole_density[dir].resize(dp3m.local_mesh.size);
      std::ranges::fill(dp3m.heffte.rs_dipole_density[dir], FloatType{0});
    }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  }

protected:
  /** Compute the k-space part of forces and energies. */
  double long_range_kernel(bool force_flag, bool energy_flag);
  double calc_average_self_energy_k_space() const override;
  void calc_energy_correction() override;
  void calc_influence_function_force() override;
  void calc_influence_function_energy() override;
  double calc_surface_term(bool force_flag, bool energy_flag) override;
  void init_cpu_kernels();
  void scaleby_box_l() override;
#ifdef ESPRESSO_NPT
  void npt_add_virial_contribution(double energy) const override;
#endif
};

#endif // ESPRESSO_DP3M
