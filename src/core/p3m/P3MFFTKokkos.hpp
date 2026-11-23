/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#ifdef ESPRESSO_KOKKOS_FFT

#include "p3m/P3MFFTBackend.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <KokkosFFT.hpp>
#include <Kokkos_Core.hpp>

#include <cassert>
#include <complex>
#include <memory>
#include <vector>

/**
 * @brief Single-MPI-rank P3M FFT backend built on kokkos-fft.
 *
 * kokkos-fft wraps the host FFT library (FFTW) behind a @c Kokkos::View API
 * and performs purely local transforms, so it is only selected when the
 * simulation runs on a single MPI rank; multi-rank runs keep the heFFTe
 * backend (@ref P3MFFTHeffte). It reproduces heFFTe's row-major r2c layout
 * (reduced last axis, size @c mesh[2]/2+1) and its unscaled convention
 * (@c scale::none in both directions; P3M folds the @c 1/N into the influence
 * function), so forces and energies agree with the heFFTe path to
 * floating-point round-off. The transform runs on the host regardless of the
 * Kokkos default device, matching the CPU heFFTe backend it replaces.
 *
 * The transforms execute in place on the caller's buffers via kokkos-fft's
 * new-array execute path (FFTW @c fftw_execute_dft_r2c / @c _c2r on the pointer
 * pair passed at call time). Because kokkos-fft plans with @c FFTW_ESTIMATE and
 * without @c FFTW_UNALIGNED, a plan is only reused on buffers whose alignment
 * matches the ones it was built with. We guarantee that by building each plan
 * from the exact buffers it will run on (keyed and cached by pointer): the P3M
 * k-space and no-halo real-space buffers are allocated once and stable, so
 * their plans are built on the first step and reused thereafter with no data
 * copies. The forward input goes through an owned, consistently aligned
 * scratch view for the same reason; P3M fills it directly (via
 * @ref forward_input_buffer) so no staging copy occurs, and only an external
 * caller passing a foreign buffer pays one.
 */
template <typename FloatType, class FFTConfig>
struct P3MFFTKokkos final : public P3MFFTBackend<FloatType, FFTConfig> {
  static_assert(FFTConfig::use_r2c,
                "kokkos-fft P3M backend implements the r2c transform only");
  static_assert(FFTConfig::r2c_dir == 2u,
                "kokkos-fft P3M backend reduces the contiguous last axis");
  static_assert(FFTConfig::r_space_order == Utils::MemoryOrder::ROW_MAJOR and
                    FFTConfig::k_space_order == Utils::MemoryOrder::ROW_MAJOR,
                "kokkos-fft P3M backend supports row-major layout only");

  using Base = P3MFFTBackend<FloatType, FFTConfig>;
  using ComplexType = typename Base::ComplexType;
  using RSpaceScalar = typename Base::RSpaceScalar;
  using KComplex = Kokkos::complex<FloatType>;
  using ExecSpace = Kokkos::DefaultHostExecutionSpace;

  using RealView =
      Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace>;
  using RealViewU =
      Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using CplxViewU =
      Kokkos::View<KComplex ***, Kokkos::LayoutRight, Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ForwardPlan = KokkosFFT::Plan<ExecSpace, RealViewU, CplxViewU, 3>;
  using BackwardPlan = KokkosFFT::Plan<ExecSpace, CplxViewU, RealViewU, 3>;

  P3MFFTKokkos(boost::mpi::communicator comm,
               Utils::Vector3i const &global_mesh,
               Utils::Vector3i const &rs_local_ld_index,
               Utils::Vector3i const &rs_local_ur_index,
               Utils::Vector3i const & /* node_grid */)
      : m_mesh{global_mesh},
        m_ks_size{global_mesh[0], global_mesh[1], global_mesh[2] / 2 + 1} {
    // kokkos-fft is local-only; single rank means the no-halo local mesh spans
    // the whole global mesh, so there is no domain decomposition to mirror.
    assert(comm.size() == 1);
    assert((rs_local_ur_index - rs_local_ld_index) == global_mesh);
    static_cast<void>(comm);
    static_cast<void>(rs_local_ld_index);
    static_cast<void>(rs_local_ur_index);

    m_real_scratch = RealView(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                 "P3MFFTKokkos::real_scratch"),
                              m_mesh[0], m_mesh[1], m_mesh[2]);
  }

  Utils::Vector3i ks_local_ld_index() const override { return {0, 0, 0}; }
  Utils::Vector3i ks_local_ur_index() const override { return m_ks_size; }
  Utils::Vector3i ks_local_size() const override { return m_ks_size; }
  Utils::Vector3i rs_local_size() const override { return m_mesh; }

  // Hand the caller our owned, aligned scratch so it can fill the input
  // directly (via extract_block); then forward() runs in place with no copy.
  RSpaceScalar *forward_input_buffer() override {
    return m_real_scratch.data();
  }

  void forward(RSpaceScalar const *in, ComplexType *out) override {
    RealViewU const scratch_view(m_real_scratch.data(), m_mesh[0], m_mesh[1],
                                 m_mesh[2]);
    // If the caller filled our scratch in place (via forward_input_buffer),
    // there is nothing to copy; otherwise stage the external input.
    if (in != m_real_scratch.data()) {
      RealViewU const in_view(const_cast<FloatType *>(in), m_mesh[0], m_mesh[1],
                              m_mesh[2]);
      Kokkos::deep_copy(m_real_scratch, in_view);
    }
    CplxViewU const out_view(reinterpret_cast<KComplex *>(out), m_ks_size[0],
                             m_ks_size[1], m_ks_size[2]);
    if (not m_forward or m_forward_out != out) {
      m_forward = std::make_unique<ForwardPlan>(
          ExecSpace{}, scratch_view, out_view, KokkosFFT::Direction::forward,
          KokkosFFT::axis_type<3>({0, 1, 2}));
      m_forward_out = out;
    }
    KokkosFFT::execute(*m_forward, scratch_view, out_view,
                       KokkosFFT::Normalization::none);
  }

  void backward(ComplexType *in, RSpaceScalar *out) override {
    // The c2r transform runs in place on the input and destroys it, which
    // the interface's non-const parameter permits (see
    // P3MFFTBackend::backward).
    CplxViewU const in_view(reinterpret_cast<KComplex *>(in), m_ks_size[0],
                            m_ks_size[1], m_ks_size[2]);
    RealViewU const out_view(out, m_mesh[0], m_mesh[1], m_mesh[2]);
    KokkosFFT::execute(backward_plan(in, out, in_view, out_view), in_view,
                       out_view, KokkosFFT::Normalization::none);
  }

private:
  /** @brief Plan bound to a fixed (input, output) buffer pair. */
  struct BackwardEntry {
    void const *in;
    void const *out;
    std::unique_ptr<BackwardPlan> plan;
  };

  /**
   * @brief Fetch (or lazily build) the backward plan for a buffer pair.
   *
   * P3M calls @ref backward with the same handful of stable (@c ks_E_fields,
   * @c rs_E_fields_no_halo) buffers every step, so building the plan on the
   * exact buffers it runs on both fixes the FFTW alignment and lets it be
   * reused indefinitely without a copy.
   */
  BackwardPlan &backward_plan(void const *in, void const *out,
                              CplxViewU const &in_view,
                              RealViewU const &out_view) {
    for (auto &entry : m_backward) {
      if (entry.in == in and entry.out == out) {
        return *entry.plan;
      }
    }
    auto plan = std::make_unique<BackwardPlan>(
        ExecSpace{}, in_view, out_view, KokkosFFT::Direction::backward,
        KokkosFFT::axis_type<3>({0, 1, 2}));
    m_backward.push_back({in, out, std::move(plan)});
    return *m_backward.back().plan;
  }

  Utils::Vector3i m_mesh;
  Utils::Vector3i m_ks_size;
  RealView m_real_scratch;
  ComplexType const *m_forward_out = nullptr;
  std::unique_ptr<ForwardPlan> m_forward;
  std::vector<BackwardEntry> m_backward;
};

#endif // ESPRESSO_KOKKOS_FFT
