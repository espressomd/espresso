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

#include "p3m/P3MFFTBackend.hpp"

#include "communication.hpp"

#include <instrumentation/fe_trap.hpp>

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <boost/mpi/communicator.hpp>

#include <heffte.h>
#include <heffte_backends.h>

#include <algorithm>
#include <array>
#include <initializer_list>
#include <memory>

#if defined(__CUDACC__)
#include "cuda/utils.cuh"
#endif

/**
 * @brief FFT manager.
 */
template <typename FloatType, Arch Architecture, class FFTConfig> class P3MFFT {
public:
  using OutputType = typename heffte::fft_output<FloatType>::type;
  using backend =
      std::conditional_t<Architecture == Arch::CPU, heffte::backend::fftw,
                         heffte::backend::cufft>;
  template <class T = OutputType>
  using buffer_container = heffte::fft3d<backend>::template buffer_container<T>;

private:
  using FFT3D =
      std::conditional_t<FFTConfig::use_r2c, heffte::fft3d_r2c<backend>,
                         heffte::fft3d<backend>>;
  using Box = heffte::box3d<>;
  using stream_type =
      heffte::backend::device_instance<heffte::tag::gpu>::stream_type;

  /* input box */
  std::unique_ptr<Box> in_box;
  /* output box */
  std::unique_ptr<Box> out_box;
  /* workspace for the FFT */
  buffer_container<OutputType> m_workspace;
  /* FFT backend */
  std::unique_ptr<FFT3D> fft3d;
  std::shared_ptr<boost::mpi::environment> m_mpi_env_lock;

  template <typename T, std::size_t N>
  static auto to_array(Utils::Vector<T, N> const &vec) {
    std::array<T, N> res{};
    std::ranges::copy(vec, res.begin());
    return res;
  }

public:
  ~P3MFFT() {
    fft3d.reset();
    m_mpi_env_lock.reset();
  }
  P3MFFT(stream_type gpu_stream, boost::mpi::communicator comm,
         Utils::Vector3i const &global_mesh,
         Utils::Vector3i const &rs_local_ld_index,
         Utils::Vector3i const &rs_local_ur_index,
         Utils::Vector3i const &node_grid) {
    auto constexpr row_major_order = std::array<int, 3>{2, 1, 0};
    auto constexpr col_major_order = std::array<int, 3>{0, 1, 2};
    auto constexpr in_box_order =
        (FFTConfig::r_space_order == Utils::MemoryOrder::ROW_MAJOR)
            ? row_major_order
            : col_major_order;
    auto constexpr out_box_order =
        (FFTConfig::k_space_order == Utils::MemoryOrder::ROW_MAJOR)
            ? row_major_order
            : col_major_order;
    auto const n_procs = Utils::product(node_grid);
    auto const high = to_array(global_mesh - Utils::Vector3i::broadcast(1));
    auto const global_out_box_full = Box({0, 0, 0}, high, out_box_order);
    auto const global_out_box =
        FFTConfig::use_r2c ? global_out_box_full.r2c(FFTConfig::r2c_dir)
                           : global_out_box_full;
    auto best_grid = node_grid;
    for (auto i : {0u, 1u, 2u}) {
      if (global_mesh[i] % (2 * n_procs) == 0) {
        best_grid = {n_procs, 1, 1};
        break;
      }
    }
    // use optimal output box decomposition based on prime factors
    auto out_boxes = heffte::split_world(global_out_box, to_array(best_grid));
    out_box = std::make_unique<Box>(out_boxes[comm.rank()]);

    in_box = std::make_unique<Box>(
        to_array(rs_local_ld_index),
        to_array(rs_local_ur_index - Utils::Vector3i::broadcast(1)),
        in_box_order);

    // at this stage we can manually adjust some HeFFTe options
    heffte::plan_options options = heffte::default_options<backend>();

    // use strided 1-D FFT operations
    // some backends work just as well when the entries of the data are not
    // contiguous then there is no need to reorder the data in the intermediate
    // stages which saves time
    options.use_reorder = true;

    // use point-to-point communications
    // collaborative all-to-all and individual point-to-point communications are
    // two alternatives one may be better than the other depending on the
    // version of MPI, the hardware interconnect, and the problem size
    options.algorithm = heffte::reshape_algorithm::p2p_plined;

    // in the intermediate steps, the data can be shapes as either 2-D slabs or
    // 1-D pencils for sufficiently large problem, it is expected that the
    // pencil decomposition is better but for smaller problems, the slabs may
    // perform better (depending on hardware and backend)
    options.use_pencils = true;
#if defined(__CUDACC__)
    if constexpr (Architecture == Arch::CUDA) {
      options.use_gpu_aware = ::communication_environment->is_mpi_gpu_aware();
    }
#endif
#ifdef ESPRESSO_FPE
    // cuFFT builds device kernels using CUDA-JIT
    // (https://docs.nvidia.com/cuda/archive/13.1.1/cufft/#plan-initialization-time)
    // but this operation is not guaranteed to succeed for all mesh sizes,
    // and in rare cases, it can send the SIGFPE signal
    auto const trap_pause = (Architecture == Arch::CUDA)
                                ? fe_trap::make_shared_pause_scoped()
                                : nullptr;
#endif

    if constexpr (FFTConfig::use_r2c) {
      fft3d = std::make_unique<FFT3D>(gpu_stream, *in_box, *out_box,
                                      FFTConfig::r2c_dir, comm, options);
    } else {
      fft3d =
          std::make_unique<FFT3D>(gpu_stream, *in_box, *out_box, comm, options);
    }
    m_workspace = decltype(m_workspace)(fft3d->size_workspace());
    // MPI communicator is needed to destroy the FFT plans
    m_mpi_env_lock = ::communication_environment->get_mpi_env();
  }

  Utils::Vector3i ks_local_ld_index() const {
    return Utils::Vector3i(out_box->low);
  }
  Utils::Vector3i ks_local_ur_index() const {
    return Utils::Vector3i(out_box->high) + Utils::Vector3i::broadcast(1);
  }
  Utils::Vector3i ks_local_size() const {
    return ks_local_ur_index() - ks_local_ld_index();
  }
  Utils::Vector3i rs_local_size() const {
    return Utils::Vector3i(in_box->high) + Utils::Vector3i::broadcast(1) -
           Utils::Vector3i(in_box->low);
  }
  void forward(auto &&in, auto &&out) {
    fft3d->forward(in, out, m_workspace.data());
  }
  void backward(auto &&in, auto &&out) {
    fft3d->backward(in, out, m_workspace.data());
  }
};

/**
 * @brief heFFTe-backed implementation of the P3M FFT interface.
 *
 * Thin adapter that exposes @ref P3MFFT through the @ref P3MFFTBackend
 * virtual interface, so the solver can hold it interchangeably with other
 * backends. This is the general backend, used for any rank count.
 */
template <typename FloatType, class FFTConfig>
struct P3MFFTHeffte final : public P3MFFTBackend<FloatType, FFTConfig> {
  using Base = P3MFFTBackend<FloatType, FFTConfig>;
  using ComplexType = typename Base::ComplexType;
  using RSpaceScalar = typename Base::RSpaceScalar;

  P3MFFTHeffte(boost::mpi::communicator comm,
               Utils::Vector3i const &global_mesh,
               Utils::Vector3i const &rs_local_ld_index,
               Utils::Vector3i const &rs_local_ur_index,
               Utils::Vector3i const &node_grid)
      : m_impl(nullptr, comm, global_mesh, rs_local_ld_index, rs_local_ur_index,
               node_grid),
        m_input(
            static_cast<std::size_t>(Utils::product(m_impl.rs_local_size()))) {}

  Utils::Vector3i ks_local_ld_index() const override {
    return m_impl.ks_local_ld_index();
  }
  Utils::Vector3i ks_local_ur_index() const override {
    return m_impl.ks_local_ur_index();
  }
  Utils::Vector3i ks_local_size() const override {
    return m_impl.ks_local_size();
  }
  Utils::Vector3i rs_local_size() const override {
    return m_impl.rs_local_size();
  }
  RSpaceScalar *forward_input_buffer() override { return m_input.data(); }
  void forward(RSpaceScalar const *in, ComplexType *out) override {
    m_impl.forward(in, out);
  }
  void backward(ComplexType *in, RSpaceScalar *out) override {
    // heFFTe preserves the input; the non-const parameter is the interface's
    // contract for backends that cannot (see P3MFFTBackend::backward).
    m_impl.backward(in, out);
  }

private:
  P3MFFT<FloatType, Arch::CPU, FFTConfig> m_impl;
  std::vector<RSpaceScalar> m_input;
};
