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
#include "p3m/common.hpp"

#include "communication.hpp"
#include "fft/vector.hpp"

#include <instrumentation/fe_trap.hpp>

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/operations.hpp>

#include <heffte.h>
#include <heffte_backends.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

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
  using buffer_container = std::conditional_t<
      Architecture == Arch::CPU, fft::vector<T>,
      typename heffte::fft3d<backend>::template buffer_container<T>>;

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
  /* batched transform size */
  unsigned m_batch_size = 0u;
  /* whether a batched transform is both possible and faster */
  bool use_batched_fft = false;

  template <typename T, std::size_t N>
  static auto to_array(Utils::Vector<T, N> const &vec) {
    std::array<T, N> res{};
    std::ranges::copy(vec, res.begin());
    return res;
  }

  /**
   * @brief Time the FFTs of one P3M time step on scratch buffers.
   *
   * One forward transform followed by @p n_fields backward transforms,
   * either as separate calls or as one batched call. The first step is
   * not timed: FFTW creates its plans on first use. This is a collective
   * operation that returns the mean time per step of the slowest rank,
   * hence all ranks get the same value.
   */
  static double time_fft_step(FFT3D const &fft,
                              boost::mpi::communicator const &comm, unsigned n,
                              bool batched, unsigned n_steps) {
    using RSpaceType =
        std::conditional_t<FFTConfig::use_r2c, FloatType, OutputType>;
    auto const rs_size = static_cast<std::size_t>(fft.size_inbox());
    auto const ks_size = static_cast<std::size_t>(fft.size_outbox());
    auto const ws_size = static_cast<std::size_t>(fft.size_workspace());
    // scratch buffers with batched size
    buffer_container<RSpaceType> rs_buffer(n * rs_size);
    buffer_container<OutputType> ks_buffer(n * ks_size);
    // scratch workspace with normal buffer size if no batching was requested
    buffer_container<OutputType> workspace((batched ? n : 1u) * ws_size);
    auto const step = [&]() {
      fft.forward(rs_buffer.data(), ks_buffer.data(), workspace.data());
      if (batched) {
        fft.backward(n, ks_buffer.data(), rs_buffer.data(), workspace.data());
      } else {
        for (std::size_t i = 0u; i < n; ++i) {
          fft.backward(ks_buffer.data() + i * ks_size,
                       rs_buffer.data() + i * rs_size, workspace.data());
        }
      }
    };
    step(); // warm-up: FFTW plans are created on first use
    comm.barrier();
    if (n_steps == 0u) {
      return 0.;
    }
    auto const tick = MPI_Wtime();
    for (auto i = 0u; i < n_steps; ++i) {
      step();
    }
    auto const tock = MPI_Wtime();
    auto const time = (tock - tick) / static_cast<double>(n_steps);
    return boost::mpi::all_reduce(comm, time, boost::mpi::maximum<double>());
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
         Utils::Vector3i const &node_grid, unsigned n_fields = 1u) {
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
    m_batch_size = n_fields;
    auto const this_rank = comm.rank();
    // initial output boxes based on prime factors; on CPU with several ranks
    // they only seed the plan and are replaced by its last pencil/slab stage
    auto out_boxes = heffte::split_world(global_out_box, to_array(best_grid));
    out_box = std::make_unique<Box>(out_boxes[this_rank]);

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
    // 1-D pencils; for sufficiently large problem, it is expected that the
    // pencil decomposition is better but for smaller problems, the slabs may
    // perform better (depending on hardware and backend)
    options.use_pencils = true;

    auto const make_fft3d = [&](Box const &output_box,
                                heffte::plan_options const &opts) {
      if constexpr (FFTConfig::use_r2c) {
        auto constexpr r2c_dir = static_cast<int>(FFTConfig::r2c_dir);
        return std::make_unique<FFT3D>(gpu_stream, *in_box, output_box, r2c_dir,
                                       comm, opts);
      } else {
        return std::make_unique<FFT3D>(gpu_stream, *in_box, output_box, comm,
                                       opts);
      }
    };
    auto const make_plan = [&](heffte::ioboxes<int> const &boxes,
                               bool use_pencils) {
      auto opts = options;
      opts.use_pencils = use_pencils;
      auto constexpr r2c_dir = static_cast<int>(FFTConfig::r2c_dir);
      return heffte::plan_operations(
          boxes, FFTConfig::use_r2c ? r2c_dir : -1,
          heffte::set_options<backend, FFTConfig::use_r2c>(opts), this_rank);
    };
    struct Candidate { // LCOV_EXCL_LINE
      std::vector<Box> out_boxes;
      unsigned int n_mpi_reshapes;
      bool valid;
    };
    auto const is_empty = [](Box const &box) { return box.empty(); };
    auto const make_candidate = [&](heffte::ioboxes<int> const &seed,
                                    bool use_pencils) {
      auto const seed_plan = make_plan(seed, use_pencils);
      auto boxes = heffte::ioboxes<int>{seed.in, {}};
      for (auto const &box : seed_plan.out_shape[2]) {
        boxes.out.emplace_back(box.low, box.high, out_box_order);
      }
      auto const plan = make_plan(boxes, use_pencils);
      auto n_mpi_reshapes = 0u;
      auto valid = true;
      for (auto stage : {0u, 1u, 2u, 3u}) {
        if (not heffte::match(plan.in_shape[stage], plan.out_shape[stage])) {
          ++n_mpi_reshapes;
        }
        valid &= std::ranges::none_of(plan.out_shape[stage], is_empty);
      }
      return Candidate{std::move(boxes.out), n_mpi_reshapes, valid};
    };
    // use the last FFT stage of a heFFTe plan for the input boxes as output
    // boxes: the final reshape of each forward transform (and the first one
    // of each backward transform) then becomes a local operation. This costs
    // extra planning steps and one MPI_Allgather of the boxes per FFT
    // construction. Pencils are the default; slabs need fewer MPI reshapes
    // unless the input boxes already are slabs, but each slab reshape is an
    // all-to-all with many more messages, so slabs are only kept when they
    // are measured to be faster.
    // The backward transforms of several fields can be batched, resulting in a
    // single MPI exchange per MPI reshape instead of @c n_fields exchanges,
    // but messages are @c n_fields times larger. Whether this pays off
    // depends on the MPI transport (e.g. larger messages can switch from the
    // eager to the rendezvous protocol), so this is measured as well, jointly
    // with the decomposition.
    if constexpr (Architecture == Arch::CPU) {
      auto const seed = heffte::mpi::gather_boxes(*in_box, *out_box, comm);
      auto const pencils = make_candidate(seed, true);
      auto const slabs = make_candidate(seed, false);
      struct Layout {
        std::unique_ptr<FFT3D> fft;
        Box box;
        bool use_pencils;
        bool batchable;
        bool simd_veto;
      };
      auto const make_layout = [&](Box const &obox, bool pencils) -> Layout {
        auto opts = options;
        opts.use_pencils = pencils;
        auto fft = make_fft3d(obox, opts);
        auto const is_local_buffer_even = obox.count() % 2ll == 0ll;
        auto const batchable =
            (n_fields > 1u) and
            boost::mpi::all_reduce(comm, is_local_buffer_even,
                                   std::logical_and<bool>());
        auto const local_simd_veto =
            std::is_same_v<FloatType, float> and
            (not is_local_buffer_even or fft->size_comm_buffers() % 2u == 1u);
        auto const simd_veto = boost::mpi::all_reduce(comm, local_simd_veto,
                                                      std::logical_or<bool>());
        return {std::move(fft), obox, pencils, batchable, simd_veto};
      };
      std::vector<Layout> layouts;
      layouts.reserve(2u);
      layouts.emplace_back(make_layout(pencils.out_boxes[this_rank], true));
      // slab decomposition only works with column-major data in C2C 2D FFT
      // https://github.com/icl-utk-edu/heffte/issues/93
      if ((FFTConfig::use_r2c or
           (FFTConfig::r_space_order == Utils::MemoryOrder::COLUMN_MAJOR and
            FFTConfig::k_space_order == Utils::MemoryOrder::COLUMN_MAJOR)) and
          slabs.valid and slabs.n_mpi_reshapes < pencils.n_mpi_reshapes) {
        layouts.emplace_back(make_layout(slabs.out_boxes[this_rank], false));
      }
      auto *best_layout = &layouts.front();
      auto best_batched = false;
      if (layouts.size() >= 2ul or layouts[0].batchable) {
        auto constexpr n_rounds = 3u;
        auto constexpr n_steps = 5u;
        auto best_time = std::numeric_limits<double>::max();
        for (auto &layout : layouts) {
          if (layout.simd_veto) {
            // skip candidates with memory alignment issues
            continue;
          }
          auto const batch_modes = layout.batchable
                                       ? std::vector<bool>{false, true}
                                       : std::vector<bool>{false};
          for (bool const batched : batch_modes) {
            auto min_time = std::numeric_limits<double>::max();
            for (auto round = 0u; round < n_rounds; ++round) {
              auto cur_time =
                  time_fft_step(*layout.fft, comm, n_fields, batched, n_steps);
              if (batched) {
                // batched FFT runtimes are skewed by the MPI barrier;
                // account for MPI load imbalance with a small penalty
                cur_time *= (comm.size() > 4) ? 1.25 : 1.1;
              }
              min_time = std::min(min_time, cur_time);
            }
            if (min_time < best_time) {
              best_time = min_time;
              best_layout = &layout;
              best_batched = batched;
            }
          }
        }
      }
      if (best_layout->simd_veto) {
        throw p3m_tuning_exception("cannot align memory for SIMD");
      }
      options.use_pencils = best_layout->use_pencils;
      out_box = std::make_unique<Box>(best_layout->box);
      fft3d = std::move(best_layout->fft);
      use_batched_fft = best_batched;
    }
    if constexpr (Architecture == Arch::CUDA) {
      options.use_gpu_aware = ::communication_environment->is_mpi_gpu_aware();
    }
#ifdef ESPRESSO_FPE
    // cuFFT builds device kernels using CUDA-JIT
    // (https://docs.nvidia.com/cuda/archive/13.1.1/cufft/#plan-initialization-time)
    // but this operation is not guaranteed to succeed for all mesh sizes,
    // and in rare cases, it can send the SIGFPE signal
    auto const trap_pause = (Architecture == Arch::CUDA)
                                ? fe_trap::make_shared_pause_scoped()
                                : nullptr;
#endif

    if (not fft3d) {
      fft3d = make_fft3d(*out_box, options);
      // warm-up: FFTW plans are created on first use
      time_fft_step(*fft3d, comm, 1, false, 0u);
    }
    m_workspace = decltype(m_workspace)((use_batched_fft ? n_fields : 1u) *
                                        fft3d->size_workspace());
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
  /**
   * Backward transform of @p n fields stored back to back, batched
   * or one by one, whichever was measured to be faster.
   */
  void backward_batch(unsigned const n, auto *const in, auto *const out) {
    if (use_batched_fft) {
      assert(n <= m_batch_size);
      fft3d->backward(static_cast<int>(n), in, out, m_workspace.data());
    } else {
      auto const ks_size = fft3d->size_outbox();
      auto const rs_size = fft3d->size_inbox();
      for (auto i = 0u; i < n; ++i) {
        fft3d->backward(in + i * ks_size, out + i * rs_size,
                        m_workspace.data());
      }
    }
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
               Utils::Vector3i const &node_grid, unsigned n_fields = 1)
      : m_impl(nullptr, comm, global_mesh, rs_local_ld_index, rs_local_ur_index,
               node_grid, n_fields),
        m_input(static_cast<std::size_t>(Utils::product(rs_local_size()))) {}

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
  void backward_batch(unsigned const n, ComplexType *const in,
                      RSpaceScalar *const out) override {
    m_impl.backward_batch(n, in, out);
  }

private:
  P3MFFT<FloatType, Arch::CPU, FFTConfig> m_impl;
  fft::vector<RSpaceScalar> m_input;
};
