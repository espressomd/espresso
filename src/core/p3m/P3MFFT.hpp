/*
 * Copyright (C) 2024-2025 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <heffte.h>
#include <heffte_backends.h>

#include <algorithm>
#include <array>
#include <initializer_list>
#include <memory>

template <typename T, std::size_t N>
auto to_array(Utils::Vector<T, N> const &vec) {
  std::array<T, N> res{};
  std::copy(vec.begin(), vec.end(), res.begin());
  return res;
};

template <typename FloatType> class P3MFFT {
private:
  using backend_tag = heffte::backend::default_backend<heffte::tag::cpu>::type;
  using Box = heffte::box3d<>;
  boost::mpi::communicator comm;
  Utils::Vector3i m_memory_layout;
  Utils::Vector3i m_global_mesh;
  std::shared_ptr<Box> in_box;
  std::shared_ptr<Box> out_box;
  std::vector<std::complex<FloatType>> m_buffer;
  heffte::fft3d<backend_tag> fft3d;

public:
  P3MFFT(boost::mpi::communicator comm, Utils::Vector3i const &global_mesh,
         Utils::Vector3i const &rs_local_ld_index,
         Utils::Vector3i const &rs_local_ur_index,
         Utils::Vector3i const &memory_layout)
      : comm(comm), m_memory_layout(memory_layout), m_global_mesh(global_mesh),
        in_box(std::make_shared<Box>(
            to_array(rs_local_ld_index),
            to_array(rs_local_ur_index - Utils::Vector3i::broadcast(1)),
            to_array(m_memory_layout))),
        out_box(std::make_shared<Box>(
            to_array(rs_local_ld_index),
            to_array(rs_local_ur_index - Utils::Vector3i::broadcast(1)),
            to_array(m_memory_layout))),
        fft3d(*in_box, *out_box, comm) {
    init_fft();
  }

  void set_preferred_kspace_decomposition(Utils::Vector3i const &node_grid) {
    auto const global_box = heffte::box3d<>(
        {0, 0, 0}, to_array(m_global_mesh - Utils::Vector3i::broadcast(1)),
        to_array(m_memory_layout));
    auto const n_procs = Utils::product(node_grid);
    auto best_grid = node_grid;
    for (auto i : {0u, 1u, 2u}) {
      if (m_global_mesh[i] % (2 * n_procs) == 0) {
        best_grid = {n_procs, 1, 1};
        break;
      }
    }
    auto all_boxes = heffte::split_world(
        global_box, {best_grid[0], best_grid[1], best_grid[2]});
    out_box = std::make_shared<Box>(all_boxes[comm.rank()]);
    init_fft();
  }

  void init_fft() {
    // at this stage we can manually adjust some HeFFTe options
    heffte::plan_options options = heffte::default_options<backend_tag>();

    // use strided 1-D FFT operations
    // some backends work just as well when the entries of the data are not
    // contiguous then there is no need to reorder the data in the intermediate
    // stages which saves time
    options.use_reorder = false;

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
    fft3d = heffte::fft3d<backend_tag>(*in_box, *out_box, comm, options);
    m_buffer.resize(fft3d.size_workspace());
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
  template <typename In, typename Out> void forward(In in, Out out) {
    fft3d.forward(in, out, m_buffer.data());
  }
  template <typename T1, typename T2> auto backward(T1 &in, T2 &out) {
    return fft3d.backward(in, out, m_buffer.data());
  }
  auto const &get_memory_layout() const { return m_memory_layout; }
};
