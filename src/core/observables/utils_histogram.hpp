/*
 * Copyright (C) 2023 The ESPResSo project
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

#include <utils/Histogram.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace Observables::detail {

/** @brief Gather data from all MPI ranks. */
template <class Pos>
auto gather(boost::mpi::communicator const &comm,
            std::vector<Pos> const &local_pos) {
  std::vector<std::vector<Pos>> global_pos{};
  global_pos.reserve(static_cast<std::size_t>(comm.size()));
  boost::mpi::gather(comm, local_pos, global_pos, 0);
  return global_pos;
}

/** @brief Gather data from all MPI ranks. */
template <class Pos, class Val>
auto gather(boost::mpi::communicator const &comm,
            std::vector<Pos> const &local_pos,
            std::vector<Val> const &local_val) {
  std::vector<std::vector<Pos>> global_pos{};
  global_pos.reserve(static_cast<std::size_t>(comm.size()));
  boost::mpi::gather(comm, local_pos, global_pos, 0);
  std::vector<std::vector<Val>> global_val{};
  global_val.reserve(static_cast<std::size_t>(comm.size()));
  boost::mpi::gather(comm, local_val, global_val, 0);
  return std::make_pair(global_pos, global_val);
}

/** @brief Accumulate histogram data gathered from multiple MPI ranks. */
template <class T, std::size_t N, std::size_t M, class U, class Pos, class Val>
void accumulate(Utils::Histogram<T, N, M, U> &histogram,
                std::vector<std::vector<Pos>> const &pos,
                std::vector<std::vector<Val>> const &val) {
  for (std::size_t rank = 0u; rank < pos.size(); ++rank) {
    auto const &pos_vec = pos[rank];
    auto const &val_vec = val[rank];
    for (std::size_t i = 0u; i < pos_vec.size(); ++i) {
      histogram.update(pos_vec[i], val_vec[i]);
    }
  }
}

struct empty_bin_exception {};

/** @brief Normalize histogram by the number of values in each bin. */
template <class T, std::size_t N, std::size_t M, class U>
auto normalize_by_bin_size(Utils::Histogram<T, N, M, U> &histogram,
                           bool allow_empty_bins = true) {
  auto hist_data = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  for (std::size_t i = 0u; i < hist_data.size(); ++i) {
    if (tot_count[i] != 0u) {
      hist_data[i] /= static_cast<double>(tot_count[i]);
    } else if (not allow_empty_bins) {
      throw empty_bin_exception{};
    }
  }
  return hist_data;
}

} // Namespace Observables::detail
