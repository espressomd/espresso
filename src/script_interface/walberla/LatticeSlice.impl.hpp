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

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/scatter.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/multi_array.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ScriptInterface::walberla {

namespace detail {

// boundary types are always boost::optional types
template <class> struct is_optional_type : public std::false_type {};
template <class T>
struct is_optional_type<boost::optional<T>> : public std::true_type {};

template <class ArrayView, typename T>
void unflatten_grid(ArrayView &view, std::vector<T> const &values) {
  using array_type = boost::multi_array<T, 4>;
  auto it = values.begin();
  auto const dim_i = static_cast<typename array_type::index>(view.shape()[0]);
  auto const dim_j = static_cast<typename array_type::index>(view.shape()[1]);
  auto const dim_k = static_cast<typename array_type::index>(view.shape()[2]);
  auto const dim_t = static_cast<typename array_type::index>(view.shape()[3]);
  for (typename array_type::index i = 0; i != dim_i; ++i) {
    for (typename array_type::index j = 0; j != dim_j; ++j) {
      for (typename array_type::index k = 0; k != dim_k; ++k) {
        for (typename array_type::index t = 0; t != dim_t; ++t) {
          view[i][j][k][t] = *it;
          ++it;
        }
      }
    }
  }
}

template <class FieldSerializer, class ArrayView, typename T>
void flatten_grid(ArrayView const &view, std::vector<T> &out,
                  double units_conversion) {
  using array_type = boost::multi_array<T, 4>;
  out.reserve(view.num_elements());
  auto const dim_i = static_cast<typename array_type::index>(view.shape()[0]);
  auto const dim_j = static_cast<typename array_type::index>(view.shape()[1]);
  auto const dim_k = static_cast<typename array_type::index>(view.shape()[2]);
  auto const dim_t = static_cast<typename array_type::index>(view.shape()[3]);
  for (typename array_type::index i = 0; i != dim_i; ++i) {
    for (typename array_type::index j = 0; j != dim_j; ++j) {
      for (typename array_type::index k = 0; k != dim_k; ++k) {
        for (typename array_type::index t = 0; t != dim_t; ++t) {
          if constexpr (std::is_floating_point_v<T>) {
            out.emplace_back(view[i][j][k][t] * units_conversion);
          } else if constexpr (is_optional_type<T>{}) {
            if (view[i][j][k][t]) {
              out.emplace_back(*(view[i][j][k][t]) * units_conversion);
            } else {
              out.emplace_back(boost::none);
            }
          } else {
            out.emplace_back(view[i][j][k][t]);
          }
        }
      }
    }
  }
}

} // namespace detail

inline auto gather_slices_topology(boost::mpi::communicator const &comm,
                                   Utils::Vector3i const &local_lower_corner,
                                   Utils::Vector3i const &local_upper_corner) {
  std::vector<Utils::Vector3i> nodes_lower_corners;
  std::vector<Utils::Vector3i> nodes_upper_corners;
  boost::mpi::gather(comm, local_lower_corner, nodes_lower_corners, 0);
  boost::mpi::gather(comm, local_upper_corner, nodes_upper_corners, 0);
  return std::make_tuple(nodes_lower_corners, nodes_upper_corners);
}

template <class FieldSerializer>
template <class LatticeModel, typename T>
Variant LatticeSlice<FieldSerializer>::gather_3d(
    VariantMap const &params, std::vector<int> const &data_dims,
    LatticeModel const &lattice_model,
    std::vector<T> (LatticeModel::*getter)(Utils::Vector3i const &,
                                           Utils::Vector3i const &) const,
    double units_conversion) const {
  auto const &comm = context()->get_comm();
  auto const [slice_lower_corner, slice_upper_corner, local_lower_corner,
              local_upper_corner] = get_slices_bounding_boxes();
  auto const [nodes_lower_corners, nodes_upper_corners] =
      gather_slices_topology(comm, local_lower_corner, local_upper_corner);
  auto const data_size = std::accumulate(data_dims.cbegin(), data_dims.cend(),
                                         1, std::multiplies<>());
  auto const local_values =
      (lattice_model.*getter)(local_lower_corner, local_upper_corner);
  std::vector<std::vector<T>> nodes_values;
  boost::mpi::gather(comm, local_values, nodes_values, 0);
  if (comm.rank() == 0) {
    auto const dims = slice_upper_corner - slice_lower_corner;
    using index_range = boost::multi_array_types::index_range;
    using array_type = boost::multi_array<T, 4>;
    array_type array(boost::extents[dims[0]][dims[1]][dims[2]][data_size]);
    // populate the 3D array with data from each node
    for (std::size_t rank = 0; rank < nodes_values.size(); ++rank) {
      if (nodes_values[rank].empty()) {
        continue;
      }
      auto const range_lower_corner =
          nodes_lower_corners[rank] - slice_lower_corner;
      auto const range_upper_corner =
          nodes_upper_corners[rank] - slice_lower_corner;
      auto const local_range = [&](int j) {
        return index_range(range_lower_corner[j], range_upper_corner[j]);
      };
      typename array_type::template array_view<4>::type view =
          array[boost::indices[local_range(0)][local_range(1)][local_range(2)]
                              [index_range()]];
      detail::unflatten_grid(view, nodes_values[rank]);
    }
    // create the output flat array
    std::vector<T> out;
    detail::flatten_grid<FieldSerializer>(array, out, units_conversion);
    std::vector<int> shape = {m_shape.begin(), m_shape.end()};
    if (not(data_dims.size() == 1ul and data_dims[0] == 1)) {
      shape.insert(shape.end(), data_dims.begin(), data_dims.end());
    }
    auto const variant = FieldSerializer::serialize(out);
    return {std::vector<Variant>{{variant, Variant{shape}}}};
  }
  return {};
}

template <class FieldSerializer>
template <class LatticeModel, typename T>
void LatticeSlice<FieldSerializer>::scatter_3d(
    VariantMap const &params, std::vector<int> const &data_dims,
    LatticeModel &lattice_model,
    void (LatticeModel::*setter)(Utils::Vector3i const &,
                                 Utils::Vector3i const &,
                                 std::vector<T> const &),
    double units_conversion) {
  auto const &comm = context()->get_comm();
  auto const [slice_lower_corner, slice_upper_corner, local_lower_corner,
              local_upper_corner] = get_slices_bounding_boxes();
  auto const [nodes_lower_corners, nodes_upper_corners] =
      gather_slices_topology(comm, local_lower_corner, local_upper_corner);
  auto const data_size = std::accumulate(data_dims.cbegin(), data_dims.cend(),
                                         1, std::multiplies<>());
  auto const sentinel = get_sentinel_index(lattice_model.get_lattice());
  std::vector<std::vector<T>> nodes_values(comm.size());
  if (comm.rank() == 0) {
    auto const values =
        FieldSerializer::template deserialize<T>(params.at("values"));
    auto const dims = slice_upper_corner - slice_lower_corner;
    using index_range = boost::multi_array_types::index_range;
    using array_type = boost::multi_array<T, 4>;
    array_type array(boost::extents[dims[0]][dims[1]][dims[2]][data_size]);
    // populate the 3D array from the input flat array
    detail::unflatten_grid(array, values);
    // partition the 3D array into individual flat arrays for each MPI rank
    for (std::size_t rank = 0; rank < nodes_lower_corners.size(); ++rank) {
      auto const range_lower = nodes_lower_corners[rank] - slice_lower_corner;
      auto const range_upper = nodes_upper_corners[rank] - slice_lower_corner;
      if (not(range_lower > Utils::Vector3i::broadcast(sentinel))) {
        continue;
      }
      auto const local_range = [&](int j) {
        return index_range(range_lower[j], range_upper[j]);
      };
      typename array_type::template array_view<4>::type view =
          array[boost::indices[local_range(0)][local_range(1)][local_range(2)]
                              [index_range()]];
      detail::flatten_grid<FieldSerializer>(view, nodes_values[rank],
                                            units_conversion);
    }
  }
  std::vector<T> local_values;
  boost::mpi::scatter(comm, nodes_values, local_values, 0);
  (lattice_model.*setter)(local_lower_corner, local_upper_corner, local_values);
  lattice_model.ghost_communication();
}

} // namespace ScriptInterface::walberla
