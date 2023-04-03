/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#ifdef WALBERLA

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include "FluidNodeWalberla.hpp"
#include "FluidWalberla.hpp"
#include "LatticeIndices.hpp"

#include "core/errorhandling.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/scatter.hpp>
#include <boost/multi_array.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

namespace ScriptInterface::walberla {

namespace detail {

template <typename T, class ArrayView>
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

template <typename T, class ArrayView>
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
          } else if constexpr (std::is_same_v<
                                   T, boost::optional<Utils::Vector3d>>) {
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

class LBSlice : public AutoParameters<LBSlice, LatticeIndices> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  std::shared_ptr<FluidWalberla> m_lb_sip;
  Utils::Vector3i m_slice_lower_corner;
  Utils::Vector3i m_slice_upper_corner;
  std::vector<int> m_shape;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;
  std::unordered_map<std::string, std::vector<int>> m_shape_val;

public:
  LBSlice() { add_parameters({}); }

  void do_construct(VariantMap const &params) override {
    try {
      m_lb_sip =
          get_value<std::shared_ptr<FluidWalberla>>(params, "parent_sip");
      m_lb_fluid = m_lb_sip->lb_fluid().lock();
      assert(m_lb_fluid);
      auto const &lb_params = m_lb_sip->lb_params().lock();
      assert(lb_params);
      auto const tau = lb_params->get_tau();
      auto const agrid = lb_params->get_agrid();
      m_conv_dens = Utils::int_pow<3>(agrid);
      m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
      m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
      m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    } catch (std::exception const &e) {
      if (context()->is_head_node()) {
        runtimeErrorMsg() << "LBSlice failed: " << e.what();
      }
      m_lb_fluid.reset();
      return;
    }
    m_shape = get_value<std::vector<int>>(params, "shape");
    m_slice_lower_corner =
        get_value<Utils::Vector3i>(params, "slice_lower_corner");
    m_slice_upper_corner =
        get_value<Utils::Vector3i>(params, "slice_upper_corner");
    m_shape_val["density"] = {1};
    m_shape_val["population"] = {static_cast<int>(m_lb_fluid->stencil_size())};
    m_shape_val["velocity"] = {3};
    m_shape_val["velocity_at_boundary"] = {1};
    m_shape_val["is_boundary"] = {1};
    m_shape_val["last_applied_force"] = {3};
    m_shape_val["pressure_tensor"] = {3, 3};
  }

private:
  auto get_slices_bounding_boxes() const {
    auto const &lattice = m_lb_fluid->get_lattice();
    auto const &slice_lower_corner = m_slice_lower_corner;
    auto const &slice_upper_corner = m_slice_upper_corner;
    assert(slice_upper_corner <= lattice.get_grid_dimensions());
    assert(slice_lower_corner >= Utils::Vector3i::broadcast(0));
    auto const sentinel = -(static_cast<int>(lattice.get_ghost_layers()) + 1);
    auto [local_lower_corner, local_upper_corner] =
        lattice.get_local_grid_range();
    for (auto const i : {0, 1, 2}) {
      if (local_lower_corner[i] >= slice_upper_corner[i] or
          slice_lower_corner[i] >= local_upper_corner[i]) {
        local_lower_corner[i] = sentinel;
        local_upper_corner[i] = sentinel;
      } else {
        if (slice_lower_corner[i] > local_lower_corner[i]) {
          local_lower_corner[i] = slice_lower_corner[i];
        }
        if (slice_upper_corner[i] < local_upper_corner[i]) {
          local_upper_corner[i] = slice_upper_corner[i];
        }
      }
    }
    return std::make_tuple(slice_lower_corner, slice_upper_corner,
                           local_lower_corner, local_upper_corner);
  }

  auto gather_slices_topology(Utils::Vector3i const &local_lower_corner,
                              Utils::Vector3i const &local_upper_corner) const {
    auto const &comm = context()->get_comm();
    std::vector<Utils::Vector3i> nodes_lower_corners;
    std::vector<Utils::Vector3i> nodes_upper_corners;
    boost::mpi::gather(comm, local_lower_corner, nodes_lower_corners, 0);
    boost::mpi::gather(comm, local_upper_corner, nodes_upper_corners, 0);
    return std::make_tuple(nodes_lower_corners, nodes_upper_corners);
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_slice_size") {
      return {m_slice_upper_corner - m_slice_lower_corner};
    }
    if (name == "get_slice_ranges") {
      return {std::vector<Variant>{m_slice_lower_corner, m_slice_upper_corner}};
    }
    if (name == "get_lb_sip") {
      return {m_lb_sip};
    }
    if (name == "get_value_shape") {
      auto const name = get_value<std::string>(params, "name");
      if (m_shape_val.count(name) == 0) {
        throw std::domain_error("Unknown fluid property '" + name + "'");
      }
      return m_shape_val.at(name);
    }
    if (name == "get_density") {
      return gather_3d(params, {1}, &::LBWalberlaBase::get_slice_density,
                       1. / m_conv_dens);
    }
    if (name == "set_density") {
      scatter_3d(params, {1}, &::LBWalberlaBase::set_slice_density,
                 m_conv_dens);
      return {};
    }
    if (name == "get_population") {
      return gather_3d(params, m_shape_val.at("population"),
                       &::LBWalberlaBase::get_slice_population);
    }
    if (name == "set_population") {
      scatter_3d(params, m_shape_val.at("population"),
                 &::LBWalberlaBase::set_slice_population);
      return {};
    }
    if (name == "get_velocity") {
      return gather_3d(params, {3}, &::LBWalberlaBase::get_slice_velocity,
                       1. / m_conv_velocity);
    }
    if (name == "set_velocity") {
      scatter_3d(params, {3}, &::LBWalberlaBase::set_slice_velocity,
                 m_conv_velocity);
      return {};
    }
    if (name == "get_is_boundary") {
      return gather_3d(params, {1}, &::LBWalberlaBase::get_slice_is_boundary);
    }
    if (name == "get_velocity_at_boundary") {
      return gather_3d(params, {1},
                       &::LBWalberlaBase::get_slice_velocity_at_boundary,
                       1. / m_conv_velocity);
    }
    if (name == "set_velocity_at_boundary") {
      scatter_3d(params, {1}, &::LBWalberlaBase::set_slice_velocity_at_boundary,
                 m_conv_velocity);
      return {};
    }
    if (name == "get_pressure_tensor") {
      return gather_3d(params, {3, 3},
                       &::LBWalberlaBase::get_slice_pressure_tensor,
                       1. / m_conv_press);
    }
    if (name == "get_last_applied_force") {
      return gather_3d(params, {3},
                       &::LBWalberlaBase::get_slice_last_applied_force,
                       1. / m_conv_force);
    }
    if (name == "set_last_applied_force") {
      scatter_3d(params, {3}, &::LBWalberlaBase::set_slice_last_applied_force,
                 m_conv_force);
      return {};
    }
    if (name == "get_lattice_speed") {
      return 1. / m_conv_velocity;
    }

    return {};
  }

private:
  template <typename T>
  Variant gather_3d(VariantMap const &params, std::vector<int> const &data_dims,
                    std::vector<T> (::LBWalberlaBase::*getter)(
                        Utils::Vector3i const &, Utils::Vector3i const &) const,
                    double units_conversion = 1.) const {
    auto const &comm = context()->get_comm();
    auto const [slice_lower_corner, slice_upper_corner, local_lower_corner,
                local_upper_corner] = get_slices_bounding_boxes();
    auto const [nodes_lower_corners, nodes_upper_corners] =
        gather_slices_topology(local_lower_corner, local_upper_corner);
    auto const data_size = std::accumulate(data_dims.cbegin(), data_dims.cend(),
                                           1, std::multiplies<>());
    auto const local_values =
        (m_lb_fluid.get()->*getter)(local_lower_corner, local_upper_corner);
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
      detail::flatten_grid(array, out, units_conversion);
      std::vector<int> shape = {m_shape.begin(), m_shape.end()};
      if (not(data_dims.size() == 1ul and data_dims[0] == 1)) {
        shape.insert(shape.end(), data_dims.begin(), data_dims.end());
      }
      if constexpr (std::is_same_v<T, int> or std::is_same_v<T, double>) {
        return {std::vector<Variant>{{Variant{out}, Variant{shape}}}};
      } else if constexpr (std::is_same_v<T,
                                          boost::optional<Utils::Vector3d>>) {
        std::vector<Variant> vec;
        vec.reserve(out.size());
        for (auto const &opt : out) {
          if (opt) {
            vec.emplace_back(Variant{*opt});
          } else {
            vec.emplace_back(Variant{None{}});
          }
        }
        return {std::vector<Variant>{{Variant{vec}, Variant{shape}}}};
      } else {
        return {std::vector<Variant>{
            {Variant{make_vector_of_variants(out)}, Variant{shape}}}};
      }
    }
    return {};
  }

  template <typename T>
  void scatter_3d(VariantMap const &params, std::vector<int> const &data_dims,
                  void (::LBWalberlaBase::*setter)(Utils::Vector3i const &,
                                                   Utils::Vector3i const &,
                                                   std::vector<T> const &),
                  double units_conversion = 1.) {
    auto const &comm = context()->get_comm();
    auto const [slice_lower_corner, slice_upper_corner, local_lower_corner,
                local_upper_corner] = get_slices_bounding_boxes();
    auto const [nodes_lower_corners, nodes_upper_corners] =
        gather_slices_topology(local_lower_corner, local_upper_corner);
    auto const data_size = std::accumulate(data_dims.cbegin(), data_dims.cend(),
                                           1, std::multiplies<>());
    std::vector<std::vector<T>> nodes_values(comm.size());
    if (comm.rank() == 0) {
      auto const &variant = params.at("values");
      std::vector<T> values;
      if constexpr (std::is_same_v<T, double>) {
        if (is_type<std::vector<int>>(variant)) {
          auto const values_int = get_value<std::vector<int>>(variant);
          values.reserve(values_int.size());
          for (auto const val : values_int) {
            values.emplace_back(static_cast<double>(val));
          }
        } else {
          values = get_value<std::vector<T>>(variant);
        }
      } else if constexpr (std::is_same_v<T,
                                          boost::optional<Utils::Vector3d>>) {
        auto const vector_variants = get_value<std::vector<Variant>>(variant);
        for (auto const &value : vector_variants) {
          if (is_none(value)) {
            values.emplace_back(boost::none);
          } else {
            values.emplace_back(get_value<Utils::Vector3d>(value));
          }
        }
      } else {
        values = get_value<std::vector<T>>(variant);
      }
      auto const dims = slice_upper_corner - slice_lower_corner;
      using index_range = boost::multi_array_types::index_range;
      using array_type = boost::multi_array<T, 4>;
      array_type array(boost::extents[dims[0]][dims[1]][dims[2]][data_size]);
      // populate the 3D array from the input flat array
      detail::unflatten_grid(array, values);
      // partition the 3D array into individual flat arrays for each MPI rank
      for (std::size_t rank = 0; rank < nodes_lower_corners.size(); ++rank) {
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
        detail::flatten_grid(view, nodes_values[rank], units_conversion);
      }
    }
    std::vector<T> local_values;
    boost::mpi::scatter(comm, nodes_values, local_values, 0);
    (m_lb_fluid.get()->*setter)(local_lower_corner, local_upper_corner,
                                local_values);
    m_lb_fluid->ghost_communication();
  }
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
