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

#include "config/config.hpp"

#ifdef WALBERLA

#include "LBFluidNode.hpp"

#include <script_interface/communication.hpp>

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/optional.hpp>

#include <string>
#include <vector>

namespace ScriptInterface::walberla {

static bool is_boundary_all_reduce(boost::mpi::communicator const &comm,
                                   boost::optional<bool> const &is_boundary) {
  return boost::mpi::all_reduce(comm, is_boundary ? *is_boundary : false,
                                std::logical_or<>());
}

Variant LBFluidNode::do_call_method(std::string const &name,
                                    VariantMap const &params) {
  if (name == "set_velocity_at_boundary") {
    if (is_none(params.at("value"))) {
      m_lb_fluid->remove_node_from_boundary(m_index);
      m_lb_fluid->ghost_communication();
    } else {
      auto const u =
          get_value<Utils::Vector3d>(params, "value") * m_conv_velocity;
      m_lb_fluid->set_node_velocity_at_boundary(m_index, u);
      m_lb_fluid->ghost_communication();
    }
    return {};
  }
  if (name == "get_velocity_at_boundary") {
    auto const boundary_opt = m_lb_fluid->get_node_is_boundary(m_index);
    if (is_boundary_all_reduce(context()->get_comm(), boundary_opt)) {
      auto const result = m_lb_fluid->get_node_velocity_at_boundary(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) /
             m_conv_velocity;
    }
    return Variant{None{}};
  }
  if (name == "get_density") {
    auto const result = m_lb_fluid->get_node_density(m_index);
    return mpi_reduce_optional(context()->get_comm(), result) / m_conv_dens;
  }
  if (name == "set_density") {
    auto const dens = get_value<double>(params, "value");
    m_lb_fluid->set_node_density(m_index, dens * m_conv_dens);
    m_lb_fluid->ghost_communication();
    return {};
  }
  if (name == "get_population") {
    auto const result = m_lb_fluid->get_node_population(m_index);
    return mpi_reduce_optional(context()->get_comm(), result);
  }
  if (name == "set_population") {
    auto const pop = get_value<std::vector<double>>(params, "value");
    m_lb_fluid->set_node_population(m_index, pop);
    m_lb_fluid->ghost_communication();
    return {};
  }
  if (name == "get_velocity") {
    auto const result = m_lb_fluid->get_node_velocity(m_index);
    return mpi_reduce_optional(context()->get_comm(), result) / m_conv_velocity;
  }
  if (name == "set_velocity") {
    auto const u =
        get_value<Utils::Vector3d>(params, "value") * m_conv_velocity;
    m_lb_fluid->set_node_velocity(m_index, u);
    m_lb_fluid->ghost_communication();
    return {};
  }
  if (name == "get_is_boundary") {
    auto const result = m_lb_fluid->get_node_is_boundary(m_index);
    return mpi_reduce_optional(context()->get_comm(), result);
  }
  if (name == "get_boundary_force") {
    auto const boundary_opt = m_lb_fluid->get_node_is_boundary(m_index);
    if (is_boundary_all_reduce(context()->get_comm(), boundary_opt)) {
      auto result = m_lb_fluid->get_node_boundary_force(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) / m_conv_force;
    }
    return Variant{None{}};
  }
  if (name == "get_pressure_tensor" or name == "get_pressure_tensor_neq") {
    auto const result = m_lb_fluid->get_node_pressure_tensor(m_index);
    auto value = boost::optional<std::vector<double>>{};
    if (result) {
      value = (*result / m_conv_press).as_vector();
    }
    auto vec = mpi_reduce_optional(context()->get_comm(), value);
    if (context()->is_head_node()) {
      if (name == "get_pressure_tensor_neq") {
        auto constexpr c_sound_sq = 1. / 3.;
        auto const density = m_lb_fluid->get_density();
        auto const diagonal_term = density * c_sound_sq / m_conv_press;
        vec[0] -= diagonal_term;
        vec[4] -= diagonal_term;
        vec[8] -= diagonal_term;
      }
      auto tensor = Utils::Matrix<double, 3, 3>{};
      std::copy(vec.begin(), vec.end(), tensor.m_data.begin());
      return std::vector<Variant>{tensor.row<0>().as_vector(),
                                  tensor.row<1>().as_vector(),
                                  tensor.row<2>().as_vector()};
    }
    return {};
  }
  if (name == "get_last_applied_force") {
    auto const result = m_lb_fluid->get_node_last_applied_force(m_index);
    return mpi_reduce_optional(context()->get_comm(), result) / m_conv_force;
  }
  if (name == "set_last_applied_force") {
    auto const f = get_value<Utils::Vector3d>(params, "value");
    m_lb_fluid->set_node_last_applied_force(m_index, f * m_conv_force);
    m_lb_fluid->ghost_communication();
    return {};
  }
  if (name == "get_lattice_speed") {
    return 1. / m_conv_velocity;
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
