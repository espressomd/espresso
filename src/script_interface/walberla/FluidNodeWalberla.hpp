/*
 * Copyright (C) 2021-2022 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaBase.hpp>

#include "FluidWalberla.hpp"
#include "optional_reduction.hpp"

#include "core/communication.hpp"
#include "core/errorhandling.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class FluidNodeWalberla : public AutoParameters<FluidNodeWalberla> {
  // TODO WALBERLA: revert commit 2f0c490b8e1bb4ab3e to use a weak_ptr
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  Utils::Vector3i m_index;
  Utils::Vector3i m_grid_size;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;

public:
  FluidNodeWalberla() {
    add_parameters(
        {{"_index", AutoParameter::read_only, [this]() { return (m_index); }}});
  }

  void do_construct(VariantMap const &params) override {
    try {
      auto const lb_sip =
          get_value<std::shared_ptr<FluidWalberla>>(params, "lb_sip");
      m_lb_fluid = lb_sip->lb_fluid().lock();
      assert(m_lb_fluid);
      auto const &lb_params = lb_sip->lb_params().lock();
      assert(lb_params);
      auto const tau = lb_params->get_tau();
      auto const agrid = lb_params->get_agrid();
      m_grid_size = m_lb_fluid->lattice().get_grid_dimensions();
      m_index = get_value<Utils::Vector3i>(params, "index");
      if (not is_index_valid(m_index)) {
        throw std::out_of_range("Index error");
      }
      m_conv_dens = Utils::int_pow<3>(agrid);
      m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
      m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
      m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    } catch (const std::exception &e) {
      runtimeErrorMsg() << "FluidNodeWalberla failed: " << e.what();
      m_lb_fluid.reset();
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "override_index") {
      // this hidden feature is used to iterate a LB slice without
      // rebuilding a FluidNodeWalberla for each node in the slice
      auto const index = get_value<Utils::Vector3i>(params, "index");
      if (not is_index_valid(index)) {
        return ES_ERROR;
      }
      m_index = index;
      return ES_OK;
    } else if (name == "set_velocity_at_boundary") {
      if (is_none(params.at("value"))) {
        m_lb_fluid->remove_node_from_boundary(m_index);
        m_lb_fluid->ghost_communication();
      } else {
        auto const u =
            get_value<Utils::Vector3d>(params, "value") * m_conv_velocity;
        m_lb_fluid->set_node_velocity_at_boundary(m_index, u);
        m_lb_fluid->ghost_communication();
      }
    } else if (name == "get_velocity_at_boundary") {
      auto result = m_lb_fluid->get_node_is_boundary(m_index);
      bool is_boundary = (result) ? *result : false;
      is_boundary =
          boost::mpi::all_reduce(comm_cart, is_boundary, std::logical_or<>());
      if (is_boundary) {
        auto result = m_lb_fluid->get_node_velocity_at_boundary(m_index);
        return optional_reduction_with_conversion(result, m_conv_velocity);
      }
      return Variant{None{}};
    } else if (name == "set_velocity") {
      auto const u =
          get_value<Utils::Vector3d>(params, "value") * m_conv_velocity;
      m_lb_fluid->set_node_velocity(m_index, u);
      m_lb_fluid->ghost_communication();
    } else if (name == "get_velocity") {
      auto result = m_lb_fluid->get_node_velocity(m_index);
      return optional_reduction_with_conversion(result, m_conv_velocity);
    } else if (name == "set_density") {
      auto const dens = get_value<double>(params, "value");
      m_lb_fluid->set_node_density(m_index, dens * m_conv_dens);
      m_lb_fluid->ghost_communication();
    } else if (name == "get_density") {
      auto const result = m_lb_fluid->get_node_density(m_index);
      return optional_reduction_with_conversion(result, m_conv_dens);
    } else if (name == "set_population") {
      auto const pop = get_value<std::vector<double>>(params, "value");
      m_lb_fluid->set_node_pop(m_index, pop);
      m_lb_fluid->ghost_communication();
    } else if (name == "get_population") {
      auto const result = m_lb_fluid->get_node_pop(m_index);
      return optional_reduction_with_conversion(result);
    } else if (name == "get_is_boundary") {
      auto const result = m_lb_fluid->get_node_is_boundary(m_index);
      return optional_reduction_with_conversion(result);
    } else if (name == "get_boundary_force") {
      auto result = m_lb_fluid->get_node_is_boundary(m_index);
      bool is_boundary = (result) ? *result : false;
      is_boundary =
          boost::mpi::all_reduce(comm_cart, is_boundary, std::logical_or<>());
      if (is_boundary) {
        auto result = m_lb_fluid->get_node_boundary_force(m_index);
        return optional_reduction_with_conversion(result, m_conv_force);
      }
      return Variant{None{}};
    } else if (name == "get_pressure_tensor") {
      auto const result = m_lb_fluid->get_node_pressure_tensor(m_index);
      auto value = boost::optional<std::vector<double>>{};
      if (result) {
        value = (*result / m_conv_press).as_vector();
      }
      auto const variant = optional_reduction_with_conversion(value);
      if (context()->is_head_node()) {
        auto const vec = get_value<std::vector<double>>(variant);
        auto tensor = Utils::Matrix<double, 3, 3>{};
        std::copy(vec.begin(), vec.end(), tensor.m_data.begin());
        return std::vector<Variant>{tensor.row<0>().as_vector(),
                                    tensor.row<1>().as_vector(),
                                    tensor.row<2>().as_vector()};
      }
      return {};
    } else if (name == "set_last_applied_force") {
      auto const f = get_value<Utils::Vector3d>(params, "value");
      m_lb_fluid->set_node_last_applied_force(m_index, f * m_conv_force);
      m_lb_fluid->ghost_communication();
    } else if (name == "get_last_applied_force") {
      auto const result = m_lb_fluid->get_node_last_applied_force(m_index);
      return optional_reduction_with_conversion(result, m_conv_force);
    } else if (name == "get_lattice_speed") {
      return 1. / m_conv_velocity;
    }

    return {};
  }

private:
  inline bool is_index_valid(Utils::Vector3i const &index) const {
    return index < m_grid_size and index >= Utils::Vector3i{};
  }
};
} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
