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

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

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
#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class FluidNodeWalberla
    : public AutoParameters<FluidNodeWalberla, LatticeIndices> {
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
        {{"_index", AutoParameter::read_only, [this]() { return m_index; }}});
  }

  void do_construct(VariantMap const &params) override {
    try {
      auto const lb_sip =
          get_value<std::shared_ptr<FluidWalberla>>(params, "parent_sip");
      m_lb_fluid = lb_sip->lb_fluid().lock();
      assert(m_lb_fluid);
      auto const &lb_params = lb_sip->lb_params().lock();
      assert(lb_params);
      auto const tau = lb_params->get_tau();
      auto const agrid = lb_params->get_agrid();
      m_conv_dens = Utils::int_pow<3>(agrid);
      m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
      m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
      m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    } catch (std::exception const &e) {
      if (context()->is_head_node()) {
        runtimeErrorMsg() << "FluidNodeWalberla failed: " << e.what();
      }
      m_lb_fluid.reset();
      return;
    }
    m_grid_size = m_lb_fluid->get_lattice().get_grid_dimensions();
    m_index = get_mapped_index(get_value<Utils::Vector3i>(params, "index"),
                               m_grid_size);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "override_index") {
      // this hidden feature is used to iterate a LB slice without
      // rebuilding a FluidNodeWalberla for each node in the slice
      auto const index = get_value<Utils::Vector3i>(params, "index");
      if (not is_index_valid(index, m_grid_size)) {
        return ES_ERROR;
      }
      m_index = index;
      return ES_OK;
    }
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
      if (is_boundary_all_reduce(m_lb_fluid->get_node_is_boundary(m_index))) {
        auto const result = m_lb_fluid->get_node_velocity_at_boundary(m_index);
        return mpi_reduce_optional(context()->get_comm(), result) /
               m_conv_velocity;
      }
      return Variant{None{}};
    }
    if (name == "set_velocity") {
      auto const u =
          get_value<Utils::Vector3d>(params, "value") * m_conv_velocity;
      m_lb_fluid->set_node_velocity(m_index, u);
      m_lb_fluid->ghost_communication();
      return {};
    }
    if (name == "get_velocity") {
      auto const result = m_lb_fluid->get_node_velocity(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) /
             m_conv_velocity;
    }
    if (name == "set_density") {
      auto const dens = get_value<double>(params, "value");
      m_lb_fluid->set_node_density(m_index, dens * m_conv_dens);
      m_lb_fluid->ghost_communication();
      return {};
    }
    if (name == "get_density") {
      auto const result = m_lb_fluid->get_node_density(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) / m_conv_dens;
    }
    if (name == "set_population") {
      auto const pop = get_value<std::vector<double>>(params, "value");
      m_lb_fluid->set_node_pop(m_index, pop);
      m_lb_fluid->ghost_communication();
      return {};
    }
    if (name == "get_population") {
      auto const result = m_lb_fluid->get_node_pop(m_index);
      return mpi_reduce_optional(context()->get_comm(), result);
    }
    if (name == "get_is_boundary") {
      auto const result = m_lb_fluid->get_node_is_boundary(m_index);
      return mpi_reduce_optional(context()->get_comm(), result);
    }
    if (name == "get_boundary_force") {
      if (is_boundary_all_reduce(m_lb_fluid->get_node_is_boundary(m_index))) {
        auto result = m_lb_fluid->get_node_boundary_force(m_index);
        return mpi_reduce_optional(context()->get_comm(), result) /
               m_conv_force;
      }
      return Variant{None{}};
    }
    if (name == "get_pressure_tensor") {
      auto const result = m_lb_fluid->get_node_pressure_tensor(m_index);
      auto value = boost::optional<std::vector<double>>{};
      if (result) {
        value = (*result / m_conv_press).as_vector();
      }
      auto const vec = mpi_reduce_optional(context()->get_comm(), value);
      if (context()->is_head_node()) {
        auto tensor = Utils::Matrix<double, 3, 3>{};
        std::copy(vec.begin(), vec.end(), tensor.m_data.begin());
        return std::vector<Variant>{tensor.row<0>().as_vector(),
                                    tensor.row<1>().as_vector(),
                                    tensor.row<2>().as_vector()};
      }
      return {};
    }
    if (name == "set_last_applied_force") {
      auto const f = get_value<Utils::Vector3d>(params, "value");
      m_lb_fluid->set_node_last_applied_force(m_index, f * m_conv_force);
      m_lb_fluid->ghost_communication();
      return {};
    }
    if (name == "get_last_applied_force") {
      auto const result = m_lb_fluid->get_node_last_applied_force(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) / m_conv_force;
    }
    if (name == "get_lattice_speed") {
      return 1. / m_conv_velocity;
    }

    return {};
  }

private:
  bool is_boundary_all_reduce(boost::optional<bool> const &is_boundary) const {
    return boost::mpi::all_reduce(context()->get_comm(),
                                  is_boundary ? *is_boundary : false,
                                  std::logical_or<>());
  }
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
