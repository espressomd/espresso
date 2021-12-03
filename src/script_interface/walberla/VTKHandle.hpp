/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_VTKWALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_VTKWALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaBase.hpp>
#include <walberla_bridge/VTKHandle.hpp>

#include "core/communication.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

class VTKHandle : public AutoParameters<VTKHandle> {
  std::shared_ptr<::VTKHandle> m_vtk_handle;
  std::weak_ptr<LBWalberlaBase> m_lb_fluid;
  int m_delta_N;
  int m_flag_obs;
  std::string m_identifier;
  std::string m_base_folder;
  std::string m_prefix;

  auto get_vtk_uid() const { return m_base_folder + '/' + m_identifier; }

  std::shared_ptr<LBWalberlaBase> get_lb_fluid() {
    auto lb_walberla_instance_handle = m_lb_fluid.lock();
    if (!lb_walberla_instance_handle) {
      throw std::runtime_error(
          "Attempted access to uninitialized LBWalberla instance.");
    }
    return lb_walberla_instance_handle;
  }

public:
  VTKHandle() {
    constexpr auto read_only = AutoParameter::read_only;
    add_parameters({
        {"enabled", read_only, [this]() { return m_vtk_handle->enabled; }},
        {"delta_N", read_only, [this]() { return m_delta_N; }},
        {"flag_obs", read_only, [this]() { return m_flag_obs; }},
        {"vtk_uid", read_only, [this]() { return get_vtk_uid(); }},
        {"identifier", read_only, [this]() { return m_identifier; }},
        {"base_folder", read_only, [this]() { return m_base_folder; }},
        {"prefix", read_only, [this]() { return m_prefix; }},
        {"observables", read_only,
         [this]() {
           std::vector<Variant> observables;
           if (m_flag_obs & static_cast<int>(OutputVTK::density)) {
             observables.emplace_back(std::string("density"));
           }
           if (m_flag_obs & static_cast<int>(OutputVTK::velocity_vector)) {
             observables.emplace_back(std::string("velocity_vector"));
           }
           if (m_flag_obs & static_cast<int>(OutputVTK::pressure_tensor)) {
             observables.emplace_back(std::string("pressure_tensor"));
           }
           return observables;
         }},
        {"execution_count", read_only,
         [this]() { return m_vtk_handle->execution_count; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    if (params.count("flag_obs")) {
      // object built from a checkpoint
      m_flag_obs = get_value<int>(params, "flag_obs");
      // TODO WALBERLA: here we assume all VTK objects belong to the
      // same LB actor, which should always be the case (users cannot
      // create multiple actors with different VTK objects, unless
      // the system is checkpointed multiple times)
    } else {
      // object built by the user
      m_flag_obs = 0;
      auto vec = get_value<std::vector<Variant>>(params, "observables");
      for (auto const &variant : vec) {
        auto const observable = boost::get<std::string>(variant);
        if (observable == "density")
          m_flag_obs |= static_cast<int>(OutputVTK::density);
        if (observable == "velocity_vector")
          m_flag_obs |= static_cast<int>(OutputVTK::velocity_vector);
        if (observable == "pressure_tensor")
          m_flag_obs |= static_cast<int>(OutputVTK::pressure_tensor);
      }
    }
    m_delta_N = get_value<int>(params, "delta_N");
    m_identifier = get_value<std::string>(params, "identifier");
    m_base_folder = get_value<std::string>(params, "base_folder");
    m_prefix = get_value<std::string>(params, "prefix");
    auto const execution_count = get_value<int>(params, "execution_count");
    try {
      std::shared_ptr<LBWalberlaBase> lb_fluid;
      if (params.count("lb_fluid")) {
        m_lb_fluid =
            get_value<std::shared_ptr<FluidWalberla>>(params, "lb_fluid")
                ->lb_fluid();
        lb_fluid = get_lb_fluid();
      } else {
        lb_fluid = ::lb_walberla();
        m_lb_fluid = std::weak_ptr<LBWalberlaBase>{lb_fluid};
      }
      m_vtk_handle =
          lb_fluid->create_vtk(m_delta_N, execution_count, m_flag_obs,
                               m_identifier, m_base_folder, m_prefix);
      if (m_delta_N and not get_value<bool>(params, "enabled")) {
        lb_fluid->switch_vtk(get_vtk_uid(), false);
      }
    } catch (std::exception const &err) {
      if (this_node == 0) {
        throw;
      }
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    try {
      if (name == "enable") {
        if (m_delta_N == 0) {
          throw std::runtime_error("Manual VTK callbacks cannot be enabled");
        }
        get_lb_fluid()->switch_vtk(get_vtk_uid(), true);
      }

      if (name == "disable") {
        if (m_delta_N == 0) {
          throw std::runtime_error("Manual VTK callbacks cannot be disabled");
        }
        get_lb_fluid()->switch_vtk(get_vtk_uid(), false);
      }

      if (name == "write") {
        if (m_delta_N) {
          throw std::runtime_error("Automatic VTK callbacks cannot be "
                                   "triggered manually");
        }
        get_lb_fluid()->write_vtk(get_vtk_uid());
      }
    } catch (std::exception const &err) {
      if (this_node == 0) {
        throw;
      }
    }

    return {};
  }
};

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
