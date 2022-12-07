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

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/VTKHandle.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include "core/grid_based_algorithms/lb_walberla_instance.hpp"

#include "EKSpecies.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <boost/algorithm/string/join.hpp>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

template <class Field>
class VTKHandleBase : public AutoParameters<VTKHandleBase<Field>> {
private:
  int m_delta_N;
  int m_flag_obs;
  std::string m_identifier;
  std::string m_base_folder;
  std::string m_prefix;
  std::shared_ptr<::VTKHandle> m_vtk_handle;

  [[nodiscard]] auto get_vtk_uid() const {
    return m_base_folder + '/' + m_identifier;
  }

protected:
  virtual void setup_field_instance(VariantMap const &params) = 0;
  [[nodiscard]] virtual std::shared_ptr<Field> get_field_instance() = 0;
  virtual std::unordered_map<std::string, int> const &get_obs_map() const = 0;

private:
  [[nodiscard]] auto get_valid_observable_names() const {
    std::vector<std::string> names;
    for (auto const &kv : get_obs_map()) {
      names.emplace_back(kv.first);
    }
    std::sort(names.begin(), names.end());
    return names;
  }

  [[nodiscard]] int
  deserialize_obs_flag(std::vector<std::string> const &names) const {
    int flag{0};
    auto const &obs_map = get_obs_map();
    for (auto const &name : names) {
      if (obs_map.count(name) == 0) {
        auto const valid_names = get_valid_observable_names();
        std::stringstream message;
        message << "Only the following VTK observables are supported: ["
                << "'" << boost::algorithm::join(valid_names, "', '") << "'"
                << "], got '" << name << "'";
        throw std::invalid_argument(message.str());
      }
      flag |= obs_map.at(name);
    }
    return flag;
  }

  [[nodiscard]] Variant serialize_obs_flag(int flag) const {
    std::vector<Variant> observables;
    for (auto const &kv : get_obs_map()) {
      if (flag & kv.second) {
        observables.emplace_back(kv.first);
      }
    }
    return observables;
  }

public:
  VTKHandleBase() {
    constexpr auto read_only = AutoParameter::read_only;
    AutoParameters<VTKHandleBase<Field>>::add_parameters({
        {"enabled", read_only, [this]() { return m_vtk_handle->enabled; }},
        {"delta_N", read_only, [this]() { return m_delta_N; }},
        {"flag_obs", read_only, [this]() { return m_flag_obs; }},
        {"vtk_uid", read_only, [this]() { return get_vtk_uid(); }},
        {"identifier", read_only, [this]() { return m_identifier; }},
        {"base_folder", read_only, [this]() { return m_base_folder; }},
        {"prefix", read_only, [this]() { return m_prefix; }},
        {"observables", read_only,
         [this]() { return serialize_obs_flag(m_flag_obs); }},
        {"execution_count", read_only,
         [this]() { return m_vtk_handle->execution_count; }},
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "enable") {
      ObjectHandle::context()->parallel_try_catch([&]() {
        if (m_delta_N == 0) {
          throw std::runtime_error("Manual VTK callbacks cannot be enabled");
        }
        get_field_instance()->switch_vtk(get_vtk_uid(), true);
      });
      return {};
    }
    if (name == "disable") {
      ObjectHandle::context()->parallel_try_catch([&]() {
        if (m_delta_N == 0) {
          throw std::runtime_error("Manual VTK callbacks cannot be disabled");
        }
        get_field_instance()->switch_vtk(get_vtk_uid(), false);
      });
      return {};
    }
    if (name == "write") {
      ObjectHandle::context()->parallel_try_catch([&]() {
        if (m_delta_N) {
          throw std::runtime_error("Automatic VTK callbacks cannot be "
                                   "triggered manually");
        }
        get_field_instance()->write_vtk(get_vtk_uid());
      });
      return {};
    }
    if (name == "get_valid_observable_names") {
      return make_vector_of_variants(get_valid_observable_names());
    }

    return {};
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
      m_flag_obs = deserialize_obs_flag(
          get_value<std::vector<std::string>>(params, "observables"));
    }
    m_delta_N = get_value<int>(params, "delta_N");
    m_identifier = get_value<std::string>(params, "identifier");
    m_base_folder = get_value<std::string>(params, "base_folder");
    m_prefix = get_value<std::string>(params, "prefix");
    auto const execution_count = get_value<int>(params, "execution_count");
    ObjectHandle::context()->parallel_try_catch([&]() {
      setup_field_instance(params);
      auto &field_instance = *get_field_instance();
      m_vtk_handle =
          field_instance.create_vtk(m_delta_N, execution_count, m_flag_obs,
                                    m_identifier, m_base_folder, m_prefix);
      if (m_delta_N and not get_value<bool>(params, "enabled")) {
        field_instance.switch_vtk(get_vtk_uid(), false);
      }
    });
  }
};

class VTKHandle : public VTKHandleBase<LBWalberlaBase> {
  static std::unordered_map<std::string, int> const obs_map;
  std::weak_ptr<LBWalberlaBase> m_lb_fluid;

  [[nodiscard]] std::shared_ptr<LBWalberlaBase> get_field_instance() override {
    auto lb_walberla_instance_handle = m_lb_fluid.lock();
    if (!lb_walberla_instance_handle) {
      throw std::runtime_error(
          "Attempted access to uninitialized LBWalberla instance.");
    }
    return lb_walberla_instance_handle;
  }

  void setup_field_instance(VariantMap const &params) override {
    if (params.count("lb_fluid")) {
      m_lb_fluid = get_value<std::shared_ptr<FluidWalberla>>(params, "lb_fluid")
                       ->lb_fluid();
    } else {
      m_lb_fluid = std::weak_ptr<LBWalberlaBase>{::lb_walberla()};
    }
  }

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

std::unordered_map<std::string, int> const VTKHandle::obs_map = {
    {"density", static_cast<int>(OutputVTK::density)},
    {"velocity_vector", static_cast<int>(OutputVTK::velocity_vector)},
    {"pressure_tensor", static_cast<int>(OutputVTK::pressure_tensor)},
};

class EKVTKHandle : public VTKHandleBase<EKinWalberlaBase> {
  static std::unordered_map<std::string, int> const obs_map;
  std::shared_ptr<EKinWalberlaBase> m_ekinstance;
  std::shared_ptr<EKSpecies> m_ekspecies;

  [[nodiscard]] std::shared_ptr<EKinWalberlaBase>
  get_field_instance() override {
    return m_ekinstance;
  }

  // TODO WALBERLA
  /*
public:
  EKVTKHandle() {
    constexpr auto read_only = AutoParameter::read_only;
    add_parameters({
        {"species", read_only, [this]() { return m_ekspecies; }},
    });
  }
  */

  void setup_field_instance(VariantMap const &params) override {
    m_ekspecies = get_value<std::shared_ptr<EKSpecies>>(params, "species");
    m_ekinstance = m_ekspecies->get_ekinstance();
  }

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

std::unordered_map<std::string, int> const EKVTKHandle::obs_map = {
    {"density", static_cast<int>(EKOutputVTK::density)},
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
