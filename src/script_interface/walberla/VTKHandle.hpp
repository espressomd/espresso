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

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/VTKHandle.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <boost/algorithm/string/join.hpp>

#include <algorithm>
#include <filesystem>
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
  std::weak_ptr<Field> m_field;
  ::LatticeModel::units_map m_units;
  std::vector<Variant> m_pending_arguments;

  [[nodiscard]] auto get_vtk_uid() const {
    return m_base_folder + '/' + m_identifier;
  }

  [[nodiscard]] std::shared_ptr<Field> get_field_instance() const {
    if (auto const field = m_field.lock()) {
      return field;
    }
    auto const has_lattice_expired = m_pending_arguments.empty();
    auto const err_expired = "Attempted access to an expired lattice object";
    auto const err_detached = "This VTK object isn't attached to a lattice";
    throw std::runtime_error(has_lattice_expired ? err_expired : err_detached);
  }

  virtual std::unordered_map<std::string, int> const &get_obs_map() const = 0;

  [[nodiscard]] auto get_valid_observable_names() const {
    std::vector<std::string> names{};
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
    std::vector<Variant> observables{};
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
        {"vtk_uid", read_only, [this]() { return get_vtk_uid(); }},
        {"identifier", read_only, [this]() { return m_identifier; }},
        {"base_folder", read_only, [this]() { return m_base_folder; }},
        {"prefix", read_only, [this]() { return m_prefix; }},
        {"observables", read_only,
         [this]() { return serialize_obs_flag(m_flag_obs); }},
        {"execution_count", read_only,
         [this]() { return m_vtk_handle->execution_count; }},
        {"units", read_only,
         [this]() { return make_unordered_map_of_variants(m_units); }},
    });
  }

private:
  void do_construct(VariantMap const &params) override {
    m_delta_N = get_value<int>(params, "delta_N");
    m_identifier = get_value<std::string>(params, "identifier");
    m_base_folder = get_value<std::string>(params, "base_folder");
    m_prefix = get_value<std::string>(params, "prefix");
    auto const is_enabled = get_value<bool>(params, "enabled");
    auto const execution_count = get_value<int>(params, "execution_count");
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_flag_obs = deserialize_obs_flag(
          get_value<std::vector<std::string>>(params, "observables"));
      if (m_delta_N < 0) {
        throw std::domain_error("Parameter 'delta_N' must be >= 0");
      }
      if (m_identifier.empty()) {
        throw std::domain_error("Parameter 'identifier' cannot be empty");
      }
      if (m_identifier.find(std::filesystem::path::preferred_separator) !=
          std::string::npos) {
        throw std::invalid_argument(
            "Parameter 'identifier' cannot be a filepath");
      }
    });
    m_pending_arguments.emplace_back(is_enabled);
    m_pending_arguments.emplace_back(execution_count);
  }

protected:
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

public:
  void detach_from_lattice() { m_field.reset(); }

  void attach_to_lattice(std::weak_ptr<Field> field,
                         ::LatticeModel::units_map const &units) {
    auto const was_attached_once = m_pending_arguments.empty();
    if (not m_field.expired()) {
      throw std::runtime_error("Cannot attach VTK object to multiple lattices");
    }
    if (was_attached_once) {
      throw std::runtime_error("Detached VTK objects cannot be attached again");
    }
    assert(m_pending_arguments.size() == 2u);
    auto const is_enabled = get_value<bool>(m_pending_arguments[0]);
    auto const execution_count = get_value<int>(m_pending_arguments[1]);
    m_units = units;
    m_field = field;
    auto instance = get_field_instance();
    m_vtk_handle =
        instance->create_vtk(m_delta_N, execution_count, m_flag_obs, m_units,
                             m_identifier, m_base_folder, m_prefix);
    if (m_delta_N and not is_enabled) {
      instance->switch_vtk(get_vtk_uid(), false);
    }
    m_pending_arguments.clear();
  }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
