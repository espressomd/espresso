/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "script_interface/Exception.hpp"
#include "script_interface/ObjectHandle.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {

/**
 * @brief Bind parameters in the script interface.
 *
 * This class implements @c ScriptInterfaceBase, binding
 * the parameters added by add_parameters or by the constructor.
 * To use it, derive from this class and add parameters. For example,
 * given a class A
 * ~~~{.cpp}
 * class A {
 * public:
 *   int i() { return m_i; }
 * private:
 *   int m_i;
 * };
 * ~~~
 * that should have @c i exposed, this can be achieved by extending it
 * like this:
 * ~~~{.cpp}
 * class A : public AutoParameters {
 * public:
 *   A() : AutoParameters({"name_for_i", i}) {}
 *   int i() { return m_i; }
 * private:
 *   int m_i;
 * };
 * ~~~
 *
 * If there is more complicated logic needed, specific setters and
 * getters can be provided. E.g. given a class B like
 * ~~~{.cpp}
 * class B {
 * public:
 *   void set_i(int);
 *   int get_i();
 * private:
 *   int m_i;
 * };
 * ~~~
 * we can use a lambdas to set and get the parameter like this:
 * ~~~{.cpp}
 * class B : public AutoParameters {
 * public:
 *   B() : AutoParameters({"name_for_i",
 *                         [this](Variant const& v) {set_i(get_value<int>(v));},
 *                         [this]() {return get_i();}
 *                        }) {}
 *   void set_i(int);
 *   int get_i();
 * private:
 *   int m_i;
 * };
 * ~~~
 * (this has to be captured in the lambdas to have access to the member
 * functions of the class).
 */
template <typename Derived, typename Base = ObjectHandle>
class AutoParameters : public Base {
  static_assert(std::is_base_of_v<ObjectHandle, Base>);

public:
  /** @brief Exception thrown when accessing an unknown parameter */
  struct UnknownParameter : public Exception {
    explicit UnknownParameter(std::string const &name)
        : Exception("Unknown parameter '" + name + "'.") {}
  };

  /** @brief Exception thrown when writing to a read-only parameter */
  struct WriteError : public Exception {
    explicit WriteError(std::string const &name)
        : Exception("Parameter '" + name + "' is read-only.") {}
  };

protected:
  // NOLINTNEXTLINE(bugprone-crtp-constructor-accessibility)
  AutoParameters() = default;
  // NOLINTNEXTLINE(bugprone-crtp-constructor-accessibility)
  explicit AutoParameters(std::vector<AutoParameter> &&params) {
    add_parameters(std::move(params));
  }
  ~AutoParameters() override = default;

  bool has_parameter(std::string const &name) const override {
    return m_parameters.contains(name);
  }

  void add_parameters(std::vector<AutoParameter> &&params) {
    for (auto const &p : params) {
      if (m_parameters.contains(p.name)) {
        m_parameters.erase(p.name);
        if (auto const it = std::ranges::find(m_key_order, p.name);
            it != m_key_order.end()) {
          m_key_order.erase(it);
        }
      }
      m_key_order.emplace_back(p.name);
      m_parameters.emplace(p.name, std::move(p));
    }
  }

  auto const &get_parameter_insertion_order() const { return m_key_order; }

public:
  /* ObjectHandle implementation */
  std::vector<std::string_view> valid_parameters() const final {
    auto const view = std::views::elements<0>(m_parameters);
    return {view.begin(), view.end()};
  }

  Variant get_parameter(const std::string &name) const final {
    try {
      return m_parameters.at(name).get();
    } catch (std::out_of_range const &) {
      throw UnknownParameter{name};
    }
  }

  void do_set_parameter(const std::string &name, const Variant &value) final {
    try {
      m_parameters.at(name).set(value);
    } catch (AutoParameter::WriteError const &) {
      throw WriteError{name};
    } catch (std::out_of_range const &) {
      throw UnknownParameter{name};
    }
  }

  std::vector<std::pair<std::string, Variant>>
  serialize_parameters() const final {
    std::vector<std::pair<std::string, Variant>> parameter_pack{};
    auto const params = this->get_parameters();
    for (auto const &key : m_key_order) {
      parameter_pack.emplace_back(key, params.at(key));
    }
    return parameter_pack;
  }

private:
  /** @brief Data structure for the stored parameters. */
  std::unordered_map<std::string, AutoParameter> m_parameters;
  /** @brief Keep track of the insertion order of parameters. */
  std::vector<std::string> m_key_order;
};
} // namespace ScriptInterface
