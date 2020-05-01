/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETERS_HPP

#include "script_interface/ScriptInterfaceBase.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <unordered_map>
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
template <typename Derived, typename Base = ScriptInterfaceBase>
class AutoParameters : public Base {
  static_assert(std::is_base_of<ScriptInterfaceBase, Base>::value, "");

public:
  /** @brief Exception thrown when accessing an unknown parameter */
  struct UnknownParameter : public std::runtime_error {
    explicit UnknownParameter(std::string const &name)
        : runtime_error("Unknown parameter '" + name + "'.") {}
  };

  /** @brief Exception thrown when writing to a read-only parameter */
  struct WriteError : public std::runtime_error {
    explicit WriteError(std::string const &name)
        : runtime_error("Parameter " + name + " is read-only.") {}
  };

protected:
  AutoParameters() = default;
  explicit AutoParameters(std::vector<AutoParameter> &&params) {
    add_parameters(std::move(params));
  }

  void add_parameters(std::vector<AutoParameter> &&params) {
    for (auto const &p : params) {
      m_parameters.emplace(std::make_pair(p.name, p));
    }
  }

public:
  /* ScriptInterfaceBase implementation */
  Utils::Span<const boost::string_ref> valid_parameters() const final {
    static std::vector<boost::string_ref> valid_params;
    valid_params.clear();

    for (auto const &p : m_parameters) {
      valid_params.emplace_back(p.first);
    }

    return valid_params;
  }

  Variant get_parameter(const std::string &name) const final {
    try {
      return m_parameters.at(name).get();
    } catch (std::out_of_range const &e) {
      throw UnknownParameter{name};
    }
  }

  void set_parameter(const std::string &name, const Variant &value) final {
    try {
      m_parameters.at(name).set(value);
    } catch (AutoParameter::WriteError const &e) {
      throw WriteError{name};
    } catch (std::out_of_range const &e) {
      throw UnknownParameter{name};
    }
  }

private:
  std::unordered_map<std::string, AutoParameter> m_parameters;
};
} // namespace ScriptInterface

#endif
