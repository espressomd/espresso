/*
 * Copyright (C) 2022 The ESPResSo project
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

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLES)

#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include "scafacos.hpp"

#include "core/scafacos/ScafacosContextBase.hpp"

#include <utils/demangle.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Scafacos {

std::vector<std::string> available_methods() {
  return ScafacosContextBase::available_methods();
}

struct ConvertToStringVector
    : public boost::static_visitor<std::vector<std::string>> {
  auto operator()(std::string const &value) const {
    return std::vector<std::string>{value};
  }

  template <typename T, typename = std::enable_if_t<!std::is_arithmetic_v<T>>>
  std::vector<std::string> operator()(T const &value) const {
    throw std::runtime_error("Cannot convert " + Utils::demangle<T>());
  }

  template <typename T, typename = std::enable_if_t<std::is_same_v<T, int>>>
  auto operator()(T const &value) const {
    return operator()(to_str(value));
  }

  auto operator()(double const &value) const {
    return operator()(to_str(value));
  }

  auto operator()(std::vector<std::string> const &values) const {
    return values;
  }

  auto operator()(std::vector<Variant> const &values) const {
    std::vector<std::string> values_str;
    for (auto const &v : values) {
      values_str.emplace_back(boost::apply_visitor(*this, v).front());
    }
    return values_str;
  }

  template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
  auto operator()(std::vector<T> const &values) const {
    std::vector<std::string> values_str;
    for (auto const &v : values) {
      values_str.emplace_back(to_str(v));
    }
    return values_str;
  }

private:
  std::string to_str(int const value) const { return std::to_string(value); }

  std::string to_str(double const value) const {
    std::ostringstream serializer;
    serializer << std::scientific << std::setprecision(17);
    serializer << value;
    return serializer.str();
  }
};

struct GetParameterList
    : public boost::static_visitor<std::unordered_map<std::string, Variant>> {
  auto operator()(std::unordered_map<std::string, Variant> const &obj) const {
    return obj;
  }

  template <typename T>
  auto operator()(std::unordered_map<T, Variant> const &obj) const {
    // handle the case of the empty dict, which can have any key type
    return (obj.empty()) ? result_type{} : invalid(obj);
  }

  template <typename T> auto operator()(T const &obj) const {
    return invalid(obj);
  }

private:
  template <typename T> auto invalid(T const &obj) const {
    return get_value<result_type>(obj);
  }
};

std::string serialize_parameters(Variant const &pack) {
  auto const parameters = boost::apply_visitor(GetParameterList(), pack);
  if (parameters.empty()) {
    throw std::invalid_argument(
        "ScaFaCoS methods require at least 1 parameter");
  }
  auto const visitor = ConvertToStringVector();
  std::string method_params = "";
  for (auto const &kv : parameters) {
    method_params += "," + kv.first;
    for (auto const &value : boost::apply_visitor(visitor, kv.second)) {
      method_params += "," + value;
    }
  }
  return method_params.substr(1);
}

template <typename T>
boost::optional<Variant> string_to_number(std::string const &s) {
  auto deserializer = std::istringstream(s);
  T result;
  deserializer >> result;
  if (deserializer.fail() or not deserializer.eof()) {
    return {};
  }
  return Variant{result};
}

std::unordered_map<std::string, Variant>
deserialize_parameters(std::string const &parameters) {
  /*
   * ScaFaCoS parameters are serialized to a comma-separated string.
   * Key-value pairs can be split with a look ahead: when the next
   * item is a string, it is a parameter name and the current list
   * of arithmetic values belong to the current parameter name.
   * The only exception is string-valued parameters; in that case
   * the current list of arithmetic values is empty.
   */
  auto const numbers = std::string("-0123456789");
  std::unordered_map<std::string, Variant> method_params{};
  std::vector<std::string> flat_array;
  // Clang 10 false positive: https://github.com/boostorg/algorithm/issues/63
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
  boost::split(flat_array, parameters, boost::is_any_of(","));
  for (auto it = flat_array.begin(); it != flat_array.end();) {
    auto const parameter_name = *it;
    auto parameter_list = std::vector<Variant>{};
    for (++it; it != flat_array.end(); ++it) {
      if ((numbers.find(it->front()) == std::string::npos) and
          not parameter_list.empty()) {
        break;
      }
      auto result = Variant{*it};
      if (auto converted = string_to_number<int>(*it)) {
        result = Variant{*converted};
      } else if (auto converted = string_to_number<double>(*it)) {
        result = Variant{*converted};
      }
      parameter_list.emplace_back(result);
    }
    assert(not parameter_list.empty());
    if (parameter_list.size() == 1ul) {
      method_params[parameter_name] = parameter_list.front();
    } else {
      method_params[parameter_name] = Variant{std::move(parameter_list)};
    }
  }
  return method_params;
}

} // namespace Scafacos
} // namespace ScriptInterface

#endif // SCAFACOS or SCAFACOS_DIPOLES
