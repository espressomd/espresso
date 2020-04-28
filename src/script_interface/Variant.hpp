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
#ifndef SCRIPT_INTERFACE_VARIANT_HPP
#define SCRIPT_INTERFACE_VARIANT_HPP

#include <boost/variant.hpp>

#include "None.hpp"
#include "utils/AutoObjectId.hpp"
#include "utils/Vector.hpp"
#include "utils/serialization/unordered_map.hpp"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

namespace ScriptInterface {
class ScriptInterfaceBase;
using ObjectId = Utils::ObjectId<ScriptInterfaceBase>;
/**
 * @brief None-"literal".
 */
constexpr const None none{};

/**
 * @brief Possible types for parameters.
 */
using Variant = boost::make_recursive_variant<
    None, bool, int, size_t, double, std::string, std::vector<int>,
    std::vector<double>, ObjectId, std::vector<boost::recursive_variant_>,
    Utils::Vector2d, Utils::Vector3d, Utils::Vector4d>::type;

using VariantMap = std::unordered_map<std::string, Variant>;

namespace detail {
template <class T> struct is_type_visitor : boost::static_visitor<bool> {
  template <class U> constexpr bool operator()(const U &) const {
    return std::is_same<T, U>::value;
  }
};
} // namespace detail

/**
 * @brief Check is a Variant holds a specific type.
 *
 * @tparam T type to check for
 * @param v Variant to check in
 * @return true, if v holds a T.
 */
template <class T> bool is_type(Variant const &v) {
  return boost::apply_visitor(detail::is_type_visitor<T>{}, v);
}

inline bool is_none(Variant const &v) { return is_type<None>(v); }

template <typename T, typename U>
std::vector<Variant> pack_pair(const std::pair<T, U> &pair) {
  return {{pair.first, pair.second}};
}

template <typename T, typename U>
const std::pair<T, U> unpack_pair(const std::vector<Variant> &v) {
  return {boost::get<T>(v.at(0)), boost::get<U>(v.at(1))};
}

/**
 * @brief Pack a map into a vector of Variants
 *        by serializing the key-value pairs.
 *
 */
template <typename K, typename V>
std::vector<Variant> pack_map(const std::unordered_map<K, V> &map) {
  std::vector<Variant> ret(map.size());

  std::transform(map.begin(), map.end(), ret.begin(),
                 [](const std::pair<K, V> &p) { return pack_pair(p); });

  return ret;
}

template <typename K, typename V>
std::unordered_map<K, V> unpack_map(const std::vector<Variant> &v) {
  std::unordered_map<K, V> ret;

  for (auto const &pair : v) {
    ret.insert(unpack_pair<K, V>(boost::get<const std::vector<Variant>>(pair)));
  }

  return ret;
}
} /* namespace ScriptInterface */

#endif
