/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "None.hpp"

#include <utils/Vector.hpp>

#include <boost/variant.hpp>

/* This <boost/serialization/library_version_type.hpp> include guards against
 * an issue in boost::serialization from boost 1.74.0 that leads to compiler
 * error "'library_version_type' is not a member of 'boost::serialization'"
 * when including <boost/serialization/unordered_map.hpp>. More details
 * in ticket https://github.com/boostorg/serialization/issues/219
 */
#include <boost/serialization/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 == 74
#include <boost/serialization/library_version_type.hpp>
#endif

#include <boost/range/algorithm/transform.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Utils {
using Vector3b = Utils::Vector<bool, 3>;
}

namespace ScriptInterface {
class ObjectHandle;
using ObjectRef = std::shared_ptr<ObjectHandle>;
/**
 * @brief None-"literal".
 */
constexpr const None none{};

/**
 * @brief Possible types for parameters.
 *
 * The visitors and packing functions need to be adapted accordingly when
 * extending this variant with new types. For the exact details, see commit
 * <a href="https://github.com/espressomd/espresso/commit/b48ab62">b48ab62</a>.
 * The number of types is limited by macro @c BOOST_MPL_LIMIT_LIST_SIZE
 * (defaults to 20).
 */
using Variant = boost::make_recursive_variant<
    None, bool, int, std::size_t, double, std::string, ObjectRef,
    Utils::Vector3b, Utils::Vector3i, Utils::Vector2d, Utils::Vector3d,
    Utils::Vector4d, std::vector<int>, std::vector<double>,
    std::vector<boost::recursive_variant_>,
    std::unordered_map<int, boost::recursive_variant_>,
    std::unordered_map<std::string, boost::recursive_variant_>>::type;

using VariantMap = std::unordered_map<std::string, Variant>;

/**
 * @brief Make a Variant from argument.
 *
 * This is a convenience function, so that rather involved constructors from
 * boost::variant are not needed in the script interfaces.
 */
template <typename T> Variant make_variant(const T &x) { return Variant(x); }

template <typename K, typename V>
auto make_unordered_map_of_variants(std::unordered_map<K, V> const &v) {
  std::unordered_map<K, Variant> ret;
  for (auto const &it : v) {
    ret.insert({it.first, Variant(it.second)});
  }
  return ret;
}

template <typename T> auto make_vector_of_variants(std::vector<T> const &v) {
  std::vector<Variant> ret;
  for (auto const &item : v) {
    ret.emplace_back(item);
  }
  return ret;
}

namespace detail {
template <class T> struct is_type_visitor : boost::static_visitor<bool> {
  template <class U> constexpr bool operator()(const U &) const {
    return std::is_same_v<T, U>;
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
  return boost::apply_visitor(detail::is_type_visitor<T>(), v);
}

inline bool is_none(Variant const &v) { return is_type<None>(v); }
} // namespace ScriptInterface

#endif
