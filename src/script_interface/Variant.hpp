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

#include "None.hpp"
#include "utils/Vector.hpp"

#include <boost/variant.hpp>

#include <boost/range/algorithm/transform.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include <unordered_map>

namespace ScriptInterface {
    enum class CreationPolicy { LOCAL, GLOBAL };

    class ObjectHandle;
using ObjectRef = std::shared_ptr<ObjectHandle>;
/**
 * @brief None-"literal".
 */
constexpr const None none{};

/**
 * @brief Possible types for parameters.
 */
using Variant = boost::make_recursive_variant<
    None, bool, int, size_t, double, std::string, std::vector<int>,
    std::vector<double>, ObjectRef, std::vector<boost::recursive_variant_>,
    Utils::Vector2d, Utils::Vector3d, Utils::Vector4d>::type;

using VariantMap = std::unordered_map<std::string, Variant>;

/**
 * @brief Make a Variant from argument.
 *
 * This is a convenience function, so that rather involved constructors from
 * boost::variant are not needed in the script interfaces.
 */
template <typename T> Variant make_variant(const T &x) { return Variant(x); }

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

template <class D, class V, class R>
struct recursive_visitor : boost::static_visitor<R> {
  std::enable_if_t<not std::is_void<R>::value, R>
  operator()(const std::vector<V> &vec) const {
    std::vector<R> ret(vec.size());

    boost::transform(vec, ret.begin(),
                     [visitor = static_cast<const D *>(this)](const V &v) {
                       return boost::apply_visitor(*visitor, v);
                     });

    return ret;
  }
};
} /* namespace ScriptInterface */

#endif
