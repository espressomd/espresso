/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_GET_VALUE_HPP
#define SCRIPT_INTERFACE_GET_VALUE_HPP

#include "Exception.hpp"
#include "ObjectHandle.hpp"
#include "Variant.hpp"

#include <utils/demangle.hpp>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/algorithm/transform.hpp>

#include <cstddef>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace detail {

/**
 * @brief Convert a demangled symbol into a human-readable name that omits
 * container allocators, key hashes and implementation-specific namespaces.
 * When the data type involves the @ref Variant type, it is recursively
 * replaced by the string "Variant{type1,type2,...}" based on the actual
 * contents of the variant.
 */
namespace demangle {

inline std::string simplify_symbol_variant(Variant const &v);

/** @brief Simplify the demangled symbol of an object. */
template <typename T> auto simplify_symbol(T const *) {
  auto constexpr is_string = std::is_same_v<T, std::string>;
  auto const symbol_for_variant = Utils::demangle<Variant>();
  auto const name_for_variant = std::string("ScriptInterface::Variant");
  auto name = (is_string) ? std::string{"std::string"} : Utils::demangle<T>();
  for (std::string::size_type pos{};
       (pos = name.find(symbol_for_variant, pos)) != name.npos;
       pos += name_for_variant.length()) {
    name.replace(pos, symbol_for_variant.length(), name_for_variant);
  }
  return name;
}

/** @overload */
template <typename T, std::size_t N>
auto simplify_symbol(Utils::Vector<T, N> const *) {
  auto const name_val = simplify_symbol(static_cast<T *>(nullptr));
  return "Utils::Vector<" + Utils::demangle<T>() + ", " + std::to_string(N) +
         ">";
}

/** @overload */
template <typename T> auto simplify_symbol(std::vector<T> const *) {
  auto const name_val = simplify_symbol(static_cast<T *>(nullptr));
  return "std::vector<" + name_val + ">";
}

/** @overload */
inline auto simplify_symbol(std::vector<Variant> const *vec) {
  auto value_type_name = std::string("ScriptInterface::Variant");
  if (vec) {
    std::set<std::string> types = {};
    for (auto const &v : *vec) {
      types.insert(simplify_symbol_variant(v));
    }
    value_type_name += "{" + boost::algorithm::join(types, ", ") + "}";
  }
  return "std::vector<" + value_type_name + ">";
}

/** @overload */
template <typename K, typename V>
auto simplify_symbol(std::unordered_map<K, V> const *) {
  auto const name_key = simplify_symbol(static_cast<K *>(nullptr));
  auto const name_val = simplify_symbol(static_cast<V *>(nullptr));
  return "std::unordered_map<" + name_key + ", " + name_val + ">";
}

/** @overload */
template <typename K>
auto simplify_symbol(std::unordered_map<K, Variant> const *map) {
  auto const name_key = simplify_symbol(static_cast<K *>(nullptr));
  auto value_type_name = std::string("ScriptInterface::Variant");
  if (map) {
    std::set<std::string> types = {};
    for (auto const &kv : *map) {
      types.insert(simplify_symbol_variant(kv.second));
    }
    value_type_name += "{" + boost::algorithm::join(types, ", ") + "}";
  }
  return "std::unordered_map<" + name_key + ", " + value_type_name + ">";
}

struct simplify_symbol_visitor : boost::static_visitor<std::string> {
  template <class T> std::string operator()(T const &t) const {
    return simplify_symbol(&t);
  }
};

/** @brief Simplify the demangled symbol of an object wrapped in a variant. */
inline std::string simplify_symbol_variant(Variant const &v) {
  return boost::apply_visitor(simplify_symbol_visitor(), v);
}

/** @brief Simplify the demangled symbol of a container @c value_type. */
template <typename T> auto simplify_symbol_containee(T const *) {
  return std::string("");
}

/** @overload */
template <typename T> auto simplify_symbol_containee(std::vector<T> const *) {
  auto const name_val = simplify_symbol(static_cast<T *>(nullptr));
  return name_val;
}

/** @overload */
template <typename K, typename V>
auto simplify_symbol_containee(std::unordered_map<K, V> const *) {
  auto const name_key = simplify_symbol(static_cast<K *>(nullptr));
  auto const name_val = simplify_symbol(static_cast<V *>(nullptr));
  return name_key + "' or '" + name_val;
}

struct simplify_symbol_containee_visitor : boost::static_visitor<std::string> {
  template <class T> std::string operator()(const T &) const {
    return simplify_symbol_containee(static_cast<T *>(nullptr));
  }
};

/**
 * @brief Simplify the demangled symbol of a container @c value_type wrapped
 * in a variant.
 */
inline auto simplify_symbol_containee_variant(Variant const &v) {
  return boost::apply_visitor(simplify_symbol_containee_visitor(), v);
}

} // namespace demangle

/*
 * Allows
 * T -> T,
 * floating point -> floating point and
 * integral -> floating point
 */
template <class To, class From>
using allow_conversion =
    std::integral_constant<bool, std::is_same_v<To, From> ||
                                     (std::is_convertible_v<To, From> &&
                                      std::is_floating_point_v<To> &&
                                      std::is_arithmetic_v<From>)>;

template <class To> struct conversion_visitor : boost::static_visitor<To> {
  template <class From>
  std::enable_if_t<allow_conversion<To, From>::value, To>
  operator()(const From &value) const {
    return To(value);
  }

  template <class From>
  std::enable_if_t<!allow_conversion<To, From>::value, To>
  operator()(const From &) const {
    throw boost::bad_get{};
  }
};

/**
 * @brief Implementation of get_value.
 *
 * Helper struct is needed because partial specialization of functions
 * is not allowed.
 */
template <typename T, typename = void> struct get_value_helper {
  T operator()(Variant const &v) const {
    return boost::apply_visitor(detail::conversion_visitor<T>(), v);
  }
};

template <class T, std::size_t N>
struct vector_conversion_visitor : boost::static_visitor<Utils::Vector<T, N>> {
  Utils::Vector<T, N> operator()(Utils::Vector<T, N> const &v) const {
    return v;
  }

  /* We try to unpack variant vectors and check if they
   * are convertible element by element. */
  auto operator()(std::vector<Variant> const &vv) const {
    if (N != vv.size()) {
      throw boost::bad_get{};
    }

    Utils::Vector<T, N> ret;
    boost::transform(vv, ret.begin(),
                     [](const Variant &v) { return get_value_helper<T>{}(v); });

    return ret;
  }

  template <typename U>
  std::enable_if_t<allow_conversion<T, U>::value, Utils::Vector<T, N>>
  operator()(std::vector<U, std::allocator<U>> const &v) const {
    if (N != v.size()) {
      throw boost::bad_get{};
    }
    return Utils::Vector<T, N>(v);
  }

  template <typename U> Utils::Vector<T, N> operator()(U const &) const {
    throw boost::bad_get{};
  }
};

/* Utils::Vector<T, N> case */
template <typename T, std::size_t N>
struct get_value_helper<Utils::Vector<T, N>> {
  Utils::Vector<T, N> operator()(Variant const &v) const {
    return boost::apply_visitor(detail::vector_conversion_visitor<T, N>(), v);
  }
};

template <typename T>
struct GetVectorOrEmpty : boost::static_visitor<std::vector<T>> {
  /* Catch all case -> wrong type. */
  template <typename U> std::vector<T> operator()(U const &) const {
    throw boost::bad_get{};
  }

  /* Standard case, correct type */
  std::vector<T> operator()(std::vector<T> const &v) const { return v; }

  template <typename V = T,
            std::enable_if_t<!std::is_same_v<V, Variant>, bool> = true>
  std::vector<T> operator()(std::vector<Variant> const &vv) const {
    std::vector<T> ret(vv.size());

    boost::transform(vv, ret.begin(),
                     [](const Variant &v) { return get_value_helper<T>{}(v); });

    return ret;
  }
};

/* std::vector cases */
template <typename T> struct get_value_helper<std::vector<T>, void> {
  std::vector<T> operator()(Variant const &v) const {
    return boost::apply_visitor(GetVectorOrEmpty<T>(), v);
  }
};

template <typename K, typename T>
struct GetMapOrEmpty : boost::static_visitor<std::unordered_map<K, T>> {
  /* Catch all case -> wrong type. */
  template <typename U> std::unordered_map<K, T> operator()(U const &) const {
    throw boost::bad_get{};
  }

  /* Standard case, correct type */
  std::unordered_map<K, T> operator()(std::unordered_map<K, T> const &v) const {
    return v;
  }

  template <typename V = T,
            std::enable_if_t<!std::is_same_v<V, Variant>, bool> = true>
  std::unordered_map<K, T>
  operator()(std::unordered_map<K, Variant> const &v) const {
    std::unordered_map<K, T> ret;
    for (auto it = v.begin(); it != v.end(); ++it) {
      ret.insert({it->first, get_value_helper<T>{}(it->second)});
    }
    return ret;
  }
};

/* std::unordered_map cases */
template <typename T>
struct get_value_helper<std::unordered_map<int, T>, void> {
  std::unordered_map<int, T> operator()(Variant const &v) const {
    return boost::apply_visitor(GetMapOrEmpty<int, T>(), v);
  }
};
template <typename T>
struct get_value_helper<std::unordered_map<std::string, T>, void> {
  std::unordered_map<std::string, T> operator()(Variant const &v) const {
    return boost::apply_visitor(GetMapOrEmpty<std::string, T>(), v);
  }
};

/** Custom error for a conversion that fails when the value is a nullptr. */
class bad_get_nullptr : public boost::bad_get {};

/* This allows direct retrieval of a shared_ptr to the object from
 * an ObjectRef variant. If the type is a derived type, the type is
 * also checked.
 */
template <typename T>
struct get_value_helper<
    std::shared_ptr<T>,
    typename std::enable_if_t<std::is_base_of_v<ObjectHandle, T>, void>> {
  std::shared_ptr<T> operator()(Variant const &v) const {
    auto so_ptr = boost::get<ObjectRef>(v);
    if (!so_ptr) {
      throw bad_get_nullptr{};
    }

    auto t_ptr = std::dynamic_pointer_cast<T>(so_ptr);

    if (t_ptr) {
      return t_ptr;
    }

    throw boost::bad_get{};
  }
};

/**
 * @brief Re-throw a @c boost::bad_get exception wrapped in an @ref Exception.
 * Write a custom error message for invalid conversions due to type mismatch
 * and due to nullptr values, possibly with context information if the variant
 * is a container.
 * @tparam T     Which type the variant was supposed to convert to
 */
template <typename T> inline void handle_bad_get(Variant const &v) {
  auto const container_name = demangle::simplify_symbol_variant(v);
  auto const containee_name = demangle::simplify_symbol_containee_variant(v);
  auto const expected_containee_name =
      demangle::simplify_symbol_containee(static_cast<T *>(nullptr));
  auto const from_container = !containee_name.empty();
  auto const to_container = !expected_containee_name.empty();
  auto const what = "Provided argument of type '" + container_name + "'";
  try {
    throw;
  } catch (bad_get_nullptr const &) {
    auto const item_error = (to_container) ? " contains a value that" : "";
    throw Exception(what + item_error + " is a null pointer");
  } catch (boost::bad_get const &) {
    auto const non_convertible = std::string(" is not convertible to ");
    auto item_error = std::string("");
    if (from_container and to_container) {
      item_error += " because it contains a value that";
      item_error += non_convertible + "'" + expected_containee_name + "'";
    }
    auto const target = demangle::simplify_symbol(static_cast<T *>(nullptr));
    throw Exception(what + non_convertible + "'" + target + "'" + item_error);
  }
}

} // namespace detail

/**
 * @brief Extract value of specific type T from a Variant.
 *
 * This is a wrapper around boost::get that allows us to
 * customize the behavior for different types. This is
 * needed e.g. to deal with containers whose elements
 * have mixed types that are implicitly convertible
 * to a requested type.
 */
template <typename T> T get_value(Variant const &v) {
  try {
    return detail::get_value_helper<T>{}(v);
  } catch (...) {
    detail::handle_bad_get<T>(v);
    throw;
  }
}

/**
 * @brief Get a value from a VariantMap by name, or throw
 *        if it does not exist or is not convertible to
 *        the target type.
 *
 */
template <typename T>
T get_value(VariantMap const &vals, std::string const &name) {
  try {
    return get_value<T>(vals.at(name));
  } catch (std::out_of_range const &) {
    throw Exception("Parameter '" + name + "' is missing.");
  }
}

/**
 * @brief Get a value from a VariantMap by name, or return a default
 *        value if it does not exist.
 */
template <typename T>
T get_value_or(VariantMap const &vals, std::string const &name,
               T const &default_) {
  if (vals.count(name)) {
    return get_value<T>(vals.at(name));
  }
  return default_;
}

/**
 * @brief Make a new std::shared_ptr<T> with arguments extracted from a
 * VariantMap.
 */
template <typename T, typename... Types, typename... ArgNames>
std::shared_ptr<T> make_shared_from_args(VariantMap const &vals,
                                         ArgNames &&...args) {
  return std::make_shared<T>(
      get_value<Types>(vals, std::forward<ArgNames>(args))...);
}

template <typename T>
void set_from_args(T &dst, VariantMap const &vals, const char *name) {
  dst = get_value<T>(vals, name);
}
} // namespace ScriptInterface

#endif
