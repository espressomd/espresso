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

#pragma once

#include "Exception.hpp"
#include "ObjectHandle.hpp"
#include "Variant.hpp"

#include <utils/demangle.hpp>

#include <boost/algorithm/string/join.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
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
template <typename T> auto simplify_symbol(std::vector<T> const *vec) {
  auto const name_val = simplify_symbol(static_cast<T *>(nullptr));
  std::string metadata{""};
  if (vec) {
    metadata += "{.size=" + std::to_string(vec->size()) + "}";
  }
  return "std::vector<" + name_val + ">" + metadata;
}

/** @overload */
inline auto simplify_symbol(std::vector<Variant> const *vec) {
  auto value_type_name = std::string("ScriptInterface::Variant");
  std::string metadata{""};
  if (vec) {
    std::set<std::string> types = {};
    for (auto const &v : *vec) {
      types.insert(simplify_symbol_variant(v));
    }
    value_type_name += "{" + boost::algorithm::join(types, ", ") + "}";
    metadata += "{.size=" + std::to_string(vec->size()) + "}";
  }
  return "std::vector<" + value_type_name + ">" + metadata;
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
    for (auto const &variant : std::views::elements<1>(*map)) {
      types.insert(simplify_symbol_variant(variant));
    }
    value_type_name += "{" + boost::algorithm::join(types, ", ") + "}";
  }
  return "std::unordered_map<" + name_key + ", " + value_type_name + ">";
}

struct simplify_symbol_visitor {
  template <class T> std::string operator()(T const &t) const {
    return simplify_symbol(&t);
  }
};

/** @brief Simplify the demangled symbol of an object wrapped in a variant. */
inline std::string simplify_symbol_variant(Variant const &v) {
  return std::visit(simplify_symbol_visitor(), v);
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

struct simplify_symbol_containee_visitor {
  template <class T> std::string operator()(const T &) const {
    return simplify_symbol_containee(static_cast<T *>(nullptr));
  }
};

/**
 * @brief Simplify the demangled symbol of a container @c value_type wrapped
 * in a variant.
 */
inline auto simplify_symbol_containee_variant(Variant const &v) {
  return std::visit(simplify_symbol_containee_visitor(), v);
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

template <class To> struct conversion_visitor {
  template <class From> To operator()(const From &value) const {
    if constexpr (allow_conversion<To, From>::value) {
      return To(value);
    }
    throw std::bad_variant_access{};
  }
};

/**
 * @brief Implementation of get_value.
 *
 * Helper struct is needed because partial specialization of functions
 * is not allowed.
 */
template <typename T> struct get_value_helper {
  T operator()(Variant const &v) const {
    return std::visit(detail::conversion_visitor<T>(), v);
  }
};

template <class T, std::size_t N> struct vector_conversion_visitor {
  /* Catch all case -> wrong type. */
  template <typename U> Utils::Vector<T, N> operator()(U const &) const {
    throw std::bad_variant_access{};
  }

  template <typename U>
    requires allow_conversion<T, U>::value
  Utils::Vector<T, N> operator()(Utils::Vector<U, N> const &v) const {
    return Utils::Vector<T, N>(v);
  }

  template <typename U>
    requires(std::is_same_v<U, Variant> or allow_conversion<T, U>::value)
  Utils::Vector<T, N> operator()(std::vector<U> const &vector) const {
    if (vector.size() != N) {
      throw std::bad_variant_access{};
    }
    if constexpr (std::is_same_v<U, Variant>) {
      Utils::Vector<T, N> ret{};
      std::ranges::transform(vector, ret.begin(), get_value_helper<T>{});
      return ret;
    } else {
      return Utils::Vector<T, N>(vector);
    }
  }
};

/* Utils::Vector<T, N> case */
template <typename T, std::size_t N>
struct get_value_helper<Utils::Vector<T, N>> {
  Utils::Vector<T, N> operator()(Variant const &v) const {
    return std::visit(detail::vector_conversion_visitor<T, N>(), v);
  }
};

template <typename T> struct VisitorVector {
  /* Catch all case -> wrong type. */
  template <typename U> std::vector<T> operator()(U const &) const {
    throw std::bad_variant_access{};
  }

  /* Standard case, correct type */
  std::vector<T> operator()(std::vector<T> const &v) const { return v; }

  std::vector<T> operator()(std::vector<Variant> const &vv) const
    requires(not std::is_same_v<T, Variant>)
  {
    std::vector<T> ret(vv.size());

    std::ranges::transform(vv, ret.begin(), get_value_helper<T>{});

    return ret;
  }
};

/* std::vector cases */
template <typename T> struct get_value_helper<std::vector<T>> {
  std::vector<T> operator()(Variant const &v) const {
    return std::visit(VisitorVector<T>(), v);
  }
};

template <typename K, typename T> struct VisitorMap {
  /* Catch all case -> wrong type. */
  template <typename U> std::unordered_map<K, T> operator()(U const &) const {
    throw std::bad_variant_access{};
  }

  /* Standard case, correct type */
  auto operator()(std::unordered_map<K, T> const &map) const { return map; }

  auto operator()(std::unordered_map<K, Variant> const &map) const
    requires(not std::is_same_v<T, Variant>)
  {
    std::unordered_map<K, T> ret;
    for (auto const &[key, variant] : map) {
      ret.emplace(key, get_value_helper<T>{}(variant));
    }
    return ret;
  }
};

/* std::unordered_map cases */
template <typename T> struct get_value_helper<std::unordered_map<int, T>> {
  std::unordered_map<int, T> operator()(Variant const &v) const {
    return std::visit(VisitorMap<int, T>(), v);
  }
};
template <typename T>
struct get_value_helper<std::unordered_map<std::string, T>> {
  std::unordered_map<std::string, T> operator()(Variant const &v) const {
    return std::visit(VisitorMap<std::string, T>(), v);
  }
};

/* std::filesystem::path case */
template <> struct get_value_helper<std::filesystem::path> {
  auto operator()(Variant const &v) const {
    if (auto const *source = std::get_if<std::string>(&v)) {
      return std::filesystem::path(*source);
    }
    return std::get<std::filesystem::path>(v);
  }
};

/** Custom error for a conversion that fails when the value is a nullptr. */
class bad_get_nullptr : public std::bad_variant_access {};

/* This allows direct retrieval of a shared_ptr to the object from
 * an ObjectRef variant. If the type is a derived type, the type is
 * also checked.
 */
template <typename T>
  requires(std::is_base_of_v<ObjectHandle, T>)
struct get_value_helper<std::shared_ptr<T>> {
  std::shared_ptr<T> operator()(Variant const &v) const {
    auto so_ptr = std::get<ObjectRef>(v);
    if (!so_ptr) {
      throw bad_get_nullptr{};
    }

    if (auto t_ptr = std::dynamic_pointer_cast<T>(so_ptr)) {
      return t_ptr;
    }

    throw std::bad_variant_access{};
  }
};

/**
 * @brief Re-throw a @c std::bad_variant_access wrapped in an @ref Exception.
 * Write a custom error message for invalid conversions due to type mismatch
 * and due to nullptr values, possibly with context information if the variant
 * is a container.
 * @tparam T     Which type the variant was supposed to convert to
 */
template <typename T>
inline void handle_bad_get(Variant const &v, std::string const &name) {
  auto const container_name = demangle::simplify_symbol_variant(v);
  auto const containee_name = demangle::simplify_symbol_containee_variant(v);
  auto const expected_containee_name =
      demangle::simplify_symbol_containee(static_cast<T *>(nullptr));
  auto const from_container = !containee_name.empty();
  auto const to_container = !expected_containee_name.empty();
  auto what = "Provided argument of type '" + container_name + "'";
  if (not name.empty()) {
    what += " for parameter '" + name + "'";
  }
  try {
    throw;
  } catch (bad_get_nullptr const &) {
    auto const item_error = (to_container) ? " contains a value that" : "";
    throw Exception(what + item_error + " is a null pointer");
  } catch (std::bad_variant_access const &) {
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

template <typename T> T get_value(Variant const &v, std::string const &name) {
  try {
    return detail::get_value_helper<T>{}(v);
  } catch (...) {
    detail::handle_bad_get<T>(v, name);
    throw;
  }
}

} // namespace detail

/**
 * @brief Extract value of specific type T from a Variant.
 *
 * This is a wrapper around std::get that allows us to
 * customize the behavior for different types. This is
 * needed e.g. to deal with containers whose elements
 * have mixed types that are implicitly convertible
 * to a requested type.
 */
template <typename T> T get_value(Variant const &v) {
  return detail::get_value<T>(v, "");
}

/**
 * @brief Get a value from a VariantMap by name, or throw
 *        if it does not exist or is not convertible to
 *        the target type.
 */
template <typename T>
T get_value(VariantMap const &vals, std::string const &name) {
  if (not vals.contains(name)) {
    throw Exception("Parameter '" + name + "' is missing.");
  }
  return detail::get_value<T>(vals.at(name), name);
}

/**
 * @brief Get a value from a VariantMap by name, or return a default
 *        value if it does not exist.
 */
template <typename T>
T get_value_or(VariantMap const &vals, std::string const &name,
               T const &default_) {
  if (vals.contains(name)) {
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
