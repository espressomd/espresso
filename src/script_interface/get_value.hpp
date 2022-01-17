/*
 * Copyright (C) 2016-2020 The ESPResSo project
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

#include <boost/range/algorithm/transform.hpp>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace detail {
struct type_label_visitor : boost::static_visitor<std::string> {
  template <class T> std::string operator()(const T &) const {
    return Utils::demangle<T>();
  }
};

inline std::string type_label(const Variant &v) {
  return boost::apply_visitor(type_label_visitor{}, v);
}

/*
 * Allows
 * T -> T,
 * floating point -> floating point and
 * integral -> floating point
 */
template <class To, class From>
using allow_conversion =
    std::integral_constant<bool, std::is_same<To, From>::value ||
                                     (std::is_convertible<To, From>::value &&
                                      std::is_floating_point<To>::value &&
                                      std::is_arithmetic<From>::value)>;

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
    return boost::apply_visitor(detail::conversion_visitor<T>{}, v);
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

  Utils::Vector<T, N>
  operator()(std::vector<T, std::allocator<T>> const &v) const {
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
    return boost::apply_visitor(detail::vector_conversion_visitor<T, N>{}, v);
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
  std::vector<T> operator()(std::vector<Variant> const &vv) const {
    std::vector<T> ret(vv.size());

    boost::transform(vv, ret.begin(),
                     [](const Variant &v) { return get_value_helper<T>{}(v); });

    return ret;
  }
};

/* std::vector cases
 * We implicitly transform an empty vector<Variant> into an empty vector<T>. */
template <> struct get_value_helper<std::vector<int>, void> {
  std::vector<int> operator()(Variant const &v) const {
    return boost::apply_visitor(GetVectorOrEmpty<int>{}, v);
  }
};

template <> struct get_value_helper<std::vector<double>, void> {
  std::vector<double> operator()(Variant const &v) const {
    return boost::apply_visitor(GetVectorOrEmpty<double>{}, v);
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
};

/* std::unordered_map cases */
template <> struct get_value_helper<std::unordered_map<int, Variant>, void> {
  std::unordered_map<int, Variant> operator()(Variant const &v) const {
    return boost::apply_visitor(GetMapOrEmpty<int, Variant>{}, v);
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
    typename std::enable_if<std::is_base_of<ObjectHandle, T>::value,
                            void>::type> {
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
 * is in a container.
 * @tparam T     Type which the variant was supposed to convert to
 * @tparam U     Type of the container the variant is contained in
 */
template <typename T, typename U = void>
inline void handle_bad_get(Variant const &v) {
  using Utils::demangle;
  auto const what = "Provided argument of type " + type_label(v);
  auto const context =
      (std::is_void<U>::value)
          ? ""
          : " (raised during the creation of a " + demangle<U>() + ")";
  try {
    throw;
  } catch (bad_get_nullptr const &) {
    throw Exception(what + " is a null pointer" + context);
  } catch (boost::bad_get const &) {
    throw Exception(what + " is not convertible to " + demangle<T>() + context);
  }
}

} // namespace detail

/**
 * @brief Extract value of specific type T from a Variant.
 *
 * This is a wrapper around boost::get that allows us to
 * customize the behavior for different types. This is
 * needed e.g. to deal with Vector types that are not
 * explicitly contained in Variant, but can easily
 * be converted.
 */
template <typename T> T get_value(Variant const &v) {
  try {
    return detail::get_value_helper<T>{}(v);
  } catch (...) {
    detail::handle_bad_get<T>(v);
    throw;
  }
}

template <typename K, typename V>
std::unordered_map<K, V> get_map(std::unordered_map<K, Variant> const &v) {
  std::unordered_map<K, V> ret;
  auto it = v.begin();
  try {
    for (; it != v.end(); ++it) {
      ret.insert({it->first, detail::get_value_helper<V>{}(it->second)});
    }
  } catch (...) {
    detail::handle_bad_get<V, std::unordered_map<K, V>>(v);
    throw;
  }
  return ret;
}

template <typename K, typename V>
std::unordered_map<K, Variant> make_map(std::unordered_map<K, V> const &v) {
  std::unordered_map<K, Variant> ret;
  for (auto const &it : v) {
    ret.insert({it.first, Variant(it.second)});
  }
  return ret;
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
                                         ArgNames &&... args) {
  return std::make_shared<T>(
      get_value<Types>(vals, std::forward<ArgNames>(args))...);
}

template <typename T>
void set_from_args(T &dst, VariantMap const &vals, const char *name) {
  dst = get_value<T>(vals, name);
}
} /* namespace ScriptInterface */

#endif
