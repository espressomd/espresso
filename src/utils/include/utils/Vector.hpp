/*
 * Copyright (C) 2014-2026 The ESPResSo project
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

/**
 * @file
 *
 * @brief Vector implementation and trait types
 * for boost qvm interoperability.
 */

#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/vec_traits.hpp>

#include "utils/Array.hpp"
#include "utils/attributes.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <ranges>
#include <span>
#include <type_traits>
#include <vector>

namespace Utils {

template <typename T, std::size_t N> class Vector : public Array<T, N> {
  using Base = Array<T, N>;

public:
  using Base::at;
  using Base::Base;
  using Base::operator[];
  using Base::back;
  using Base::begin;
  using Base::cbegin;
  using Base::cend;
  using Base::data;
  using Base::empty;
  using Base::end;
  using Base::fill;
  using Base::front;
  using Base::max_size;
  using Base::size;

  template <class U> struct is_vector : std::false_type {};
  template <class U, std::size_t Np>
  struct is_vector<Vector<U, Np>> : std::true_type {};

  Vector() noexcept = default;
  Vector(Vector const &) = default;
  Vector &operator=(Vector const &) = default;

  void swap(Vector &rhs) { std::ranges::swap_ranges(*this, rhs); }

private:
  constexpr void copy_init(T const *values) noexcept {
    for (std::size_t i{0}; i != N; ++i) {
      (*this)[i] = values[i];
    }
  }

public:
  // range-based ctor that excludes Vector<T,N> to avoid ambiguous calls with
  // the copy ctor, move ctor and cast operator; std::ranges::input_range is
  // not used due to conflicts with move assignment in recursive variant types
  // and T[N] must be excluded to avoid shadowing the noexcept ctor
  template <class Range>
    requires(not is_vector<std::remove_cvref_t<Range>>::value and
             not std::is_same_v<std::remove_cvref_t<Range>, T[N]>)
  explicit constexpr Vector(Range &&rng)
      : Vector(std::begin(rng), std::end(rng)) {}

#if __cpp_lib_containers_ranges
  template <std::ranges::input_range Range>
  Vector(std::from_range_t, Range &&rng)
      : Vector(std::begin(rng), std::end(rng)) {}
#endif

  explicit constexpr Vector(T const (&v)[N]) noexcept : Base() {
    if constexpr (N != 0) {
      copy_init(std::cbegin(v));
    }
  }

  constexpr Vector(std::initializer_list<T> v) : Base() {
    if (N != v.size()) {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }
    if constexpr (N != 0) {
      copy_init(v.begin());
    }
  }

  template <typename InputIterator>
  Vector(InputIterator first, InputIterator last) : Base() {
    if (std::distance(first, last) == N) {
      std::copy_n(first, N, begin());
    } else {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }
  }

  /** @brief Create a vector that has all entries set to the same value. */
  DEVICE_QUALIFIER static constexpr Vector<T, N>
  broadcast(typename Base::value_type const &value) noexcept {
    Vector<T, N> ret;
    for (std::size_t i = 0u; i != N; ++i) {
      ret[i] = value;
    }
    return ret;
  }

  std::vector<T> as_vector() const { return std::vector<T>(begin(), end()); }

  operator std::vector<T>() const { return as_vector(); }

  constexpr std::span<T, N> as_span() const {
    return std::span<T, N>(const_cast<T *>(begin()), size());
  }

  constexpr operator std::span<T, N>() const { return as_span(); }

  template <class U> explicit operator Vector<U, N>() const {
    Vector<U, N> ret;

    std::ranges::transform(*this, ret.begin(),
                           [](T const &e) { return static_cast<U>(e); });

    return ret;
  }

  constexpr T norm2() const { return (*this) * (*this); }
  T norm() const { return std::sqrt(norm2()); }

  /*
   * @brief Normalize the vector.
   *
   * Normalize the vector by its length,
   * if not zero, otherwise the vector is unchanged.
   */
  Vector &normalize() {
    auto const l = norm();
    if (l != T(0)) {
      *this /= l;
    }

    return *this;
  }

  /*
   * @brief Return a normalized copy of the vector.
   *
   * Normalize the vector by its length,
   * if not zero, otherwise the vector is unchanged.
   */
  Vector normalized() const {
    auto const l = norm();
    return (l != T(0)) ? (*this) / l : *this;
  }
};

template <class T> using Vector3 = Vector<T, 3>;

template <std::size_t N> using VectorXd = Vector<double, N>;
using Vector2d = VectorXd<2>;
using Vector3d = VectorXd<3>;
using Vector4d = VectorXd<4>;
using Vector6d = VectorXd<6>;
using Vector9d = VectorXd<9>;

template <std::size_t N> using VectorXf = Vector<float, N>;
using Vector3f = VectorXf<3>;

template <std::size_t N> using VectorXi = Vector<int, N>;
using Vector3i = VectorXi<3>;

namespace detail {
template <std::size_t N, typename T, typename U, typename Op>
auto binary_op(Vector<T, N> const &a, Vector<U, N> const &b, Op op) {
  // we must use the non-range version std::transform for Vector<T, N> because:
  // GCC 12 cannot inspect -Wmaybe-uninitialized through std::ranges::transform
  // Clang/Xcode libc++ cannot select the binary form of std::ranges::transform
  using R = decltype(std::declval<T>() + std::declval<U>());
  Vector<R, N> ret;

  // NOLINTNEXTLINE(modernize-use-ranges)
  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 op);

  return ret;
}

template <std::size_t N, typename T, typename Op>
constexpr bool all_of(Vector<T, N> const &a, Vector<T, N> const &b, Op op) {
  for (std::size_t i = 0u; i < N; ++i) {
    if (not op(a[i], b[i])) {
      return false;
    }
  }
  return true;
}
} // namespace detail

template <std::size_t N, typename T>
constexpr bool operator<(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less<T>());
}

template <std::size_t N, typename T>
constexpr bool operator>(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater<T>());
}

template <std::size_t N, typename T>
constexpr bool operator<=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less_equal<T>());
}

template <std::size_t N, typename T>
constexpr bool operator>=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater_equal<T>());
}

template <std::size_t N, typename T>
constexpr bool operator==(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::equal_to<T>());
}

template <std::size_t N, typename T>
constexpr bool operator!=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return not(a == b);
}

template <std::size_t N, typename T, typename U>
auto operator+(Vector<T, N> const &a, Vector<U, N> const &b) {
  return detail::binary_op(a, b, std::plus<>());
}

template <std::size_t N, typename T>
ESPRESSO_ATTR_ALWAYS_INLINE inline auto &operator+=(Vector<T, N> &a,
                                                    Vector<T, N> const &b) {
  std::ranges::transform(a, b, std::begin(a), std::plus<T>());
  return a;
}

template <std::size_t N, typename T, typename U>
auto operator-(Vector<T, N> const &a, Vector<U, N> const &b) {
  return detail::binary_op(a, b, std::minus<>());
}

template <std::size_t N, typename T>
Vector<T, N> operator-(Vector<T, N> const &a) {
  Vector<T, N> ret;
  std::ranges::transform(a, std::begin(ret), std::negate<T>());
  return ret;
}

template <std::size_t N, typename T>
ESPRESSO_ATTR_ALWAYS_INLINE inline Vector<T, N> &
operator-=(Vector<T, N> &a, Vector<T, N> const &b) {
  std::ranges::transform(a, b, std::begin(a), std::minus<T>());
  return a;
}

/* Scalar multiplication */
template <std::size_t N, typename T, class U>
  requires(std::is_arithmetic_v<U>)
constexpr auto operator*(U const &a, Vector<T, N> const &b) {
  using R = decltype(a * std::declval<T>());
  Vector<R, N> ret;
  std::ranges::transform(b, std::begin(ret), [a](T const &v) { return a * v; });
  return ret;
}

template <std::size_t N, typename T, class U>
  requires(std::is_arithmetic_v<U>)
constexpr auto operator*(Vector<T, N> const &a, U const &b) {
  using R = decltype(std::declval<T>() * b);
  Vector<R, N> ret;
  std::ranges::transform(a, std::begin(ret), [b](T const &v) { return b * v; });
  return ret;
}

template <std::size_t N, typename T>
auto &operator*=(Vector<T, N> &b, T const &a) {
  std::ranges::transform(b, std::begin(b), [a](T const &v) { return a * v; });
  return b;
}

/* Scalar division */
template <std::size_t N, typename T, class U>
auto operator/(Vector<T, N> const &a, U const &b) {
  using R = decltype(std::declval<T>() / b);
  Vector<R, N> ret;
  std::ranges::transform(a, std::begin(ret), [b](T const &v) { return v / b; });
  return ret;
}

template <std::size_t N, typename T, class U>
auto operator/(U const &a, Vector<T, N> const &b) {
  using R = decltype(a / std::declval<T>());
  Vector<R, N> ret;
  std::ranges::transform(b, std::begin(ret), [a](T const &v) { return a / v; });
  return ret;
}

template <std::size_t N, typename T>
auto &operator/=(Vector<T, N> &a, T const &b) {
  std::ranges::transform(a, std::begin(a), [b](T const &v) { return v / b; });
  return a;
}

namespace detail {
template <class T> using is_vector = Vector<int, 1>::is_vector<T>;
} // namespace detail

/* Scalar product */
template <std::size_t N, typename T, class U>
  requires(not(detail::is_vector<T>::value or detail::is_vector<U>::value))
auto constexpr operator*(Vector<T, N> const &a, Vector<U, N> const &b) {
  using R = decltype(std::declval<T>() * std::declval<U>());
  // std::inner_product isn't always inlined on Intel CPUs even with -O3,
  // but a for loop can be inlined
  R acc{};
  for (std::size_t i = 0u; i < N; ++i) {
    acc += a[i] * b[i];
  }
  return acc;
}

template <std::size_t N, typename T, class U>
  requires(std::is_integral_v<T> and std::is_integral_v<U>)
auto operator%(Vector<T, N> const &a, Vector<U, N> const &b) {
  using R = decltype(std::declval<T>() % std::declval<U>());
  Vector<R, N> ret;
  std::ranges::transform(a, b, std::begin(ret), std::modulus<>());
  return ret;
}

/* Componentwise square root */
template <std::size_t N, typename T> auto sqrt(Vector<T, N> const &a) {
  using std::sqrt;
  using R = decltype(sqrt(std::declval<T>()));
  Vector<R, N> ret;
  std::ranges::transform(a, ret.begin(), [](T const &v) { return sqrt(v); });
  return ret;
}

template <class T>
Vector<T, 3> vector_product(Vector<T, 3> const &a, Vector<T, 3> const &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// Product of array elements.
template <class T, std::size_t N> T product(Vector<T, N> const &v) {
  return std::accumulate(v.cbegin(), v.cend(), T{1}, std::multiplies<T>());
}

template <class T, class U, std::size_t N>
auto hadamard_product(Vector<T, N> const &a, Vector<U, N> const &b) {
  using R = decltype(std::declval<T>() * std::declval<U>());
  auto constexpr proj = std::identity{}; // required by Clang/Xcode libc++
  Vector<R, N> ret;
  std::ranges::transform(a, b, ret.begin(), std::multiplies<>(), proj, proj);
  return ret;
}

// specialization for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <typename T, typename U>
  requires(not(detail::is_vector<T>::value and detail::is_vector<U>::value))
auto hadamard_product(T const &a, U const &b) {
  return a * b;
}

template <class T, class U, std::size_t N>
auto hadamard_division(Vector<T, N> const &a, Vector<U, N> const &b) {
  using R = decltype(std::declval<T>() / std::declval<U>());
  auto constexpr proj = std::identity{}; // required by Clang/Xcode libc++
  Vector<R, N> ret;
  std::ranges::transform(a, b, std::begin(ret), std::divides<>(), proj, proj);
  return ret;
}

// specialization for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <typename T, typename U>
  requires(not(detail::is_vector<T>::value and detail::is_vector<U>::value))
auto hadamard_division(T const &a, U const &b) {
  return a / b;
}

template <typename T> Vector<T, 3> unit_vector(unsigned int i) {
  if (i == 0u)
    return {T{1}, T{0}, T{0}};
  if (i == 1u)
    return {T{0}, T{1}, T{0}};
  if (i == 2u)
    return {T{0}, T{0}, T{1}};
  throw std::domain_error("coordinate out of range");
}

/**
 * @brief Meta function to turn a Vector<T, 1> into T.
 */
template <typename T> struct decay_to_scalar {};
template <typename T, std::size_t N> struct decay_to_scalar<Vector<T, N>> {
  using type = Vector<T, N>;
};

template <typename T> struct decay_to_scalar<Vector<T, 1>> {
  using type = T;
};

template <std::size_t I, class T, std::size_t N>
T &get(Vector<T, N> &a) noexcept {
  return a[I];
}

template <std::size_t I, class T, std::size_t N>
T const &get(Vector<T, N> const &a) noexcept {
  return a[I];
}

} // namespace Utils

template <std::size_t I, class T, std::size_t N>
struct std::tuple_element<I, Utils::Vector<T, N>> {
  static_assert(I < N, "Utils::Vector index must be in range");
  using type = T;
};

template <class T, std::size_t N>
struct std::tuple_size<Utils::Vector<T, N>>
    : std::integral_constant<std::size_t, N> {};

namespace boost::qvm {

template <class T, std::size_t N> struct vec_traits<::Utils::Vector<T, N>> {

  static constexpr std::size_t dim = N;
  using scalar_type = T;

  template <std::size_t I>
  static constexpr inline scalar_type &write_element(::Utils::Vector<T, N> &v) {
    return v[I];
  }

  template <std::size_t I>
  static constexpr inline scalar_type
  read_element(::Utils::Vector<T, N> const &v) {
    return v[I];
  }

  static inline scalar_type read_element_idx(std::size_t i,
                                             ::Utils::Vector<T, N> const &v) {
    return v[i];
  }
  static inline scalar_type &write_element_idx(std::size_t i,
                                               ::Utils::Vector<T, N> &v) {
    return v[i];
  }
};

template <typename T> struct deduce_vec<Utils::Vector<T, 3>, 3> {
  using type = Utils::Vector<T, 3>;
};

} // namespace boost::qvm

UTILS_ARRAY_BOOST_MPI_T(Utils::Vector, N)
UTILS_ARRAY_BOOST_BIT_S(Utils::Vector, N)
UTILS_ARRAY_BOOST_CLASS(Utils::Vector, N, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::Vector, N, track_never)
