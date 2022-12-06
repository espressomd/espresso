/*
 * Copyright (C) 2014-2022 The ESPResSo project
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
#ifndef SRC_UTILS_INCLUDE_UTILS_VECTOR_HPP
#define SRC_UTILS_INCLUDE_UTILS_VECTOR_HPP

/**
 * @file
 *
 * @brief Vector implementation and trait types
 * for boost qvm interoperability.
 */

#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/vec_traits.hpp>

#include "utils/Array.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

namespace Utils {

template <typename T, std::size_t N> class Vector : public Array<T, N> {
  using Base = Array<T, N>;

public:
  using Base::Base;
  using Array<T, N>::at;
  using Array<T, N>::operator[];
  using Array<T, N>::front;
  using Array<T, N>::back;
  using Array<T, N>::data;
  using Array<T, N>::begin;
  using Array<T, N>::cbegin;
  using Array<T, N>::end;
  using Array<T, N>::cend;
  using Array<T, N>::empty;
  using Array<T, N>::size;
  using Array<T, N>::max_size;
  using Array<T, N>::fill;
  using Array<T, N>::broadcast;
  Vector() = default;
  Vector(Vector const &) = default;
  Vector &operator=(Vector const &) = default;

  void swap(Vector &rhs) { std::swap_ranges(begin(), end(), rhs.begin()); }

private:
  constexpr void copy_init(T const *first, T const *last) {
    auto it = begin();
    while (first != last) {
      *it++ = *first++;
    }
  }

public:
  template <class Range>
  explicit Vector(Range const &rng) : Vector(std::begin(rng), std::end(rng)) {}
  explicit constexpr Vector(T const (&v)[N]) : Base() {
    copy_init(std::begin(v), std::end(v));
  }

  constexpr Vector(std::initializer_list<T> v) : Base() {
    if (N != v.size()) {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }

    copy_init(v.begin(), v.end());
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

  /**
   * @brief Create a vector that has all entries set to
   *         one value.
   */
  static Vector<T, N> broadcast(T const &s) {
    Vector<T, N> ret;
    std::fill(ret.begin(), ret.end(), s);

    return ret;
  }

  std::vector<T> as_vector() const { return std::vector<T>(begin(), end()); }

  operator std::vector<T>() const { return as_vector(); }

  template <class U> explicit operator Vector<U, N>() const {
    Vector<U, N> ret;

    std::transform(begin(), end(), ret.begin(),
                   [](auto e) { return static_cast<U>(e); });

    return ret;
  }

  T norm2() const { return (*this) * (*this); }
  T norm() const { return std::sqrt(norm2()); }

  /*
   * @brief Normalize the vector.
   *
   * Normalize the vector by its length,
   * if not zero, otherwise the vector is unchanged.
   */

  Vector &normalize() {
    auto const l = norm();
    if (l > T(0)) {
      for (int i = 0; i < N; i++)
        this->operator[](i) /= l;
    }

    return *this;
  }

  Vector normalized() const { return (*this) / (*this).norm(); }
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
  using std::declval;

  using R = decltype(op(declval<T>(), declval<U>()));
  Vector<R, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 op);

  return ret;
}

template <std::size_t N, typename T, typename Op>
Vector<T, N> &binary_op_assign(Vector<T, N> &a, Vector<T, N> const &b, Op op) {
  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(a), op);
  return a;
}

template <std::size_t N, typename T, typename Op>
constexpr bool all_of(Vector<T, N> const &a, Vector<T, N> const &b, Op op) {
  for (int i = 0; i < a.size(); i++) {
    /* Short circuit */
    if (!static_cast<bool>(op(a[i], b[i]))) {
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
Vector<T, N> &operator+=(Vector<T, N> &a, Vector<T, N> const &b) {
  return detail::binary_op_assign(a, b, std::plus<T>());
}

template <std::size_t N, typename T, typename U>
auto operator-(Vector<T, N> const &a, Vector<U, N> const &b) {
  return detail::binary_op(a, b, std::minus<>());
}

template <std::size_t N, typename T>
Vector<T, N> operator-(Vector<T, N> const &a) {
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(ret),
                 [](T const &v) { return -v; });

  return ret;
}

template <std::size_t N, typename T>
Vector<T, N> &operator-=(Vector<T, N> &a, Vector<T, N> const &b) {
  return detail::binary_op_assign(a, b, std::minus<T>());
}

/* Scalar multiplication */
template <std::size_t N, typename T, class U,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(U const &a, Vector<T, N> const &b) {
  using R = decltype(a * std::declval<T>());
  Vector<R, N> ret;

  std::transform(std::begin(b), std::end(b), std::begin(ret),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <std::size_t N, typename T, class U,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(Vector<T, N> const &b, U const &a) {
  using R = decltype(std::declval<T>() * a);
  Vector<R, N> ret;

  std::transform(std::begin(b), std::end(b), std::begin(ret),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <std::size_t N, typename T>
Vector<T, N> &operator*=(Vector<T, N> &b, T const &a) {
  std::transform(std::begin(b), std::end(b), std::begin(b),
                 [a](T const &val) { return a * val; });
  return b;
}

/* Scalar division */
template <std::size_t N, typename T>
Vector<T, N> operator/(Vector<T, N> const &a, T const &b) {
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), ret.begin(),
                 [b](T const &val) { return val / b; });
  return ret;
}

template <std::size_t N, typename T>
Vector<T, N> &operator/=(Vector<T, N> &a, T const &b) {
  std::transform(std::begin(a), std::end(a), std::begin(a),
                 [b](T const &val) { return val / b; });
  return a;
}

namespace detail {
template <class T> struct is_vector : std::false_type {};
template <class T, std::size_t N>
struct is_vector<Vector<T, N>> : std::true_type {};
} // namespace detail

/* Scalar product */
template <std::size_t N, typename T, class U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto operator*(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  return std::inner_product(std::begin(a), std::end(a), std::begin(b), R{});
}

template <
    std::size_t N, typename T, class U,
    class = std::enable_if_t<std::is_integral_v<T> and std::is_integral_v<U>>>
auto operator%(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() % declval<U>());
  Vector<R, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 [](T const &ai, U const &bi) { return ai % bi; });

  return ret;
}

/* Componentwise square root */
template <std::size_t N, typename T> Vector<T, N> sqrt(Vector<T, N> const &a) {
  using std::sqrt;
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), ret.begin(),
                 [](T const &v) { return sqrt(v); });

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
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(a.cbegin(), a.cend(), b.cbegin(), ret.begin(),
                 [](auto ai, auto bi) { return ai * bi; });

  return ret;
}

// specializations for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <class T, class U, std::size_t N,
          class = std::enable_if_t<not(detail::is_vector<T>::value)>>
auto hadamard_product(T const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a * b;

  return ret;
}

template <class T, class U, std::size_t N,
          class = std::enable_if_t<not(detail::is_vector<U>::value)>>
auto hadamard_product(Vector<T, N> const &a, U const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a * b;

  return ret;
}

template <typename T, typename U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto hadamard_product(T const &a, U const &b) {
  return a * b;
}

template <class T, class U, std::size_t N>
auto hadamard_division(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(a.cbegin(), a.cend(), b.cbegin(), ret.begin(),
                 [](auto ai, auto bi) { return ai / bi; });

  return ret;
}

// specializations for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <class T, class U, std::size_t N,
          class = std::enable_if_t<not(detail::is_vector<U>::value)>>
auto hadamard_division(Vector<T, N> const &a, U const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a / b;

  return ret;
}

template <class T, class U, std::size_t N,
          class = std::enable_if_t<not(detail::is_vector<T>::value)>>
auto hadamard_division(T const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(std::begin(b), std::end(b), ret.begin(),
                 [a](T const &bi) { return a / bi; });
  return ret;
}

template <typename T, typename U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto hadamard_division(T const &a, U const &b) {
  return a / b;
}

template <typename T> Vector<T, 3> unit_vector(int i) {
  if (i == 0)
    return {T{1}, T{0}, T{0}};
  if (i == 1)
    return {T{0}, T{1}, T{0}};
  if (i == 2)
    return {T{0}, T{0}, T{1}};
  throw std::domain_error("coordinate out of range");
}

/**
 * @brief Meta function to turns a Vector<1, T> into T.
 */
template <typename T> struct decay_to_scalar {};
template <typename T, std::size_t N> struct decay_to_scalar<Vector<T, N>> {
  using type = Vector<T, N>;
};

template <typename T> struct decay_to_scalar<Vector<T, 1>> { using type = T; };

template <std::size_t I, class T, std::size_t N>
struct tuple_element<I, Vector<T, N>> {
  using type = T;
};

template <class T, std::size_t N>
struct tuple_size<Vector<T, N>> : std::integral_constant<std::size_t, N> {};

template <std::size_t I, class T, std::size_t N>
auto get(Vector<T, N> const &a) -> std::enable_if_t<(I < N), const T &> {
  return a[I];
}
} // namespace Utils
namespace boost {
namespace qvm {

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
  using type = typename Utils::Vector<T, 3>;
};

} // namespace qvm
} // namespace boost

UTILS_ARRAY_BOOST_MPI_T(Utils::Vector, N)
UTILS_ARRAY_BOOST_BIT_S(Utils::Vector, N)
UTILS_ARRAY_BOOST_CLASS(Utils::Vector, N, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::Vector, N, track_never)

#endif
