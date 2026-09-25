/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>

/**
 * @brief Write-through reference to one particle's 3-vector stored in a
 * strided column (any Kokkos layout).
 *
 * Component @c j of the referenced vector lives at <tt>base[j * stride]</tt>.
 * The stride is supplied by the caller (@c view.stride(1) for a Kokkos column:
 * 1 for particle-major LayoutRight, number-of-rows for component-major
 * LayoutLeft), so the proxy is layout-agnostic.
 * The proxy is a cheap value type (pointer + stride); copying it copies the
 * reference, not the data. Reads convert to @ref Utils::Vector<T,3> by value.
 *
 * Note: because @ref Utils::Vector arithmetic operators are function
 * templates, a BasicVectorReference does not participate in template argument
 * deduction — expressions like <tt>proxy - vector</tt> require an explicit
 * conversion: <tt>Utils::Vector<T,3>(proxy) - vector</tt>.
 */
template <class T> class BasicVectorReference {
  T *m_base;
  std::size_t m_stride;

public:
  BasicVectorReference(T *base, std::size_t stride)
      : m_base(base), m_stride(stride) {
    assert(base != nullptr);
  }

  operator Utils::Vector<T, 3>() const {
    return {m_base[0u], m_base[m_stride], m_base[2u * m_stride]};
  }

  BasicVectorReference &operator=(Utils::Vector<T, 3> const &value) {
    m_base[0u] = value[0u];
    m_base[m_stride] = value[1u];
    m_base[2u * m_stride] = value[2u];
    return *this;
  }

  /** Copying the proxy rebinds the reference; assignment writes through. */
  BasicVectorReference(BasicVectorReference const &) = default;

  // Copy-assignment writes VALUES through (like std::vector<bool>::reference),
  // not the pointer/stride. The copy CONSTRUCTOR remains implicit and rebinds.
  // Self-assignment is benign for a write-through proxy (elements are copied
  // onto themselves); the guard makes that explicit (clang-tidy
  // bugprone-unhandled-self-assignment).
  BasicVectorReference &operator=(BasicVectorReference const &other) {
    if (this == &other) {
      return *this;
    }
    return *this = Utils::Vector<T, 3>(other);
  }

  BasicVectorReference &operator+=(Utils::Vector<T, 3> const &value) {
    m_base[0u] += value[0u];
    m_base[m_stride] += value[1u];
    m_base[2u * m_stride] += value[2u];
    return *this;
  }

  BasicVectorReference &operator-=(Utils::Vector<T, 3> const &value) {
    m_base[0u] -= value[0u];
    m_base[m_stride] -= value[1u];
    m_base[2u * m_stride] -= value[2u];
    return *this;
  }

  BasicVectorReference &operator*=(T const factor) {
    m_base[0u] *= factor;
    m_base[m_stride] *= factor;
    m_base[2u * m_stride] *= factor;
    return *this;
  }

  T &operator[](std::size_t const j) {
    assert(j < 3u);
    return m_base[j * m_stride];
  }
  T const &operator[](std::size_t const j) const {
    assert(j < 3u);
    return m_base[j * m_stride];
  }

  T norm2() const {
    auto const x = m_base[0u];
    auto const y = m_base[m_stride];
    auto const z = m_base[2u * m_stride];
    return x * x + y * y + z * z;
  }
  T norm() const { return std::sqrt(norm2()); }
};

using VectorReference = BasicVectorReference<double>;
using IntegerVectorReference = BasicVectorReference<int>;

/**
 * @brief Write-through reference to one particle's quaternion stored in a
 * strided column (any Kokkos layout).
 *
 * Component @c j of the referenced quaternion lives at
 * <tt>base[j * stride]</tt> (stride 1 for particle-major LayoutRight). The
 * proxy is a cheap value type (pointer + stride); copying it copies the
 * reference, not the data. Reads convert to @ref Utils::Quaternion<double> by
 * value.
 */
class QuaternionReference {
  double *m_base;
  std::size_t m_stride;

public:
  QuaternionReference(double *base, std::size_t stride)
      : m_base(base), m_stride(stride) {
    assert(base != nullptr);
  }

  operator Utils::Quaternion<double>() const {
    Utils::Quaternion<double> q;
    q[0] = m_base[0u];
    q[1] = m_base[m_stride];
    q[2] = m_base[2u * m_stride];
    q[3] = m_base[3u * m_stride];
    return q;
  }

  QuaternionReference &operator=(Utils::Quaternion<double> const &value) {
    m_base[0u] = value[0u];
    m_base[m_stride] = value[1u];
    m_base[2u * m_stride] = value[2u];
    m_base[3u * m_stride] = value[3u];
    return *this;
  }

  /** Copying the proxy rebinds the reference; assignment writes through. */
  QuaternionReference(QuaternionReference const &) = default;

  // Copy-assignment writes VALUES through (like std::vector<bool>::reference),
  // not the pointer/stride. The copy CONSTRUCTOR remains implicit and rebinds.
  // Self-assignment guard: see BasicVectorReference::operator=.
  QuaternionReference &operator=(QuaternionReference const &other) {
    if (this == &other) {
      return *this;
    }
    return *this = Utils::Quaternion<double>(other);
  }

  QuaternionReference &operator+=(Utils::Quaternion<double> const &value) {
    m_base[0u] += value[0u];
    m_base[m_stride] += value[1u];
    m_base[2u * m_stride] += value[2u];
    m_base[3u * m_stride] += value[3u];
    return *this;
  }

  QuaternionReference &operator/=(double const divisor) {
    m_base[0u] /= divisor;
    m_base[m_stride] /= divisor;
    m_base[2u * m_stride] /= divisor;
    m_base[3u * m_stride] /= divisor;
    return *this;
  }

  double &operator[](std::size_t const j) {
    assert(j < 4u);
    return m_base[j * m_stride];
  }
  double const &operator[](std::size_t const j) const {
    assert(j < 4u);
    return m_base[j * m_stride];
  }

  double norm2() const {
    auto const w = m_base[0u];
    auto const x = m_base[m_stride];
    auto const y = m_base[2u * m_stride];
    auto const z = m_base[3u * m_stride];
    return w * w + x * x + y * y + z * z;
  }
  double norm() const { return std::sqrt(norm2()); }
};
