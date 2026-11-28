/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include <utils/matrix.hpp>

#include <concepts>
#include <cstddef>

namespace FieldCoupling {
namespace Fields {
namespace detail {
template <class T, std::size_t codim> struct jacobian_type_impl {
  using type = Utils::Matrix<T, codim, 3>;
};

template <class T> struct jacobian_type_impl<T, 1> {
  using type = Utils::Vector<T, 3>;
};

/**
 * @brief Deduce type for Jacobian from codim.
 *
 * Small helper that returns `Utils::Vector3<T>` if codim = 1,
 * and `Utils::Matrix<T, codim, 3>` otherwise to avoid
 * using vectors of size one, where scalars would do.
 */
template <std::floating_point T, std::size_t codim>
using jacobian_type = typename jacobian_type_impl<T, codim>::type;
} // namespace detail
} // namespace Fields
} // namespace FieldCoupling
