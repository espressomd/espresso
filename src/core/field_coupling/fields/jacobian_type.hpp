/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_FIELD_COUPLING_GRADIENT_TYPE_HPP
#define CORE_FIELD_COUPLING_GRADIENT_TYPE_HPP

#include "utils/Vector.hpp"

namespace FieldCoupling {
namespace Fields {
namespace detail {
template <class T, size_t codim> struct jacobian_type_impl {
  using type = Utils::Vector<Utils::Vector<T, 3>, codim>;
};

template <class T> struct jacobian_type_impl<T, 1> {
  using type = Utils::Vector<T, 3>;
};

/**
 * @brief Deduce type for jacobian from codim.
 *
 * Small helper that returns Vector3d if codim = 1,
 * and Utils::Vector<codim, Utils::Vector<3, T>> otherwise to avoid
 * using Vectors of size one, where scalars would do.
 */
template <class T, size_t codim>
using jacobian_type = typename jacobian_type_impl<T, codim>::type;
} // namespace detail
} // namespace Fields
} // namespace FieldCoupling

#endif
