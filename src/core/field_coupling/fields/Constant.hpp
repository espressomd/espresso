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
#ifndef CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP

#include "jacobian_type.hpp"
#include <utils/Vector.hpp>

namespace FieldCoupling {
namespace Fields {
/**
 * @brief A vector field that is constant in space.
 */
template <typename T, size_t codim> class Constant {
public:
  using value_type =
      typename Utils::decay_to_scalar<Utils::Vector<T, codim>>::type;
  using jacobian_type = detail::jacobian_type<T, codim>;

private:
  value_type m_value;

public:
  Constant(const value_type &value) : m_value(value) {}

  value_type &value() { return m_value; }

  value_type operator()(const Utils::Vector3d &, double = {}) const {
    return m_value;
  }
  static constexpr jacobian_type jacobian(const Utils::Vector3d &) {
    return jacobian_type{};
  }

  bool fits_in_box(const Utils::Vector3d &) const { return true; }
};
} // namespace Fields
} // namespace FieldCoupling

#endif
