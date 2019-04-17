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
#ifndef CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP

#include "utils/interpolation/bspline_3d.hpp"
#include "utils/interpolation/bspline_3d_gradient.hpp"
#include "utils/math/tensor_product.hpp"

#include "jacobian_type.hpp"
#include "utils/Vector.hpp"

/* Turn off range checks if release build. */
#if defined(NDEBUG) && !defined(BOOST_DISABLE_ASSERTS)
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>

#include <array>

namespace FieldCoupling {
namespace Fields {
namespace detail {
template <typename T>
void deep_copy(boost::multi_array<T, 3> &dst,
               const boost::multi_array<T, 3> &src) {
  auto *s = src.shape();
  dst.resize(boost::extents[s[0]][s[1]][s[2]]);
  dst = src;

  auto *b = src.index_bases();
  dst.reindex(std::array<typename boost::multi_array<T, 3>::index, 3>{
      {b[0], b[1], b[2]}});
}
} // namespace detail

/**
 * @brief A vector field interpolated from a regular grid.
 *
 * This is an interpolation wrapper around a boost::multi_array,
 * which can be evaluated on any point in space by spline interpolation.
 */
template <typename T, size_t codim> class Interpolated {
public:
  using value_type =
      typename Utils::decay_to_scalar<Utils::Vector<T, codim>>::type;
  using jacobian_type = detail::jacobian_type<T, codim>;
  using storage_type = boost::multi_array<value_type, 3>;

private:
  storage_type m_global_field;
  Utils::Vector3d m_grid_spacing;
  Utils::Vector3d m_origin;

public:
  Interpolated(const boost::const_multi_array_ref<value_type, 3> &global_field,
               const Utils::Vector3d &grid_spacing,
               const Utils::Vector3d &origin)
      : m_global_field(global_field), m_grid_spacing(grid_spacing),
        m_origin(origin) {}

private:
  void copy(const Interpolated &rhs) {
    detail::deep_copy(m_global_field, rhs.m_global_field);

    m_grid_spacing = rhs.m_grid_spacing;
    m_origin = rhs.m_origin;
  }

public:
  Interpolated(const Interpolated &rhs) { copy(rhs); }
  Interpolated &operator=(const Interpolated &rhs) {
    copy(rhs);
    return *this;
  }

  Utils::Vector3d grid_spacing() const { return m_grid_spacing; }
  storage_type const &field_data() const { return m_global_field; }
  Utils::Vector3d origin() const { return m_origin; }
  Utils::Vector3i shape() const {
    return {m_global_field.shape(), m_global_field.shape() + 3};
  }

  /*
   * @brief Evaluate f at pos with the field value as argument.
   */
  value_type operator()(const Utils::Vector3d &pos, double = {}) const {
    using Utils::Interpolation::bspline_3d_accumulate;
    return bspline_3d_accumulate<2>(
        pos,
        [this](const std::array<int, 3> &ind) { return m_global_field(ind); },
        m_grid_spacing, m_origin, value_type{});
  }

  /*
   * @brief Evaluate f at pos with the jacobian field value as argument.
   */
  jacobian_type jacobian(const Utils::Vector3d &pos, double = {}) const {
    using Utils::Interpolation::bspline_3d_gradient_accumulate;
    return bspline_3d_gradient_accumulate<2>(
        pos,
        [this](const std::array<int, 3> &ind) { return m_global_field(ind); },
        m_grid_spacing, m_origin, jacobian_type{});
  }

  bool fits_in_box(const Utils::Vector3d &box) const {
    const Utils::Vector3d grid_size = {m_grid_spacing[0] * shape()[0],
                                       m_grid_spacing[1] * shape()[1],
                                       m_grid_spacing[2] * shape()[2]};
    return (m_origin < Utils::Vector3d::broadcast(0.)) &&
           ((m_origin + grid_size) >= box);
  }
};
} // namespace Fields
} // namespace FieldCoupling

#endif
