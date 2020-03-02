/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP

#include "core/field_coupling/fields/AffineMap.hpp"
#include "core/field_coupling/fields/Constant.hpp"
#include "core/field_coupling/fields/Interpolated.hpp"
#include "core/field_coupling/fields/PlaneWave.hpp"

#include "script_interface/ScriptInterface.hpp"

#include <boost/multi_array.hpp>

namespace ScriptInterface {
namespace Constraints {

namespace detail {
using namespace ::FieldCoupling::Fields;

/**
 * @brief ScriptInterface implementations for the
 *        various fields provided.
 *
 * These are separated from the Constraints because
 * they can be reused together with the fields themselves.
 */
template <typename Field> struct field_params_impl;

template <typename T, size_t codim>
struct field_params_impl<Constant<T, codim>> {
  static Constant<T, codim> make(const VariantMap &params) {
    return Constant<T, codim>{
        get_value<typename Constant<T, codim>::value_type>(params, "value")};
  }
  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"value", AutoParameter::read_only,
             [this_]() { return this_().value(); }}};
  }
};

template <typename T, size_t codim>
struct field_params_impl<AffineMap<T, codim>> {
  using jacobian_type = typename AffineMap<T, codim>::jacobian_type;
  using value_type = typename AffineMap<T, codim>::value_type;

  static AffineMap<T, codim> make(const VariantMap &params) {
    return AffineMap<T, codim>{
        get_value<jacobian_type>(params, "A"),
        get_value_or<value_type>(params, "b", value_type{})};
  }

  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"A", AutoParameter::read_only, [this_]() { return this_().A(); }},
            {"b", AutoParameter::read_only, [this_]() { return this_().b(); }}};
  }
};

template <typename T, size_t codim>
struct field_params_impl<PlaneWave<T, codim>> {
  using jacobian_type = typename PlaneWave<T, codim>::jacobian_type;
  using value_type = typename PlaneWave<T, codim>::value_type;

  static PlaneWave<T, codim> make(const VariantMap &params) {
    return PlaneWave<T, codim>{get_value<value_type>(params, "amplitude"),
                               get_value<value_type>(params, "wave_vector"),
                               get_value<T>(params, "frequency"),
                               get_value_or<T>(params, "phase", 0.)};
  }

  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"amplitude", AutoParameter::read_only,
             [this_]() { return this_().amplitude(); }},
            {"wave_vector", AutoParameter::read_only,
             [this_]() { return this_().k(); }},
            {"frequency", AutoParameter::read_only,
             [this_]() { return this_().omega(); }},
            {"phase", AutoParameter::read_only,
             [this_]() { return this_().phase(); }}};
  }
};

template <typename T, size_t codim>
struct field_params_impl<Interpolated<T, codim>> {
  static Interpolated<T, codim> make(const VariantMap &params) {
    auto const field_data =
        get_value<std::vector<double>>(params, "_field_data");
    auto const field_shape = get_value<Utils::Vector3i>(params, "_field_shape");
    auto const field_codim = get_value<int>(params, "_field_codim");

    if (field_codim != codim) {
      throw std::runtime_error(
          "Field data has the wrong dimensions, needs to be [n, m, o, " +
          std::to_string(codim) + ']');
    }

    if (*std::min_element(field_shape.begin(), field_shape.end()) < 1) {
      throw std::runtime_error("Field is too small, needs to be at least "
                               "one in all directions.");
    }

    auto const grid_spacing =
        get_value<Utils::Vector3d>(params, "grid_spacing");
    auto const origin = -0.5 * grid_spacing;

    using field_data_type =
        typename Utils::decay_to_scalar<Utils::Vector<T, codim>>::type;
    auto array_ref = boost::const_multi_array_ref<field_data_type, 3>(
        reinterpret_cast<const field_data_type *>(field_data.data()),
        field_shape);

    return Interpolated<T, codim>{array_ref, grid_spacing, origin};
  }

  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"grid_spacing", AutoParameter::read_only,
             [this_]() { return this_().grid_spacing(); }},
            {"origin", AutoParameter::read_only,
             [this_]() { return this_().origin(); }},
            {"_field_shape", AutoParameter::read_only,
             [this_]() { return this_().shape(); }},
            {"_field_codim", AutoParameter::read_only,
             []() { return static_cast<int>(codim); }},
            {"_field_data", AutoParameter::read_only,
             [this_]() { return this_().field_data_flat(); }}};
  }
};

template <typename Field, typename T>
static std::vector<AutoParameter> field_parameters(const T &this_) {
  return field_params_impl<Field>::params(this_);
}

template <typename Field> Field make_field(const VariantMap &params) {
  return field_params_impl<Field>::make(params);
}
} // namespace detail
} // namespace Constraints
} // namespace ScriptInterface

#endif
