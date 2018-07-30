#ifndef SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP

#include "core/field_coupling/fields/AffineMap.hpp"
#include "core/field_coupling/fields/Constant.hpp"
#include "core/field_coupling/fields/Interpolated.hpp"

#include "ScriptInterface.hpp"

#include <boost/multi_array.hpp>

namespace ScriptInterface {
namespace Constraints {

namespace detail {
using namespace ::FieldCoupling::Fields;

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
  using gradient_type = typename AffineMap<T, codim>::gradient_type;
  using value_type = typename AffineMap<T, codim>::value_type;

  static AffineMap<T, codim> make(const VariantMap &params) {
    return AffineMap<T, codim>{
        get_value<gradient_type>(params, "A"),
        get_value_or<value_type>(params, "b", value_type{})};
  }

  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"A", AutoParameter::read_only, [this_]() { return this_().A(); }},
            {"b", AutoParameter::read_only, [this_]() { return this_().b(); }}};
  }
};

template <typename T, size_t codim>
struct field_params_impl<Interpolated<T, codim>> {
  static Interpolated<T, codim> make(const VariantMap &params) {
    auto const field_data =
        get_value<std::vector<double>>(params, "_field_data");
    auto const field_shape = get_value<Vector<3, int>>(params, "_field_shape");
    auto const field_codim = get_value<int>(params, "_field_codim");

    if (field_codim != codim) {
      throw std::runtime_error(
          "Field data has the wrong dimensions, needs to be [n, m, o, " +
          std::to_string(codim) + ']');
    }

    auto const order = get_value<int>(params, "interpolation_order");
    if (*std::min_element(field_shape.begin(), field_shape.end()) <=
        (order / 2)) {
      throw std::runtime_error("Field is to small, needs to be at least " +
                               std::to_string(order / 2 + 1) +
                               " in all directions.");
    }

    auto const grid_spacing = get_value<Vector3d>(params, "grid_spacing");
    auto const halo_points = static_cast<double>(order / 2);
    auto const origin = -(halo_points + 0.5) * grid_spacing;

    using field_data_type = typename decay_to_scalar<Vector<codim, T>>::type;
    auto array_ref = boost::const_multi_array_ref<field_data_type, 3>(
        reinterpret_cast<const field_data_type *>(field_data.data()),
        field_shape);

    return Interpolated<T, codim>{array_ref, order, grid_spacing, origin};
  }

  template <typename This>
  static std::vector<AutoParameter> params(const This &this_) {
    return {{"interpolation_order", AutoParameter::read_only,
             [this_]() { return this_().interpolation_order(); }},
            {"grid_spacing", AutoParameter::read_only,
             [this_]() { return this_().grid_spacing(); }},
            {"origin", AutoParameter::read_only,
             [this_]() { return this_().origin(); }},
            {"_field_shape", AutoParameter::read_only,
             [this_]() { return this_().shape(); }},
            {"_field_codim", AutoParameter::read_only,
             []() { return static_cast<int>(codim); }},
            {"_field_data", AutoParameter::read_only, [this_]() {
               auto &field_data = this_().field_data();
               auto data_ptr =
                   reinterpret_cast<const double *>(field_data.data());
               return std::vector<double>(
                   data_ptr, data_ptr + codim * field_data.num_elements());
             }}};
  }
};

template <typename Field, typename T>
static std::vector<AutoParameter> field_parameters(const T &this_) {
  return field_params_impl<Field>::params(this_);
}

template <typename Field> Field make_field(const VariantMap &params) {
  return field_params_impl<Field>::make(params);
}
}
}
}

#endif
