#ifndef SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_DETAIL_FIELDS_HPP

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
    return {{"value", this_().value()}};
  }
};

template <typename T, size_t codim>
struct field_params_impl<Interpolated<T, codim>> {
  static Interpolated<T, codim> make(const VariantMap &params) {
    auto const field = get_value<std::vector<double>>(params, "field");
    auto const shape = get_value<Vector<3, int>>(params, "shape");
    auto const order = get_value<int>(params, "interpolation_order");
    auto const origin = get_value<Vector3d>(params, "origin");
    auto const grid_spacing = get_value<Vector3d>(params, "grid_spacing");

    using field_data_type = typename decay_to_scalar<Vector<codim, T>>::type;

    auto array_ref = boost::const_multi_array_ref<field_data_type, 3>(
        reinterpret_cast<const field_data_type *>(field.data()), shape);

    array_ref.reindex(origin);

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
            {"shape", AutoParameter::read_only,
             [this_]() { return this_().shape(); }},
            {"field", AutoParameter::read_only, [this_]() {
               auto &field_data = this_().field_data();
               double *data_ptr = reinterpret_cast<double *>(field_data.data());
               return std::vector<double>(data_ptr,
                                          data_ptr + field_data.num_elements());
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
