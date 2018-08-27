#ifndef CORE_FIELD_COUPLING_GRADIENT_TYPE_HPP
#define CORE_FIELD_COUPLING_GRADIENT_TYPE_HPP

#include "Vector.hpp"

namespace FieldCoupling {
namespace Fields {
namespace detail {
template <class T, size_t codim> struct gradient_type_impl {
  using type = Vector<codim, Vector<3, T>>;
};

template <class T> struct gradient_type_impl<T, 1> {
  using type = Vector<3, T>;
};

/**
 * @brief Deduce type for gradient from codim.
 *
 * Small helper that returns Vector3d if codim = 1,
 * and Vector<codim, Vector<3, T>> otherwise to avoid
 * using Vectors of size one, where scalars would do.
 */
template <class T, size_t codim>
using gradient_type = typename gradient_type_impl<T, codim>::type;
} // namespace detail
} // namespace Fields
} // namespace FieldCoupling

#endif
