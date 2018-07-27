#ifndef CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP

#include "Vector.hpp"
#include "gradient_type.hpp"

namespace FieldCoupling {
namespace Fields {
/**
 * @brief A vector field that is constant in space.
 */
template <typename T, size_t codim> class Constant {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using gradient_type = detail::gradient_type<T, codim>;

private:
  value_type m_value;

public:
  Constant(const value_type &value) : m_value(value) {}

  value_type &value() { return m_value; }

  template <typename F>
  value_type operator()(const F &f, const Vector3d &) const {
    return f(m_value);
  }

  template <typename F>
  static constexpr gradient_type gradient(const F &f, const Vector3d &) {
    return f(gradient_type{});
  }

  bool fits_in_box(const Vector3d &) const { return true; }
};
}
}

#endif
