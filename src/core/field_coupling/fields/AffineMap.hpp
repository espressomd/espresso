#ifndef CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_CONSTANT_HPP

#include "Vector.hpp"

namespace FieldCoupling {
namespace Fields {
/**
 * @brief Affine transform of an vector field.
 *
 * Returns A * x + b, where * is matrix multiplication.
 */
template <typename T, size_t codim> class AffineMap {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using gradient_type = Vector<3, Vector3d>;

private:
  gradient_type m_A;
  value_type m_b;

public:
  AffineMap(const gradient_type &A, const value_type &b) : m_A(A), m_b(b) {}

  gradient_type &A() { return m_A; }
  value_type &b() { return m_b; }

  template <typename F>
  value_type operator()(const F &f, const Vector3d &pos) const {

    value_type ret;

    std::transform(m_A.begin(), m_A.end(), m_b.begin(), ret.begin(),
                   [&pos](const value_type &A_i, const T &b_i) -> T {
                     return A_i * pos + b_i;
                   });

    return f(ret);
  }

  template <typename F>
  gradient_type gradient(const F &f, const Vector3d &) const {
    return f(m_A);
  }

  bool fits_in_box(const Vector3d &) const { return true; }
};
}
}

#endif
