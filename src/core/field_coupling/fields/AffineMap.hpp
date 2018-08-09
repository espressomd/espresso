#ifndef CORE_EXTERNAL_FIELD_FIELDS_AFFINE_MAP_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_AFFINE_MAP_HPP

#include "Vector.hpp"
#include "gradient_type.hpp"

namespace FieldCoupling {
namespace Fields {

namespace detail {

template <typename T, size_t codim> struct matrix_vector_impl {
  Vector<codim, T> operator()(const Vector<codim, Vector<3, T>> &A,
                              Vector<3, T> const &v) const {
    Vector<codim, T> ret;

    for (int i = 0; i < codim; i++)
      ret[i] = A[i] * v;

    return ret;
  }
};

template <typename T> struct matrix_vector_impl<T, 1> {
  T operator()(const Vector<3, T> &A, Vector<3, T> const &v) const {
    return A * v;
  }
};
}

/**
 * @brief Affine transform of an vector field.
 *
 * Returns A * x + b, where * is matrix multiplication.
 */
template <typename T, size_t codim> class AffineMap {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using gradient_type = detail::gradient_type<T, codim>;

private:
  gradient_type m_A;
  value_type m_b;

public:
  AffineMap(const gradient_type &A, const value_type &b) : m_A(A), m_b(b) {}

  gradient_type &A() { return m_A; }
  value_type &b() { return m_b; }

  template <typename F>
  value_type operator()(const F &f, const Vector3d &pos) const {
    auto const axb = detail::matrix_vector_impl<T, codim>{}(m_A, pos) + m_b;

    return f(axb);
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
