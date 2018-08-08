#ifndef CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP

#include "utils/interpolation/bspline_3d.hpp"
#include "utils/interpolation/bspline_3d_gradient.hpp"
#include "utils/math/tensor_product.hpp"

#include "Vector.hpp"
#include "gradient_type.hpp"

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
      b[0], b[1], b[2]});
}

template <typename F, class = void>
struct is_linear : public std::false_type {};

template <typename F>
struct is_linear<F, typename std::enable_if<F::is_linear>::type>
    : std::true_type {};
}

/**
* @brief A vector field interpolated from a regular grid.
*
* This is an interpolation wrapper around a boost::multi_array,
* which can be evaluated on any point in space by spline interpolation.
*/
template <typename T, size_t codim> class Interpolated {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using gradient_type = detail::gradient_type<T, codim>;
  using storage_type = boost::multi_array<value_type, 3>;

private:
  storage_type m_global_field;
  Vector3d m_grid_spacing;
  Vector3d m_origin;

public:
  Interpolated(const boost::const_multi_array_ref<value_type, 3> &global_field,
               const Vector3d &grid_spacing, const Vector3d &origin)
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

  Vector3d grid_spacing() const { return m_grid_spacing; }
  storage_type const &field_data() const { return m_global_field; }
  Vector3d origin() const { return m_origin; }
  Vector<3, int> shape() const {
    return {m_global_field.shape(), m_global_field.shape() + 3};
  }

  /*
   * @brief Evaluate f at pos with the field value as argument.
   */
  template <typename F>
  value_type operator()(const F &f, const Vector3d &pos) const {
    using Utils::Interpolation::bspline_3d_accumulate;

    /* If F is linear we can first interpolate the field value and
       then evaluate the function once on the result. */
    if (detail::is_linear<F>::value) {
      return f(bspline_3d_accumulate(
          pos,
          [this](const std::array<int, 3> &ind) { return m_global_field(ind); },
          m_grid_spacing, m_origin, 2, value_type{}));
    } else {
      return bspline_3d_accumulate(pos,
                                   [this, &f](const std::array<int, 3> &ind) {
                                     return f(m_global_field(ind));
                                   },
                                   m_grid_spacing, m_origin, 2, value_type{});
    }
  }

  /*
   * @brief Evaluate f at pos with the gradient field value as argument.
   */
  template <typename F>
  gradient_type gradient(const F &f, const Vector3d &pos) const {
    using Utils::Interpolation::bspline_3d_gradient_accumulate;

    /* If F is linear we can first interpolate the field value and
       then evaluate the function once on the result. */
    if (detail::is_linear<F>::value) {
      return f(bspline_3d_gradient_accumulate(
          pos,
          [this](const std::array<int, 3> &ind) { return m_global_field(ind); },
          m_grid_spacing, m_origin, 2, gradient_type{}));
    } else {
      return bspline_3d_gradient_accumulate(
          pos,
          [this, &f](const std::array<int, 3> &ind) {
            return f(m_global_field(ind));
          },
          m_grid_spacing, m_origin, 2, gradient_type{});
    }
  }

  bool fits_in_box(const Vector3d &box) const {
    const Vector3d grid_size = {m_grid_spacing[0] * shape()[0],
                                m_grid_spacing[1] * shape()[1],
                                m_grid_spacing[2] * shape()[2]};
    return (m_origin < Vector3d::broadcast(0.)) &&
           ((m_origin + grid_size) >= box);
  }
};
}
}

#endif
