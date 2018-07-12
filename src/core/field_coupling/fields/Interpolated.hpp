#ifndef CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_INTERPOLATED_HPP

#include "Vector.hpp"

#include "utils/interpolation/bspline_3d.hpp"
#include "utils/interpolation/bspline_3d_gradient.hpp"
#include "utils/math/tensor_product.hpp"

/* Turn off range checks if release build. */
#if defined(NDEBUG) && !defined(BOOST_DISABLE_ASSERTS)
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>

#include <array>

#include "utils/print.hpp"

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
  using gradient_type = Vector<3, value_type>;
  using storage_type = boost::multi_array<value_type, 3>;

private:
  storage_type m_global_field;
  int m_order = 0;
  Vector3d m_grid_spacing = {0, 0, 0};
  Vector3d m_origin = {0, 0, 0};

public:
  Interpolated(const boost::const_multi_array_ref<value_type, 3> &global_field,
               int interpolation_order, const Vector3d &grid_spacing,
               const Vector3d &origin)
      : m_global_field(global_field), m_order(interpolation_order),
        m_grid_spacing(grid_spacing), m_origin(origin) {}

private:
  void copy(const Interpolated &rhs) {
    detail::deep_copy(m_global_field, rhs.m_global_field);

    m_order = rhs.m_order;
    m_grid_spacing = rhs.m_grid_spacing;
    m_origin = rhs.m_origin;
  }

public:
  Interpolated(const Interpolated &rhs) { copy(rhs); }
  Interpolated &operator=(const Interpolated &rhs) {
    copy(rhs);
    return *this;
  }

  int interpolation_order() const { return m_order; }
  Vector3d grid_spacing() const { return m_grid_spacing; }
  storage_type &field_data() { return m_global_field; }
  Vector3d origin() { return m_origin; }
  Vector<3, int> shape() {
    return {m_global_field.shape(), m_global_field.shape() + 3};
  }

  /*
   * @brief Evaluate f at pos with the field value as argument.
   */
  template <typename F>
  value_type operator()(const F &f, const Vector3d &pos) const {
    value_type value{};

    /* If F is linear we can first interpolate the field value and
       then evaluate the function once on the result. */
    if (F::is_linear) {
      auto const kernel = [&value, this](const std::array<int, 3> &ind,
                                         double weight) {
        value += weight * m_global_field(ind);
      };
      Utils::Interpolation::bspline_3d(pos, kernel, m_grid_spacing, -m_origin,
                                       m_order);

      return f(value);
    } else {
      auto const kernel = [&value, &f, this](const std::array<int, 3> &ind,
                                             double weight) {
        value += weight * f(m_global_field(ind));
      };

      Utils::Interpolation::bspline_3d(pos, kernel, m_grid_spacing, {},
                                       m_order);

      return value;
    }
  }

  /*
   * @brief Evaluate f at pos with the gradient field value as argument.
   */
  template <typename F>
  gradient_type gradient(const F &f, const Vector3d &pos) const {
    gradient_type gradient{};
    auto kernel = [&gradient, &f, this](const std::array<int, 3> &ind,
                                        const Vector3d &weight) {
      auto const field_val = m_global_field(ind);
      auto const f_val = f(Utils::tensor_product(field_val, weight));
      gradient += f_val;
    };

    Utils::Interpolation::bspline_3d_gradient(pos, kernel, m_grid_spacing, {},
                                              m_order);

    return gradient;
  }

  bool fits_in_box(const Vector3d &box) const {
    auto const halo_size = (0.5 * m_order) * m_grid_spacing;

    /* Not enough halo on the left side */
    if (m_origin > -halo_size)
      return false;

    return true;
  }
};
}
}

#endif
