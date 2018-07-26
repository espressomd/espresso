#ifndef UTILS_INTERPOLATION_GRADIENT_HPP
#define UTILS_INTERPOLATION_GRADIENT_HPP

#include "Vector.hpp"

#include "detail/ll_and_dist.hpp"
#include "detail/pos_shift.hpp"
#include "utils/math/bspline.hpp"
#include "utils/math/tensor_product.hpp"

#include <array>

namespace Utils {
namespace Interpolation {

/**
 * @brief cardinal B-spline interpolation with internal iteration.
 */
template <size_t order, typename Kernel>
void bspline_3d_gradient(const Vector3d &pos, const Kernel &kernel,
                         const Vector3d &grid_spacing, const Vector3d &offset) {
  /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
  std::array<double, 3> dist;
  /* Index of the lower left corner of the assignment cube */
  std::array<int, 3> ll;

  std::tie(ll, dist) = detail::ll_and_dist<order>(pos, grid_spacing, offset);

  /* Precalc weights that are used multiple times. */
  std::array<double, order> w_y;
  std::array<double, order> w_z;
  std::array<double, order> dw_y;
  std::array<double, order> dw_z;
  for (int i = 0; i < order; i++) {
    w_y[i] = Utils::bspline<order>(i, dist[1]);
    w_z[i] = Utils::bspline<order>(i, dist[2]);
    dw_y[i] = Utils::bspline_d<order>(i, dist[1]) / grid_spacing[1];
    dw_z[i] = Utils::bspline_d<order>(i, dist[2]) / grid_spacing[2];
  }

  std::array<int, 3> ind;
  for (int i = 0; i < order; i++) {
    ind[0] = ll[0] + i;
    const auto w_x = Utils::bspline<order>(i, dist[0]);
    const auto dw_x = Utils::bspline_d<order>(i, dist[0]) / grid_spacing[0];
    for (int j = 0; j < order; j++) {
      ind[1] = ll[1] + j;
      for (int k = 0; k < order; k++) {
        ind[2] = ll[2] + k;
        kernel(ind, Vector3d{dw_x * w_y[j] * w_z[k], w_x * dw_y[j] * w_z[k],
                             w_x * w_y[j] * dw_z[k]});
      }
    }
  }
}

/**
 * @brief cardinal B-spline interpolation with internal iteration.
 */
template <typename Kernel>
void bspline_3d_gradient(const Vector3d &pos, const Kernel &kernel,
                         const Vector3d &grid_spacing, const Vector3d &offset,
                         int order) {
  switch (order) {
  case 1:
    bspline_3d_gradient<1>(pos, kernel, grid_spacing, offset);
    break;
  case 2:
    bspline_3d_gradient<2>(pos, kernel, grid_spacing, offset);
    break;
  case 3:
    bspline_3d_gradient<3>(pos, kernel, grid_spacing, offset);
    break;
  case 4:
    bspline_3d_gradient<4>(pos, kernel, grid_spacing, offset);
    break;
  case 5:
    bspline_3d_gradient<5>(pos, kernel, grid_spacing, offset);
    break;
  case 6:
    bspline_3d_gradient<6>(pos, kernel, grid_spacing, offset);
    break;
  case 7:
    bspline_3d_gradient<7>(pos, kernel, grid_spacing, offset);
    break;
  default:
    throw std::runtime_error("Interpolation order not implemented.");
  }
}

/**
 * @brief cardinal B-spline weighted sum.
 */
template <typename T, typename Kernel>
T bspline_3d_gradient_accumulate(const Vector3d &pos, const Kernel &kernel,
                                 const Vector3d &grid_spacing,
                                 const Vector3d &offset, int order,
                                 T const &init) {
  T value = init;
  bspline_3d_gradient(
      pos,
      [&value, &kernel](const std::array<int, 3> &ind, const Vector3d &w) {
        value += Utils::tensor_product(kernel(ind), w);
      },
      grid_spacing, offset, order);

  return value;
}
}
}

#endif
