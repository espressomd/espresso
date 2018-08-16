#ifndef UTILS_INTERPOLATION_HPP
#define UTILS_INTERPOLATION_HPP

#include "Vector.hpp"

#include "detail/ll_and_dist.hpp"

#include "utils/math/bspline.hpp"

namespace Utils {
namespace Interpolation {

/**
 * @brief cardinal B-spline interpolation with internal iteration.
 */
template <size_t order, typename Kernel>
void bspline_3d(const Vector3d &pos, const Kernel &kernel,
                const Vector3d &grid_spacing, const Vector3d &offset) {
  using Utils::bspline;

  /* The coordinates and relative distance of the assigment cube. */
  const auto block = detail::ll_and_dist<order>(pos, grid_spacing, offset);

  /* Precalc weights that are used multiple times. */
  std::array<double, order> w_y;
  std::array<double, order> w_z;
  for (int i = 0; i < order; i++) {
    w_y[i] = bspline<order>(i, block.distance[1]);
    w_z[i] = bspline<order>(i, block.distance[2]);
  }

  std::array<int, 3> ind;
  for (int i = 0; i < order; i++) {
    ind[0] = block.corner[0] + i;
    const auto wx = bspline<order>(i, block.distance[0]);
    for (int j = 0; j < order; j++) {
      ind[1] = block.corner[1] + j;
      const auto wxy = wx * w_y[j];
      for (int k = 0; k < order; k++) {
        ind[2] = block.corner[2] + k;
        kernel(ind, wxy * w_z[k]);
      }
    }
  }
}

/**
 * @brief cardinal B-spline interpolation with internal iteration.
 */
template <typename Kernel>
void bspline_3d(const Vector3d &pos, const Kernel &kernel,
                const Vector3d &grid_spacing, const Vector3d &offset,
                int order) {

  switch (order) {
  case 1:
    bspline_3d<1>(pos, kernel, grid_spacing, offset);
    break;
  case 2:
    bspline_3d<2>(pos, kernel, grid_spacing, offset);
    break;
  case 3:
    bspline_3d<3>(pos, kernel, grid_spacing, offset);
    break;
  case 4:
    bspline_3d<4>(pos, kernel, grid_spacing, offset);
    break;
  case 5:
    bspline_3d<5>(pos, kernel, grid_spacing, offset);
    break;
  case 6:
    bspline_3d<6>(pos, kernel, grid_spacing, offset);
    break;
  case 7:
    bspline_3d<7>(pos, kernel, grid_spacing, offset);
    break;

  default:
    throw std::runtime_error("Interpolation order not implemented.");
  }
}

/**
 * @brief cardinal B-spline weighted sum.
 */
template <typename T, typename Kernel>
T bspline_3d_accumulate(const Vector3d &pos, const Kernel &kernel,
                        const Vector3d &grid_spacing, const Vector3d &offset,
                        int order, T const &init) {
  T value = init;
  bspline_3d(pos, [&value, &kernel](const std::array<int, 3> &ind,
                                    double w) { value += w * kernel(ind); },
             grid_spacing, offset, order);

  return value;
}
}
}

#endif
