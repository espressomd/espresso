#ifndef UTILS_INTERPOLATION_HPP
#define UTILS_INTERPOLATION_HPP

#include "Vector.hpp"

#include "detail/ll_and_dist.hpp"

#include "utils/math/bspline.hpp"

namespace Utils {
namespace Interpolation {
/**
 * @brief cardinal B-spline interpolation with internal iteration.
 *
 * @tparam order The number of interpolation points in each direction.
 * @tparam Kernel Callable that can be called at every point with the
 *         index of the point and it's weight as arguments.
 *
 * @param pos The point at which the interpolation should be run.
 * @param kernel Function to run for every point.
 * @param grid_spacing The distance between the grid points.
 * @param offset Shift of the grid relative to the origin.
 */
template <size_t order, typename Kernel>
void bspline_3d(const Vector3d &pos, const Kernel &kernel,
                const Vector3d &grid_spacing, const Vector3d &offset) {
  using Utils::bspline;

  /* The coordinates and relative distance of the assignment cube. */
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
 * @brief cardinal B-spline weighted sum.
 */
template <size_t order, typename T, typename Kernel>
T bspline_3d_accumulate(const Vector3d &pos, const Kernel &kernel,
                        const Vector3d &grid_spacing, const Vector3d &offset,
                        T const &init) {
  T value = init;
  bspline_3d<order>(pos,
                    [&value, &kernel](const std::array<int, 3> &ind, double w) {
                      value += w * kernel(ind);
                    },
                    grid_spacing, offset);

  return value;
}
} // namespace Interpolation
} // namespace Utils

#endif
