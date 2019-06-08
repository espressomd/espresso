#ifndef UTILS_RASTER_HPP
#define UTILS_RASTER_HPP

#include <boost/multi_array.hpp>

#include "utils/Vector.hpp"

namespace Utils {

/**
 * @brief Raster a function over a regular 3d grid.
 *
 * This evaluates a function over a regular grid and
 * returns the function values at the grid points.
 *
 * @param offset Positon of the lowest grid point
 * @param grid_spacing Grid constant
 * @param size Grid size
 * @param f Function to evaluate
 * @return Function values at the grid points.
 */
template <class T, class F>
auto raster(Vector3d const &offset, Vector3d const &grid_spacing, Vector3i size,
            F f) {

  using R = decltype(f((offset)));
  boost::multi_array<R, 3> res(size);

  Vector3i ind;
  for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[2] = 0; ind[2] < size[2]; ind[2]++) {
        auto const x = offset + Vector3d{ind[0] * grid_spacing[0],
                                         ind[1] * grid_spacing[1],
                                         ind[2] * grid_spacing[2]};

        res(ind) = f(x);
      }

  return res;
}
} // namespace Utils

#endif // UTILS_RASTER_HPP
