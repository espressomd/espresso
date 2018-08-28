#ifndef UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP
#define UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP

#include <cmath>

namespace Utils {
namespace Interpolation {
namespace detail {
/**
 * @brief Calculates the shift of the aufpunkt relative to the mesh.
 *
 * For even meshes, distances are calculated relative to the middle
 * of two grid points, hence the distances are shifted by half a grid
 * constant. For odd size, distances are relative to the central point,
 * so there is no shift.
 */
template <unsigned order> constexpr double pos_shift() {
  return 0.5 * (order % 2);
}
} // namespace detail
} // namespace Interpolation
} // namespace Utils

#endif
