#ifndef UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP
#define UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP

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
template <size_t order> constexpr double pos_shift() {
  return (order % 2 == 0) ? 0.0 : 0.5;
}
}
}
}

#endif
