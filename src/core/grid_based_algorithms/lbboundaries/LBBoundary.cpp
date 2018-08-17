#include "LBBoundary.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/lbboundaries.hpp"

namespace LBBoundaries {
  // Vector3d LBBoundary::get_force() {
  //   double tmp_force[3];
  //   int n = 0;
  //   for (auto it = lbboundaries.begin(); it != lbboundaries.end(); ++it, ++n) {
  //     if (&(**it) == this)
  // 	break;
  //   }
  //   if (n == lbboundaries.size())
  //     throw std::runtime_error("You probably tried to get the force of an lbboundary that was not added to system.lbboundaries.");
  //   lbboundary_get_force(n, tmp_force);
  //   return Vector3d{tmp_force[0], tmp_force[1], tmp_force[2]};
  // }
}
