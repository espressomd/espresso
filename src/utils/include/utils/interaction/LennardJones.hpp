#ifndef UTILS_LENNARD_JONES
#define UTILS_LENNARD_JONES

#include "utils/Vector.hpp"
#include "utils/math/int_pow.hpp"
#include "utils/math/sqr.hpp"

namespace Utils {
namespace Interaction {

struct LennardJones {
  double epsilon = 0.0;
  double sigma = 0.0;

  double U(double r) const {
    auto const temp = Utils::int_pow<6>(sigma / r);
    return 4.0 * epsilon * (Utils::sqr(temp) - temp);
  }

  Utils::Vector3d F(double r, Utils::Vector3d const &r_ij) const {
    auto const sigma6 = Utils::int_pow<6>(sigma);
    return (48.0 * epsilon * Utils::int_pow<2>(sigma6) / Utils::int_pow<13>(r) -
            24.0 * epsilon * sigma6 / Utils::int_pow<7>(r)) *
           r_ij / r_ij.norm();
  }
};

} // namespace Interaction
} // namespace Utils
#endif
