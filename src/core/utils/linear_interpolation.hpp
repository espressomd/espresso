#ifndef UTILS_LINEAR_INTERPOLATION_HPP
#define UTILS_LINEAR_INTERPOLATION_HPP

#include <cassert>

namespace Utils {
template <typename T, typename Container>
T linear_interpolation(Container const &table, T hi, T offset, T x) {
  auto const dind = (x - offset) * hi;
  auto const ind = static_cast<int>(dind);
  assert(ind <= dind);
  auto const dx = dind - ind;

  /* linear interpolation between data points */
  return table[ind] * (T{1} - dx) + table[ind + 1] * dx;
}
}

#endif
