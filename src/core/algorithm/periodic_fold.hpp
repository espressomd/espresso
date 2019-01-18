#ifndef CORE_ALGORITHM_PERIODIC_FOLD_HPP
#define CORE_ALGORITHM_PERIODIC_FOLD_HPP

#include <limits>
#include <utility>

namespace Algorithm {
template <typename T, typename I>
std::pair<T, I> periodic_fold(T x, I i, T const &l) {
  using limits = std::numeric_limits<I>;

  while ((x < 0) && (i > limits::min())) {
    x += l;
    --i;
  }

  while ((x >= l) && (i < limits::max())) {
    x -= l;
    ++i;
  }

  return {x, i};
}

template <typename T> T periodic_fold(T x, T const &l) {
  while (x < 0) {
    x += l;
  }

  while (x >= l) {
    x -= l;
  }

  return x;
}
} // namespace Algorithm

#endif
