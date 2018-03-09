#ifndef UTILS_MATH_TENSOR_PRODUCT_HPP
#define UTILS_MATH_TENSOR_PRODUCT_HPP

#include <Vector.hpp>

#include <algorithm>

namespace Utils {
template <typename T, size_t N, size_t M>
Vector<N, Vector<M, T>> tensor_product(const Vector<N, T> &x,
                                       const Vector<M, T> &y) {
  Vector<N, Vector<M, T>> ret;

  std::transform(x.begin(), x.end(), ret.begin(),
                 [&y](const T &x_i) { return x_i * y; });

  return ret;
}
}

#endif
