#ifndef UTILS_INNER_PRODUCT_HPP
#define UTILS_INNER_PRODUCT_HPP

#include <array>

namespace Utils {

template <int I, std::size_t N, typename T, bool end>
struct inner_product_impl {
  double operator()(std::array<int, N> const &left_array,
                    std::array<T, N> const &right_array) const {
    if (left_array[I] == 0) {
      return inner_product_impl<I + 1, N, T, I + 1 == N - 1>{}(left_array,
                                                               right_array);
    } else {
      return left_array[I] * right_array[I] +
             inner_product_impl<I + 1, N, T, I + 1 == N - 1>{}(left_array,
                                                               right_array);
    }
  }
};

template <int I, std::size_t N, typename T>
struct inner_product_impl<I, N, T, true> {
  double operator()(std::array<int, N> const &left_array,
                    std::array<T, N> const &right_array) const {
    if (left_array[N - 1] == 0) {
      return 0.0;
    } else {
      return left_array[N - 1] * right_array[N - 1];
    }
  }
};

template <typename T, std::size_t N>
double inner_product(const std::array<int, N> &left_array,
                     const std::array<T, N> &right_array) {
  return inner_product_impl<0, N, T, N == 0>{}(left_array, right_array);
}

} // namespace Utils

#endif
