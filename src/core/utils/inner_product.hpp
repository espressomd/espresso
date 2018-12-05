#ifndef UTILS_INNER_PRODUCT_HPP
#define UTILS_INNER_PRODUCT_HPP

#include <array>

namespace Utils {

template <int I, std::size_t N, typename T, int end> struct inner_product_impl {
  double operator()(std::array<int, N> const &left_array,
                    std::array<T, N> const &right_array) const {
    if (left_array[I] == 0) {
      return inner_product_impl<I + 1, N, T, end>{}(left_array, right_array);
    } else {
      return left_array[I] * right_array[I] +
             inner_product_impl<I + 1, N, T, end>{}(left_array, right_array);
    }
  }
};

template <int I, std::size_t N, typename T>
struct inner_product_impl<I, N, T, I> {
  double operator()(std::array<int, N> const &,
                    std::array<T, N> const &) const {
    return 0.0;
  }
};

template <typename T, std::size_t N>
double inner_product(const std::array<int, N> &left_array,
                     const std::array<T, N> &right_array) {
  if (N > 18) {
    auto const s1 =
        inner_product_impl<0, N, T, N / 4>{}(left_array, right_array);
    auto const s2 =
        inner_product_impl<N / 4, N, T, N / 2>{}(left_array, right_array);
    auto const s3 =
        inner_product_impl<N / 2, N, T, 3 * N / 4>{}(left_array, right_array);
    auto const s4 =
        inner_product_impl<3 * N / 4, N, T, N>{}(left_array, right_array);
    return s1 + s2 + s3 + s4;
  }
  return inner_product_impl<0, N, T, N>{}(left_array, right_array);
}

} // namespace Utils

#endif
