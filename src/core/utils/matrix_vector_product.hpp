#ifndef UTILS_MATRIX_VECTOR_PRODUCT_HPP
#define UTILS_MATRIX_VECTOR_PRODUCT_HPP

#include <array>
#include <utility>

namespace Utils {

namespace detail {

template <int c, typename T> struct mul {
  constexpr T operator()(const T a) const { return c * a; }
};

template <typename T> struct mul<0, T> {
  constexpr T operator()(const T a) const { return T{}; }
};

template <typename T> struct mul<1, T> {
  constexpr T operator()(const T a) const { return a; }
};

template <int I, typename T, std::size_t N, int c, int... cs>
struct inner_product_template_impl {
  constexpr T operator()(std::array<T, N> const &a) const {
    return mul<c, T>{}(std::get<I>(a)) +
           inner_product_template_impl<I + 1, T, N, cs...>{}(a);
  }
};

template <int I, typename T, std::size_t N, int c>
struct inner_product_template_impl<I, T, N, c> {
  constexpr T operator()(std::array<T, N> const &a) const {
    return mul<c, T>{}(std::get<I>(a));
  }
};

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix,
          std::size_t row_index, std::size_t... column_indices>
constexpr T inner_product_helper(const std::array<T, N> &vec,
                                 std::index_sequence<column_indices...>) {
  return inner_product_template_impl<
      0, T, N, std::get<column_indices>(std::get<row_index>(matrix))...>{}(vec);
}

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix,
          std::size_t row_index>
constexpr T inner_product_template(const std::array<T, N> &vec) {
  return detail::inner_product_helper<T, N, matrix, row_index>(
      vec, std::make_index_sequence<N>{});
}

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix,
          std::size_t... column_indices>
constexpr std::array<T, N>
matrix_vector_product_helper(const std::array<T, N> &vec,
                             std::index_sequence<column_indices...>) {
  return std::array<T, N>{
      {inner_product_template<T, N, matrix, column_indices>(vec)...}};
}

} // namespace detail

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix>
constexpr std::array<T, N> matrix_vector_product(const std::array<T, N> &vec) {
  return detail::matrix_vector_product_helper<T, N, matrix>(
      vec, std::make_index_sequence<N>{});
}

} // namespace Utils

#endif
