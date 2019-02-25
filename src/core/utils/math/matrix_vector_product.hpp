#ifndef UTILS_MATRIX_VECTOR_PRODUCT_HPP
#define UTILS_MATRIX_VECTOR_PRODUCT_HPP

#include <array>
#include <utility>

#include "utils/get.hpp"

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

template <int I, typename T, typename Container, std::size_t N, int c,
          int... cs>
struct inner_product_template_impl {
  constexpr T operator()(const Container &vec) const {
    return mul<c, T>{}(get<I>(vec)) +
           inner_product_template_impl<I + 1, T, Container, N, cs...>{}(vec);
  }
};

template <int I, typename T, typename Container, std::size_t N, int c>
struct inner_product_template_impl<I, T, Container, N, c> {
  constexpr T operator()(const Container &vec) const {
    return mul<c, T>{}(get<I>(vec));
  }
};

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix, typename Container,
          std::size_t row_index, std::size_t... column_indices>
constexpr T inner_product_helper(const Container &vec,
                                 std::index_sequence<column_indices...>) {
  return inner_product_template_impl<
      0, T, Container, N, get<column_indices>(get<row_index>(matrix))...>{}(
      vec);
}

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix, typename Container,
          std::size_t row_index>
constexpr T inner_product_template(const Container &vec) {
  return detail::inner_product_helper<T, N, matrix, Container, row_index>(
      vec, std::make_index_sequence<N>{});
}

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix, typename Container,
          std::size_t... column_indices>
constexpr std::array<T, N>
matrix_vector_product_helper(const Container &vec,
                             std::index_sequence<column_indices...>) {
  return std::array<T, N>{
      {inner_product_template<T, N, matrix, Container, column_indices>(
          vec)...}};
}

} // namespace detail

/**
 * @brief Calculate the matrix-vector product for a statically given (square)
 * matrix.
 * @tparam T data type for the vector.
 * @tparam N size of the vector.
 * @tparam matrix const reference to a static integer array (size N) of arrays
 * (each of size N).
 * @tparam Container container type for the vector.
 * @param vec Container with data of type T and length N with the vector data.
 * @retval An array with the result of the matrix-vector product.
 */
template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &matrix, typename Container>
constexpr std::array<T, N> matrix_vector_product(const Container &vec) {
  return detail::matrix_vector_product_helper<T, N, matrix, Container>(
      vec, std::make_index_sequence<N>{});
}

} // namespace Utils

#endif
