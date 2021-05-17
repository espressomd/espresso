/*
 * Copyright (C) 2010-2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef SRC_UTILS_INCLUDE_UTILS_MATRIX_HPP
#define SRC_UTILS_INCLUDE_UTILS_MATRIX_HPP

#include "utils/Array.hpp"
#include "utils/Vector.hpp"
#include "utils/flatten.hpp"

#include <array>

// These includes need to come first due to ADL reasons.
// clang-format off
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/vec_operations.hpp>
// clang-format on

#include <boost/qvm/deduce_mat.hpp>
#include <boost/qvm/deduce_scalar.hpp>
#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/map_mat_mat.hpp>
#include <boost/qvm/map_mat_vec.hpp>
#include <boost/qvm/map_vec_mat.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_traits.hpp>

/**
 * @file matrix.hpp
 *
 * @brief This file contains a matrix implementation and the trait types needed
 * for the boost qvm interoperability.
 */

namespace Utils {

/**
 * @brief Matrix representation with static size.
 * @tparam T The data type.
 * @tparam Rows Number of rows.
 * @tparam Cols Number of columns.
 */
template <typename T, std::size_t Rows, std::size_t Cols> struct Matrix {
  using container = Utils::Array<T, Cols * Rows>;
  using pointer = typename container::pointer;
  using const_pointer = typename container::const_pointer;
  using iterator = typename container::iterator;
  using const_iterator = typename container::const_iterator;
  using value_type = typename container::value_type;
  using reference = typename container::reference;
  using const_reference = typename container::const_reference;

  container m_data;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &m_data;
  }

  Matrix() = default;
  Matrix(std::initializer_list<T> init_list) {
    assert(init_list.size() == Rows * Cols);
    std::copy(init_list.begin(), init_list.end(), begin());
  }
  Matrix(std::initializer_list<std::initializer_list<T>> init_list) {
    assert(init_list.size() == Rows);
    Utils::flatten(init_list, begin());
  }

  /**
   * @brief Element access (const).
   * @param row The row used for access.
   * @param col The column used for access.
   * @return The matrix element at row @p row and column @p col.
   */
  constexpr value_type operator()(std::size_t row, std::size_t col) const {
    assert(row < Rows);
    assert(col < Cols);
    return m_data[Cols * row + col];
  }
  /**
   * @brief Element access (non const).
   * @param row The row used for access.
   * @param col The column used for access.
   * @return The matrix element at row @p row and column @p col.
   */
  constexpr reference operator()(std::size_t row, std::size_t col) {
    assert(row < Rows);
    assert(col < Cols);
    return m_data[Cols * row + col];
  }

  /**
   * @brief Access to the underlying data pointer (non const).
   * @return Pointer to first element of the data.
   */
  constexpr pointer data() { return m_data.data(); }
  /**
   * @brief Access to the underlying data pointer (non const).
   * @return Pointer to first element of the data.
   */
  constexpr const_pointer data() const noexcept { return m_data.data(); }
  /**
   * @brief Iterator access (non const).
   * @return Returns an iterator to the first element of the matrix.
   */
  constexpr iterator begin() noexcept { return m_data.begin(); };
  /**
   * @brief Iterator access (const).
   * @return Returns an iterator to the first element of the matrix.
   */
  constexpr const_iterator begin() const noexcept { return m_data.begin(); };
  /**
   * @brief Iterator access (non const).
   * @return Returns an iterator to the element following the last element of
   * the matrix.
   */
  constexpr iterator end() noexcept { return m_data.end(); };
  /**
   * @brief Iterator access (non const).
   * @return Returns an iterator to the element following the last element of
   * the matrix.
   */
  constexpr const_iterator end() const noexcept { return m_data.end(); };
  /**
   * @brief Retrieve an entire matrix row.
   * @tparam R The row index.
   * @return A vector containing the elements of row @p R.
   */
  template <std::size_t R> Vector<T, Cols> row() const {
    static_assert(R < Rows, "Invalid row index.");
    return boost::qvm::row<R>(*this);
  }
  /**
   * @brief Retrieve an entire matrix column.
   * @tparam C The column index.
   * @return A vector containing the elements of column @p C.
   */
  template <std::size_t C> Vector<T, Rows> col() const {
    static_assert(C < Cols, "Invalid column index.");
    return boost::qvm::col<C>(*this);
  }
  /**
   * @brief Retrieve the diagonal.
   * @return Vector containing the diagonal elements of the matrix.
   */
  Vector<T, Cols> diagonal() const {
    static_assert(Rows == Cols,
                  "Diagonal can only be retrieved from square matrices.");
    return boost::qvm::diag(*this);
  }
  /**
   * @brief Retrieve the trace.
   * @return Vector containing the sum of diagonal matrix elements.
   */
  T trace() const {
    auto const d = diagonal();
    return std::accumulate(d.begin(), d.end(), T{}, std::plus<T>{});
  }

  /**
   * @brief Retrieve a transposed copy of the matrix.
   * @return Transposed matrix.
   */
  Matrix<T, Cols, Rows> transposed() const {
    return boost::qvm::transposed(*this);
  }

  /**
   * @brief Retrieve an inverted copy of the matrix.
   * @return Inverted matrix.
   */
  Matrix<T, Rows, Cols> inversed() const {
    static_assert(Rows == Cols,
                  "Inversion of a non-square matrix not implemented.");
    return boost::qvm::inverse(*this);
  }
  /**
   * @brief Retrieve the shape of the matrix.
   * @return Pair containing number of rows and number of columns of the matrix.
   */
  constexpr std::pair<std::size_t, std::size_t> shape() const noexcept {
    return {Rows, Cols};
  }
};

using boost::qvm::operator+;
using boost::qvm::operator+=;
using boost::qvm::operator-;
using boost::qvm::operator-=;
using boost::qvm::operator*;
using boost::qvm::operator*=;
using boost::qvm::operator==;

template <typename T, std::size_t M, std::size_t N>
Utils::Vector<T, M * N> flatten(Matrix<T, M, N> const &m) {
  return Utils::Vector<T, M * N>(m.begin(), m.end());
}

template <typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> diagonal_mat(Utils::Vector<T, Rows> const &v) {
  static_assert(Rows == Cols, "Diagonal matrix has to be a square matrix.");
  return boost::qvm::diag_mat(v);
}

template <typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> identity_mat() {
  static_assert(Rows == Cols,
                "Identity matrix only defined for square matrices.");
  return boost::qvm::identity_mat<T, Rows>();
}

} // namespace Utils

namespace boost {
namespace qvm {

template <typename T, std::size_t Rows, std::size_t Cols>
struct mat_traits<Utils::Matrix<T, Rows, Cols>> {
  using mat_type = typename Utils::Matrix<T, Rows, Cols>;
  static int const rows = Rows;
  static int const cols = Cols;
  using scalar_type = T;

  template <std::size_t R, std::size_t C>
  static inline scalar_type read_element(mat_type const &m) {
    static_assert(R < Rows, "Invalid row index.");
    static_assert(C < Cols, "Invalid column index.");
    return m(R, C);
  }

  template <std::size_t R, std::size_t C>
  static inline scalar_type &write_element(mat_type &m) {
    static_assert(R < Rows, "Invalid row index.");
    static_assert(C < Cols, "Invalid column index.");
    return m(R, C);
  }

  static inline scalar_type read_element_idx(std::size_t r, std::size_t c,
                                             mat_type const &m) {
    assert(r < Rows);
    assert(c < Cols);
    return m(r, c);
  }
  static inline scalar_type &write_element_idx(std::size_t r, std::size_t c,
                                               mat_type &m) {
    assert(r < Rows);
    assert(c < Cols);
    return m(r, c);
  }
};

template <typename T, typename U>
struct deduce_vec2<Utils::Matrix<T, 2, 2>, Utils::Vector<U, 2>, 2> {
  using type = Utils::Vector<std::common_type_t<T, U>, 2>;
};

template <typename T, typename U>
struct deduce_vec2<Utils::Matrix<T, 3, 3>, Utils::Vector<U, 3>, 3> {
  using type = Utils::Vector<std::common_type_t<T, U>, 3>;
};

template <typename T, typename U>
struct deduce_vec2<Utils::Matrix<T, 4, 4>, Utils::Vector<U, 4>, 4> {
  using type = Utils::Vector<std::common_type_t<T, U>, 4>;
};

template <typename T, typename U>
struct deduce_vec2<Utils::Matrix<T, 2, 3>, Utils::Vector<U, 3>, 2> {
  using type = Utils::Vector<std::common_type_t<T, U>, 2>;
};

template <typename T, typename U>
struct deduce_mat2<Utils::Matrix<T, 3, 3>, Utils::Matrix<U, 3, 3>, 3, 3> {
  using type = Utils::Matrix<std::common_type_t<T, U>, 3, 3>;
};

} // namespace qvm
} // namespace boost

#endif
