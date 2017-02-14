/*
  Copyright (C) 2017 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef UTILS_TENSOR_HPP
#define UTILS_TENSOR_HPP

#include <algorithm>
#include <array>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace Utils {

namespace {
/**
 * @brief Helper function needed because std::rbegin is not in c++11.
 */
template <typename Iterator>
auto make_reverse_iterator(Iterator it) -> std::reverse_iterator<Iterator> {
  return std::reverse_iterator<Iterator>(it);
}
}

template <typename T> class Tensor {
public:
  /* Concept Container requirements */
  using size_type = typename std::vector<T>::size_type;
  using difference_type = typename std::vector<T>::difference_type;
  using value_type = T;
  using referece = T &;
  using const_reference = const T &;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  /**
   * @brief Construct an empty Tensor.
   */
  Tensor() {}

  /**
   * @brief Construct a Tensor with given extents.
   */
  template <typename Extents> explicit Tensor(Extents extents) {
    resize(extents);
  }

  /**
   * @brief Construct a Tensor with given extents,
   * and fill the elements with value.
   */
  template <typename Extents> Tensor(Extents extents, T const &value) {
    resize(extents);

    std::fill(begin(), end(), value);
  }

  template <typename Indices> size_type linear_index(Indices indices) const {
    return std::inner_product(std::begin(indices), std::end(indices),
                              m_strides.begin(), 0);
  }

  /**
   * @brief Iterator pointing to the first element.
   */
  iterator begin() { return m_data.begin(); }

  /**
   * @brief Iterator pointing one past the last element.
   */
  iterator end() { return m_data.end(); }

  /**
   * @brief Constant iterator pointing to the first element.
   */
  const_iterator cbegin() const { return m_data.cbegin(); }

  /**
   * @brief Constant iterator pointing one past the last element.
   */
  const_iterator cend() const { return m_data.cend(); }

  /**
   * @brief Iterator to a specific element. If Indices is
   * shorter than the tensor rank, the remaining indices
   * are assumed to be 0. E.g. if you have a tensor of
   * rank 2, begin({1}) points to the element (1,0),
   * the beginning of the 2nd row.
   */
  template <typename Indices> iterator begin(Indices indices) {
    return begin() + linear_index(indices);
  }

  /**
   * @brief Iterator one past a specific element. If Indices is
   * shorter than the tensor rank, the iterator points to the
   * next row of the last index that was given: If you have
   * a tensor of rank 2, end({1}) points to the beginning
   * of row 3, such that the range [begin{1}, end{1}) spans
   * row 2.
   */
  template <typename Indices> iterator end(Indices indices) {
    std::vector<size_type> next(indices);
    next.back()++;

    return begin(next);
  }

  iterator end(std::initializer_list<size_type> indices) {
    return end<std::initializer_list<size_type>>(indices);
  }

  iterator begin(std::initializer_list<size_type> indices) {
    return begin<std::initializer_list<size_type>>(indices);
  }

  /**
   * @brief Access specific element. If Indices is shorter
   * thank the rank, the remaining indices are assumed to
   * be 0.
   */
  template <typename Indices> T &operator()(Indices const &indices) {
    return *(begin(indices));
  }

  template <typename Indices>
  T const &operator()(Indices const &indices) const {
    return *(begin(indices));
  }

  T const *data() const { return m_data.data(); }
  size_type size() const { return m_data.size(); }
  size_type max_size() const { return m_data.max_size(); }
  bool empty() const { return m_data.empty(); }

  /**
   * @brief Compare Tensors for equality.
   *
   * Two Tensors are equal iff they have the same dimensions,
   * and all the elements are equal.
   */
  bool operator==(Tensor const &rhs) const {
    return (m_extents == rhs.m_extents) &&
           std::equal(begin(), end(), rhs.begin());
  }

  /**
   * @brief Compare Tensors for inequality.
   *
   * Returns true iff the Tensors are not equal.
   */
  bool operator!=(Tensor const &rhs) const { return !(this->operator==(rhs)); }

  T &operator()(std::initializer_list<size_type> indices) {
    return *(begin(indices));
  }

  /**
   * @brief Swap two Tensors.
   */
  void swap(Tensor const &rhs) {
    /* Tensor is movable, so std::swap is fine. */
    std::swap(*this, rhs);
  }

  /**
   * @brief Change the rank and extents of the Tensor.
   *
   * Examples:
   * resize({}) makes the Tensor a scalar (rank 0),
   * resize({5}) a vector of length 5,
   * resize({3,2}) a 3 x 2 matrix and so on.
   *
   * May reallocate and invalidates all iterators.
   * The values of the elements are undefiened after
   * this operation.
   */
  template <typename Size> void resize(Size new_size) {
    auto const rank = std::distance(std::begin(new_size), std::end(new_size));
    m_strides.resize(rank);
    m_extents.resize(rank);

    std::copy(std::begin(new_size), std::end(new_size), m_extents.begin());

    auto const total_size =
        std::accumulate(std::begin(new_size), std::end(new_size), 1,
                        std::multiplies<size_type>());
    m_data.resize(total_size);

    size_type stride = 1;
    auto jt = make_reverse_iterator(m_strides.end());
    for (auto it = make_reverse_iterator(new_size.end());
         it != make_reverse_iterator(new_size.begin()); it++) {
      *jt = stride;
      ++jt;
      stride *= *it;
    }
  }

  void resize(std::initializer_list<size_type> new_size) {
    resize<std::initializer_list<size_type>>(new_size);
  }

  /**
   * @brief Reduce over a tensor axis.
   *
   * This function reduces the tensor rank by applying a binary op
   * along an axis. E.g. with std::plus<T>(), the elements along
   * the axis are summed up. This is equivalent to the sum() function
   * from numpy, but with arbitrary operations.
   *
   * The tensor rank is unchainged by this, but the extent for the axis
   * is set to 1.
   * The storage size of the object can change bis this function,
   * and reallocation may occur. It invalidates all iterators
   * into the Tensor.
   */
  template <typename BinaryOp>
  Tensor &reduce(size_type axis = 0, BinaryOp op = BinaryOp()) {
    if (rank() <= axis)
      return *this;

    auto const stride =
        (axis == 0) ? m_extents[0] * m_strides[0] : strides()[axis - 1];
    /* sum up in the first element for each row of the axis */
    for (auto it = begin(); it != end(); it += stride) {
      for (size_type i = 1; i < m_extents[axis]; i++) {
        std::transform(it + i * m_strides[axis], it + (i + 1) * m_strides[axis],
                       it, it, op);
      }
    }
    /* copy the elements to front, if axis == 0,
    * the new values are already at the front,
    * and there is nothing to copy. */
    if (axis > 0) {
      for (auto it = begin(), jt = begin(); it != end();
           it += stride, jt += m_strides[axis])
        std::copy_n(it, m_strides[axis], jt);
    }

    auto new_size = extents();
    new_size[axis] = 1;
    resize(new_size);
    return *this;
  }

  /**
   * @brief Sum over a tensor axis.
   *
   * This function reduces the tensor rank by summing
   * along an axis. This is equivalent to the sum() function
   * from numpy.
   *
   * The tensor rank is unchainged by this, but the extent for the axis
   * is set to 1.
   * The storage size of the object can change bis this function,
   * and reallocation may occur. It invalidates all iterators
   * into the Tensor.
   */
  Tensor &sum(size_type axis = 0) { return reduce(axis, std::plus<T>()); }

  /**
   * @brief the tensor rank.
   *
   * The dimension of the index set.
   */
  size_type rank() const { return m_extents.size(); }

  std::vector<size_type> const &strides() const { return m_strides; }
  std::vector<size_type> const &extents() const { return m_extents; }

private:
  std::vector<size_type> m_strides;
  std::vector<size_type> m_extents;
  std::vector<T> m_data;
};

} /* namespace Utils */

#endif
