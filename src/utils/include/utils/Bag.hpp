/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef ESPRESSO_UTILS_BAG_HPP
#define ESPRESSO_UTILS_BAG_HPP

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace Utils {
/**
 * @brief Bag of elements.
 *
 * A bag is a container in which the elements do not have
 * a fixed order. It can be considered an unordered variant
 * of vector, and implements the Container named requirement
 * specified in C++11.
 *
 * Elements in the container do not have a stable position and
 * removing elements can change the order of the other elements.
 * The looser contract (compared to a vector) allows removing any
 * element in the container in constant time.
 *
 * @tparam T Element type, needs to be Swappable.
 */
template <class T> class Bag {
  static_assert(std::is_swappable_v<T>);

  /** Storage backend */
  using storage_type = std::vector<T>;

public:
  using value_type = T;
  using iterator = T *;
  using const_iterator = const T *;
  using pointer = T *;
  using reference = T &;

  /**
   * @brief Construct an empty container.
   */
  Bag() = default;

private:
  /** Underlying storage of the container */
  storage_type m_storage;

  friend boost::serialization::access;
  /**
   * @brief Serialize the container.
   *
   * Serialization requires T to be serializable.
   */
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &m_storage;
  }

public:
  iterator begin() { return m_storage.data(); }
  iterator end() { return m_storage.data() + size(); }
  const_iterator begin() const { return m_storage.data(); }
  const_iterator end() const { return m_storage.data() + size(); }

  /**
   * @brief Number of elements in the container.
   */
  std::size_t size() const { return m_storage.size(); }

  /**
   * @brief Is the container empty?
   * @return True if there are no elements.
   */
  bool empty() const { return m_storage.empty(); }

  /**
   * @brief Capacity of the container.
   *
   * Number of elements the container can at least hold
   * without reallocating.
   */
  std::size_t capacity() const { return m_storage.capacity(); }

  /**
   * @brief Maximum number of elements the container can hold.
   */
  std::size_t max_size() const { return m_storage.max_size(); }

  /**
   * @brief Reserve storage.
   *
   * Increase capacity to at least the specified value.
   *
   * @param new_capacity New minimum capacity.
   */
  void reserve(std::size_t new_capacity) { m_storage.reserve(new_capacity); }

  /**
   * @brief Resize container.
   *
   * Newly added Ts are default-initialized.
   * If the new size is larger than the capacity, all
   * iterators into the container are invalidated.
   *
   * @param new_size Size to resize to.
   */
  void resize(std::size_t new_size) { m_storage.resize(new_size); }

  /**
   * @brief Remove all elements form container.
   */
  void clear() { m_storage.clear(); }

  /**
   * @brief Insert an element into the container.
   *
   * If before the call size() >= capacity(),
   * this may reallocate, in which case all
   * iterators into the container are invalidated.
   * Otherwise only the end iterator is invalidated.
   *
   * @param v Element to add.
   * @return Reference to the added element.
   */
  T &insert(T const &v) {
    m_storage.push_back(v);

    return m_storage.back();
  }
  /** @overload */
  T &insert(T &&v) {
    m_storage.push_back(std::move(v));

    return m_storage.back();
  }

  /**
   * @brief Remove element from the list.
   *
   * @param it Iterator pointing to the element to remove.
   * @return An iterator past the element that was removed.
   */
  iterator erase(iterator it) {
    *it = std::move(m_storage.back());

    m_storage.pop_back();

    return it;
  }

  /**
   * @brief Swap two Bags.
   *
   * Efficiently swap to bags by swapping
   * their contained storage.
   */
  friend void swap(Bag &lhs, Bag &rhs) {
    using std::swap;

    swap(lhs.m_storage, rhs.m_storage);
  }
};
} // namespace Utils
#endif
