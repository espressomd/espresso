/*
  Copyright (C) 2015-2018 The ESPResSo project

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

#ifndef UTILS_ENUMERATED_CONTAINER_HPP
#define UTILS_ENUMERATED_CONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class.
 */

#include <cassert>
#include <set>
#include <unordered_map>

namespace Utils {

/**
 * @brief Container for objects that are identified by a numeric id.
 *
 * This class implements a container that holds T instances and references
 * them by an integral index. New elements get the lowest free index
 * and indexes are reused. The Container keeps a copy of the objects added.
 */
template <class T, typename index_type = int> class NumeratedContainer {
public:
  typedef typename std::unordered_map<index_type, T>::iterator iterator;
  typedef
      typename std::unordered_map<index_type, T>::const_iterator const_iterator;
  typedef typename std::unordered_map<index_type, T>::value_type value_type;

  NumeratedContainer() {
    m_free_indices.insert(0);
    m_free_indices.insert(1);
  }

  /**
   * @brief Construct from list of key-value pairs.
   * The keys have to be unique.
   */
  explicit NumeratedContainer(std::initializer_list<value_type> l)
      : NumeratedContainer() {
    for (auto const &e : l) {
      m_container.insert(e);
      /* Remove the index from the index set if it exists. */
      m_free_indices.erase(m_free_indices.find(e.first), m_free_indices.end());
    }

    /* Refill the index set */
    for (index_type it(0); m_free_indices.size() < 2; ++it) {
      if (m_container.find(it) == m_container.end()) {
        m_free_indices.insert(it);
      }
    }
  }

  /**
   * @brief Copy c into the container.
   *
   * Assign a free id to c and copy it into the container.
   *
   * @param c The object to add.
   */
  index_type add(const T &c) {
    const index_type ind = get_index();
    m_container.emplace(std::make_pair(ind, c));
    return ind;
  }

  /** @overload */
  index_type add(T &&c) {
    const index_type ind = get_index();
    m_container.emplace(std::make_pair(ind, std::move(c)));
    return ind;
  }

  /**
   * @brief Remove element from container
   *
   * Remove element i and add i to the free
   * indices.
   *
   * @param i The object to remove.
   */
  void remove(index_type i) {
    /** Check that the object actually exists */
    assert(m_container.find(i) != m_container.end());

    m_container.erase(i);
    m_free_indices.insert(i);
  }

  /**
   * @brief Get element from container
   *
   * @param i The object to get.
   * @return Reference to the object with id i.
   */
  T &operator[](index_type i) { return m_container.at(i); }

  /**
   * @brief Get element from container
   *
   * @param i The object to get.
   * @return Reference to the object with id i.
   */
  const T &operator[](index_type i) const { return m_container.at(i); }

  /**
   * @brief Get iterator to beginning of the container.
   */
  iterator begin() { return m_container.begin(); }

  /**
   * @brief Get a const iterator to beginning of the container.
   */
  const_iterator begin() const { return m_container.begin(); }

  /**
   * @brief Get an iterator to end of the container.
   */
  iterator end() { return m_container.end(); }

  /**
   * @brief Get a const iterator to end of the container.
   */
  const_iterator end() const { return m_container.end(); }

  /**
   * @brief find object by id.
   *
   * @param id id of the object to find.
   * @return Iterator to the object if found or end().
   */
  iterator find(int id) { return m_container.find(id); }

  /**
   * @brief find object by id.
   *
   * @param id id of the object to find.
   * @return Iterator to the object if found or end().
   */
  const_iterator find(int id) const { return m_container.find(id); }

  size_t size() const { return m_container.size(); }

private:
  /** Data storage */
  std::unordered_map<index_type, T> m_container;
  /** Set for free index bookkeeping */
  std::set<index_type> m_free_indices;

  /**
   * @brief Get the next free index.
   *
   * This function gets the lowest free index by
   * popping the first element from the ordered set
   * of free indices.
   * If there is only 1 or less elements in the free
   * index set, we reused all indices that were ever
   * freed and we add a new one at the end of the set.
   *
   * @return Free index.
   */
  index_type get_index() {
    /** Get lowest free index */
    /** If we don't have an free index, sth went wrong. */
    assert(!m_free_indices.empty());
    const index_type index = *m_free_indices.begin();
    /** and remove it from the list */
    m_free_indices.erase(index);

    /** If there is only on left, it is the highest ever seen, so we can safely
     * add +1 */
    if (m_free_indices.size() == 1) {
      m_free_indices.insert(*(--m_free_indices.end()) + 1);
    }

    return index;
  }
};
} // namespace Utils

#endif
