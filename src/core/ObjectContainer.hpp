#ifndef __OBJECTCONTAINER_HPP
#define __OBJECTCONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class. 
*/

#include <unordered_map>
#include <set>
#include <cassert>

/**
 * \brief Container for objects that are identified by a numeric id.
 * The Container keeps a copy of the objects added.
 */
template<class T>
class ObjectContainer {
public:
  typedef typename std::unordered_map<int, T>::iterator iterator;
  ObjectContainer() {
    m_free_indices.insert(0);
    m_free_indices.insert(1);
  }
  /**
   * @brief Copy c into the container.
   *
   * Asign a free id to c and copy it into the container.
   *
   * @param c The object to add.
   */
  int add(T &c) {
    const int ind = get_index();
    m_container.insert(std::pair<int, T>(ind, c));
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
  void remove(int i) {
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
  T &operator[](int i) { return m_container[i]; }

  /**
   * @brief Get iterator to beginning of the container.
   */
  iterator begin() {
    return m_container.begin();
  }

  /**
   * @brief Get iterator to end of the container.
   */
  iterator end() {
    return m_container.end();
  }
  
private:
  /** Data storage */
  std::unordered_map<int, T> m_container;
  /** Set for free index bookkeeping */
  std::set<int> m_free_indices;
  /** @brief Get the next free index.
   *
   * This function gets the lowest free index by
   * poping the first element from the ordered set
   * of free indices.
   * If there is only 1 or less elements in the free
   * index set, we reused all indices that were ever
   * freeed and we add a new one at the end of the set.
   *
   * @return Free index.
   */
  int get_index() {
    /** Get lowest free index */
    /** If we don't have an free index, sth went wrong. */
    assert(!m_free_indices.empty());
    const int index = *m_free_indices.begin();
    /** and remove it from the list */
    m_free_indices.erase(index);

    /** If there is only on left, it is the highest ever seen, so we can savely add +1 */
    if(m_free_indices.size() == 1) {
      m_free_indices.insert(*(--m_free_indices.end())+1);
    }

    return index;
  }
};

#endif
