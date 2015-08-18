#ifndef __OBJECTCONTAINER_HPP
#define __OBJECTCONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class. 
    T needs to be copy-constructible.
 */

#include <map>
#include <set>

template<class T>
class ObjectContainer : public std::map<int, T> {
public:
  ObjectContainer() {
    m_free_indices.insert(0);
    m_free_indices.insert(1);
  }
  int add(const T& c) {
    const int ind = get_index();
    std::map<int, T>::insert(std::pair<int, T>(ind, c));
    return ind;
  }
  void remove(int i) {
    std::map<int, T>::erase(i);
    m_free_indices.insert(i);
  }
private:
  std::set<int> m_free_indices;
  int get_index() {
    /** Get lowest free index */
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
