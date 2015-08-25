#ifndef __OBJECTCONTAINER_HPP
#define __OBJECTCONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class. 
*/

#include <map>
#include <set>

template<class T>
class ObjectContainer {
public:
  ObjectContainer() {
    m_free_indices.insert(0);
    m_free_indices.insert(1);
  }
  int add(T* c) {
    const int ind = get_index();
    m_container.insert(std::pair<int, T*>(ind, c));
    return ind;
  }
  void remove(int i) {
    m_container.erase(i);
    m_free_indices.insert(i);
  }
  T* operator[](int i) { return m_container[i]; }
private:
  std::map<int, T*> m_container;
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
