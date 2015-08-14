#ifndef __OBJECTCONTAINER_HPP
#define __OBJECTCONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class. 
    T needs to be copy-constructible.
 */

#include <map>

template<class T>
class ObjectConstainer : public std::map<int, T> {
public:
  ObjectContainer() : m_next_id(0) {};
  int add(const &T c) {
    insert(std::pair<int, T>(m_next_id, c));
    return m_next_id++;
  }
  void remove(int i) {
    erase(i);
  }
  private:
    int m_next_id;
  }
};

#endif
