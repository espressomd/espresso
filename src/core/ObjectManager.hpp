#ifndef __OBJECTMANAGER_HPP
#define __OBJECTMANAGER_HPP

#include "utils/ParallelFactory.hpp"
#include "ObjectContainer.hpp"

template<class T>
class ObjectManager {
public:
  typedef typename ObjectContainer<T>::iterator iterator;
  int add(std::string name) {
    T *p =  Utils::ParallelFactory<T>::make(name);
    p->call_slaves();
    const int id = m_objects.add(p);
    m_names[id] = name;
    return id;
  }
  void remove(int i) {
    delete m_objects[i];
    m_objects.remove(i);
    m_names.erase(i);
  }
  T* operator[](int i) { return m_objects[i]; }
  iterator begin() {
    return m_objects.begin();
  };
  iterator end() {
    return m_objects.end();
  };

  const std::string& name(int i) { return m_names[i]; }
private:
  std::map<int, std::string> m_names;
  ObjectContainer<T> m_objects;
};

#endif
