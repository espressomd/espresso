#ifndef __OBJECTMANAGER_HPP
#define __OBJECTMANAGER_HPP

#include "Factory.hpp"
#include "ObjectContainer.hpp"

#include <stdio.h>

template<class T>
class ObjectManager {
public:
  int add(std::string name) {
    T *p =  Utils::Factory<T>::Instance().make(name);
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
  const std::string& name(int i) { return m_names[i]; }
private:
  std::map<int, std::string> m_names;
  ObjectContainer<T> m_objects;
};

#endif
