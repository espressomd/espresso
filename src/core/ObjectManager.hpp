#ifndef __OBJECTMANAGER_HPP
#define __OBJECTMANAGER_HPP

#include "Factory.hpp"
#include "ObjectContainer.hpp"

#include <stdio.h>

template<class T>
class ObjectManager {
public:
  int add(std::string name) {
    T *p =  Factory<T>::Instance().make(name);
    return m_objects.add(p);
  }
  void remove(int i) {
    delete m_objects[i];
    m_objects.remove(i);
  }
  T* operator[](int i) { return m_objects[i]; }
private:
  ObjectContainer<T> m_objects;
};

#endif
