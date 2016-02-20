#ifndef __OBJECTMANAGER_HPP
#define __OBJECTMANAGER_HPP

#include "utils/ParallelFactory.hpp"
#include "ObjectContainer.hpp"

#include <memory>

template<class T, class Factory = Utils::ParallelFactory<T> >
class ObjectManager {
 public:
  typedef typename ObjectContainer<std::shared_ptr<T>>::iterator iterator;
  typedef Factory factory_type;
  
  /** Construct and add a T */
  int add(std::string name) {
    std::shared_ptr<T> p = Factory::make(name);
    const int id = m_objects.add(p);
    m_names[id] = name;
    return id;
  }
  
  /** Add an externally constructed T */
  int add(std::shared_ptr<T> o) {    
    const int id = m_objects.add(o);

    return id;
  }
  
  void remove(int i) {
    /** remove calls the destructor of the shared_ptr, which in turn
        destoys the referenced object iff it is referenced nowhere else. */    
    m_objects.remove(i);
    m_names.erase(i);
  }
  std::shared_ptr<T> operator[](int i) { return m_objects[i]; }
  
  iterator begin() {
    return m_objects.begin();
  };
  iterator end() {
    return m_objects.end();
  };

  const std::string& name(int i) { return m_names[i]; }
 private:
  std::map<int, std::string> m_names;
  ObjectContainer<std::shared_ptr<T>> m_objects;
};

#endif
