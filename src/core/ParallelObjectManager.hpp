#ifndef __PARALLEL_OBJECT_MANAGER_HPP
#define __PARALLEL_OBJECT_MANAGER_HPP

#ifdef HAVE_CXX11
#include <type_traits>
#endif

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

#include "ObjectManager.hpp"
#include "ScriptObject.hpp"

#include "communication.hpp"

template<class T>
class ParallelObjectManager {
 public:
  ParallelObjectManager() {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<ScriptObject, T>::value, "Type has to be subclass of ScriptObject");
    static_assert(std::is_default_constructible<T>::value, "Type has to be default constructible");
#endif
  }
    int add(std::string name) {
    const int id = m_om.add(name);
    
    return id;    
  }
    int remove(std::string name) {
    
    
  }
    
 private:
    static void add_slave(int id, int) {
    std::string name;
    
    m_om.add(name);
  }
    static void remove_slave(int id, int) {
    std::string name;
    
    m_om.remove(name);
  }
    ObjectManager<T> m_om;
  
}

#endif
