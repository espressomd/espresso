#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <map>
#include <string>
#include <functional>
#ifdef HAVE_CXX11
#include <type_traits>
#endif

using namespace std;

/** Factory template */

template<class T>
class Factory {
  typedef std::function<T *()> Builder;
public:
  template<class Derived>
  static T *builder() {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<T, Derived>::value, "Class to build needs to be a subclass of the class the factory is for.");
#endif
    return new Derived();
  }

  static Factory &Instance() {
    static Factory *S = new Factory();
    return *S;
  }

  T *make(string name) {
    if(m_map.find(name) == m_map.end()) {
      throw std::string().append(name).append(std::string(" not found."));
    }

    if(m_map[name]) {
      return m_map[name]();
    }
    else {
      throw std::string("invalid function pointer");
    }
  }

  bool register_new(string name, Builder b)  {
    m_map[name] = b;
    return true;
  }
private:
  Factory() {};
  Factory(Factory &rhs) {};
  Factory &operator=(Factory &rhs) {};
  map<string, Builder> m_map;
};

#endif
