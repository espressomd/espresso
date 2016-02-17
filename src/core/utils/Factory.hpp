/*
  Copyright (C) 2010,2011,2012,2013,2014,2015 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
  Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef __UTILS_FACTORY_HPP
#define __UTILS_FACTORY_HPP

#include <map>
#include <string>
#include <functional>
#include <type_traits>
#include <exception>

namespace Utils {

/** Factory template */

template<class T>
class Factory {
 public:
  typedef T* pointer_type;
  typedef std::function<pointer_type()> Builder;
  
  template<class Derived>
  static pointer_type builder() {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<T, Derived>::value,
                  "Class to build needs to be a subclass of the class the factory is for.");
#endif
    return pointer_type(new Derived());
  }

  static Factory &Instance() {
    static Factory *S = new Factory();
    return *S;
  }

  pointer_type make(std::string name) {
    if (m_map.find(name) == m_map.end()) {
      throw std::domain_error("Class '" + name + "' not found.");
    }

    if (m_map[name]) {
      return m_map[name]();
    } else {
      throw std::out_of_range("Invalid function pointer");
    }
  }

  bool has_builder(const std::string &name) {
    return not (m_map.find(name) == m_map.end());
  }
  
  bool register_new(const std::string &name, const Builder &b)  {
    m_map[name] = b;
    return true;
  }
  
 private:
  Factory() {}
  Factory(const Factory &rhs) {}
  Factory &operator=(const Factory &rhs) {}
  std::map<std::string, Builder> m_map;
};

}

#endif
