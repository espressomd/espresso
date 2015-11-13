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

#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <map>
#include <string>
#include <functional>
#ifdef HAVE_CXX11
#include <type_traits>
#endif

namespace Utils {

/** Factory template */

template<class T>
class Factory {
 public:
  typedef std::function<T *()> Builder;
  template<class Derived>
  static T *builder() {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<T, Derived>::value,
                  "Class to build needs to be a subclass of the class the factory is for.");
#endif
    return new Derived();
  }

  static Factory &Instance() {
    static Factory *S = new Factory();
    return *S;
  }

  T *make(std::string name) {
    if (m_map.find(name) == m_map.end()) {
      throw std::string().append(name).append(std::string(" not found."));
    }

    if (m_map[name]) {
      return m_map[name]();
    } else {
      throw std::string("invalid function pointer");
    }
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
