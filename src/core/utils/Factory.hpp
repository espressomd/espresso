/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#include <memory>

namespace Utils {

/**
 * @brief Factory template.
 *
 * Can be used to construct registered classes by name.
 * One registry per base type (T). To get a new one,
 * use new type ( struct NewT : public T {}; ).
 */
template<class T>
class Factory {
 public:
  /** The returned pointer type */
  typedef std::unique_ptr<T> pointer_type;
  /** Type of the constructor functions */
  typedef std::function<T* ()> builder_type;

  /** Default constructor function, which just
   * calls new with the Derived type.
   */
  template<class Derived>
  static T* builder() {
    static_assert(std::is_base_of<T, Derived>::value,
                  "Class to build needs to be a subclass of the class the factory is for.");
    static_assert(std::is_default_constructible<Derived>::value,
                  "Derived class needs to be default constructible.");
    return new Derived();
  }

  static pointer_type make(const std::string &name) {
    if (m_map.find(name) == m_map.end()) {
      throw std::domain_error("Class '" + name + "' not found.");
    }

    if (m_map[name]) {
      return pointer_type(m_map[name]());
    } else {
      throw std::out_of_range("Invalid function pointer");
    }
  }

  static bool has_builder(const std::string &name) {
    return not (m_map.find(name) == m_map.end());
  }
  
  static void register_new(const std::string &name, const builder_type &b)  {
    m_map[name] = b;
  }

  template<typename Derived>
  static void register_new(const std::string &name) {
    register_new(name, builder<Derived>);
  }
  
 private:
  static std::map<std::string, builder_type> m_map;
};

template<class T>
std::map<std::string, typename Factory<T>::builder_type> Factory<T>::m_map;

  template<class T>
  auto factory_make(const std::string &name) -> typename Factory<T>::pointer_type {
    return Factory<T>::make(name);
  }

} /* namespace Utils */

#endif
