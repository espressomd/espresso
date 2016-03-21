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

#ifndef __UTILS_PARALLEL_FACTORY_HPP
#define __UTILS_PARALLEL_FACTORY_HPP

/** For is_base_of */
#include <type_traits>
/** For function */
#include <functional>
/** For unique_ptr */
#include <memory>
/** For ParallelObject */
#include "ParallelObject.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

namespace Utils {

/** Implements the same concept as Utils::Factory, but creates
    the instance on all nodes.
*/

template<class T>
class ParallelFactory {
 public:
  static_assert(std::is_base_of<ParallelObject, T>::value, "Class needs to be an ParallelObject for use with the parallel Factory.");
  typedef std::unique_ptr<T> pointer_type;
  typedef std::function<T *()> Builder;

  template<class Derived>
  static T *builder() {
    static_assert(std::is_base_of<T, Derived>::value,
                  "Class to build needs to be a subclass of the class the factory is for.");    
    return new Derived();
  }
  
  static pointer_type make(const std::string &name) {
    mpi_call(ParallelFactory<T>::mpi_slave, name.size(), 0);

    MPI_Bcast(const_cast<char *>(&(*name.begin())), name.size(), MPI_CHAR, 0, comm_cart);
    
    return pointer_type(do_make(name));
  }

  static void mpi_slave(int name_size, int) {
    std::string name;

    name.resize(name_size);
    MPI_Bcast(&(*name.begin()), name_size, MPI_CHAR, 0, comm_cart);
    
    try {      
      T *p = do_make(name);
      assert(p != nullptr);
    } catch(std::exception &e) {
      runtimeErrorMsg() << e.what();
    }
  }

  static void register_new(const std::string &name, const Builder &b) {
    static bool callback_added = false;
    /** Make sure that we have a callback registered as soon as
        there are classes that could be constructed. */
    if(!callback_added) {
      mpi_add_callback(mpi_slave);
      callback_added = true;
    }

    m_map[name] = b;
  }
  
  static bool has_builder(const std::string &name) {
    return not (m_map.find(name) == m_map.end());
  }
  
 private:
  static T *do_make(std::string name) {
    if (m_map.find(name) == m_map.end()) {
      throw std::domain_error("Class '" + name + "' not found.");
    }

    if (m_map[name]) {
      return m_map[name]();
    } else {
      throw std::out_of_range("Invalid function pointer");
    }
  }
  static std::map<std::string, Builder> m_map;
};

template<class T>
std::map<std::string, typename ParallelFactory<T>::Builder> ParallelFactory<T>::m_map;

}

#endif






