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

#include <type_traits>
#include <functional>
#include <memory>
#include <exception>
#include <algorithm>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

#include "Factory.hpp"
#include "errorhandling.hpp"

#include "MpiCallbacks.hpp"

namespace Utils {

/**
 * @brief Factory for parallel object construction.
 *
 * Implements the same concept as Utils::Factory, but creates
 * the instance on all nodes.
*/
template<class T, /* Class the factory is for. */
         class Factory = typename Utils::Factory<T> /* The factory used to locally create the objects */
         >
class ParallelFactory {
 public:
  typedef std::shared_ptr<T> pointer_type;
  typedef Factory factory_type;

  ParallelFactory() {
    /** Register the callback for object creation. */
    register_callback();
  }
  
  ~ParallelFactory() {
    Communication::mpiCallbacks().remove(m_callback_id);
  }
  
  pointer_type make(std::string name) {
    Communication::mpiCallbacks().call(m_callback_id, MAKE, 0);

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), name, 0);

    auto p = do_make(name);
    assert(p != nullptr);

    /*
     * We set an emty deleter [](T *){}, which
     * let's us use types with private destructor.
     * This is ok because the objects should only
     * be deleted by the factory, and we always
     * keep a copy around here. So we use shared_ptr
     * only as a reference counter.
     */
    pointer_type shared(std::move(p), [](T *){});
    m_instances.add(shared);
    
    return shared;
  }

  /**
   * @brief Delete the object referenced by shared.
   * Shared is invalid after the call.
   * The object referenced by shared must have been
   * created by this instance of ParallelFactory.
   */
  void remove(pointer_type &shared) {
    /* Search for the instance. */
    auto it = std::find_if(m_instances.begin(), m_instances.end(),
                           [&shared](const std::pair<const int, pointer_type>  &j) {
                             return shared == j.second;
                           });

    if(it == m_instances.end()) {
      throw std::out_of_range("Instance not found.");
    }

    /* Extra reference check on the master.
     * if this fails, we can stay in a defined
     * state because no node has detroyed anything.
     */
    if(shared.use_count() != 2) {
      throw std::runtime_error(
          "Tried to delete an object that is still referenced.");
    }
    
    /* Reset the pointer from the caller to decrease
     *  the refcount. */
    shared.reset();

    /* Actually remove the object */
    Communication::mpiCallbacks().call(m_callback_id, DELETE, it->first);
    
    do_remove(it->first);
  }
      
 private:
  enum CallbackAction { MAKE, DELETE };
  
  T* do_make(const std::string &name) {
    auto builder = Factory::get_builder(name);
    return builder();
  }

  void do_remove(int id) {
    auto sp = m_instances[id];

    /* There should be no copies out if
     * the instance is to be deleted.
     * The two copies left are the one in
     * m_instances and sp.
     */
    if(sp.use_count() != 2) {
      throw std::runtime_error("Tried to delete an object that is still referenced.");
    }

    /* Delete the shared_ptr, does not delete the object because
     * the deleter is empty */
    m_instances.remove(id);
    
    /* Delete the object */
    delete sp.get();
  }
  /** Register the callback for object creation. */
  void register_callback() {
    using std::placeholders::_1;
    using std::placeholders::_2;
    Communication::MpiCallbacks::function_type f;

    /* Bind member function to this instance */
    auto cb = std::bind(&ParallelFactory::mpi_slave, this, _1, _2);

    m_callback_id = Communication::mpiCallbacks().add(cb);
  }

  friend Communication::MpiCallbacks;
  void mpi_slave(int action, int id) {
    switch(action) {
      case MAKE:
        {
          std::string name;

          boost::mpi::broadcast(Communication::mpiCallbacks().comm(), name, 0);    
    
          try {      
            auto p = do_make(name);
            assert(p != nullptr);
      
            pointer_type shared(std::move(p), [](T *){});  
            m_instances.add(shared);      
          } catch(std::exception &e) {
            runtimeErrorMsg() << e.what();
          }
          break;
        }
      case DELETE:
        {
          try {
            do_remove(id);
          } catch(std::exception &e) {
            runtimeErrorMsg() << e.what();
          }
          break;
        }
    }
  }
  
  int m_callback_id;
  Utils::NumeratedContainer<pointer_type> m_instances;
};

}

#endif
