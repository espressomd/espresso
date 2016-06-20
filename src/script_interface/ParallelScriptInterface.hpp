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

#ifndef SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_HPP
#define SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_HPP

#include <utility>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "errorhandling.hpp"

#include "ScriptInterface.hpp"
#include "InstanceCallback.hpp"
#include "utils/ParallelFactory.hpp"

namespace ScriptInterface {

template<typename Factory = Utils::Factory<ScriptInterfaceBase> >
class ParallelScriptInterface : public ScriptInterfaceBase, public Communication::InstanceCallback {
 public:
  void init(const std::string &name) {
    std::string n_name(name);
    
    call(INIT, 0);
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), n_name, 0);
    do_init(name);
  }
  
  const std::string name() const override {
    return m_p->name();
  }

  void set_parameter(const std::string &name, const Variant &value) override {
    auto d = std::make_pair(name, value);

    call(SET_PARAMETER, 0);
    
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);
    m_p->set_parameter(name, value);
  }

  void set_parameters(const std::map<std::string, Variant> &parameters) override {
    call(SET_PARAMETERS, 0);

    /* Copy parameters into a non-const buffer, needed by boost::mpi */
    std::map<std::string, Variant> p(parameters);    
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);
  }
  
  std::map<std::string, Parameter> all_parameters() const override {
    return m_p->all_parameters();
  }
  
  std::map<std::string, Variant> get_parameters() const override {
    return m_p->get_parameters();
  }
  
 private:
  enum CallbackActions{ INIT, SET_PARAMETER, SET_PARAMETERS };
  typename Factory::pointer_type m_p;

  void do_init(const std::string &name) {
    try {
      m_p = Factory::make(name);
    } catch(std::exception &e) {
      runtimeErrorMsg() << e.what();
    }
  }
  
  void mpi_slave(int action, int) override {
    switch(action) {
      case INIT:
        {
          std::string name;
          boost::mpi::broadcast(Communication::mpiCallbacks().comm(), name, 0);
          do_init(name);
        }
        break;
      case SET_PARAMETER:
        {
          std::pair<std::string, Variant> d;
          boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);
          m_p->set_parameter(d.first, d.second);
        }
        break;
      case SET_PARAMETERS:
        {
          std::map<std::string, Variant> parameters;
          boost::mpi::broadcast(Communication::mpiCallbacks().comm(), parameters, 0);
          m_p->set_parameters(parameters);
        }
        break;
    }
  }
};

} /* namespace ScriptInterface */

#endif
