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

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#include "ScriptInterface.hpp"

#include "core/errorhandling.hpp"
#include "core/utils/parallel/InstanceCallback.hpp"
#include "core/utils/parallel/ParallelObject.hpp"

#include <iostream>

namespace ScriptInterface {

namespace {
template <typename T>
class ParallelScriptInterfaceSlave : public Communication::InstanceCallback {
protected:
  enum class CallbackAction { SET_PARAMETER, SET_PARAMETERS, DELETE };
  T m_p;

private:
  void mpi_slave(int action, int) override {
    switch (CallbackAction(action)) {
    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);
      m_p.set_parameter(d.first, d.second);
    } break;
    case CallbackAction::SET_PARAMETERS: {
      std::map<std::string, Variant> parameters;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), parameters,
                            0);
      m_p.set_parameters(parameters);
    } break;
    case CallbackAction::DELETE: {
      delete this;
      break;
    }
    }
  }
};
}

template <typename T>
class ParallelScriptInterface
    : public ScriptInterfaceBase,
      public ParallelScriptInterfaceSlave<T>,
      public Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave<T>> {
public:
  using ParallelScriptInterfaceSlave<T>::m_p;
  using ParallelScriptInterfaceSlave<T>::call;
  using typename ParallelScriptInterfaceSlave<T>::CallbackAction;

  const std::string name() const override { return m_p.name(); }

  void set_parameter(const std::string &name, const Variant &value) override {
    std::cout << __PRETTY_FUNCTION__ << " name = " << name << std::endl;
    auto d = std::make_pair(name, value);

    call(static_cast<int>(CallbackAction::SET_PARAMETER), 0);

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);
    m_p.set_parameter(name, value);
  }

  void
  set_parameters(const std::map<std::string, Variant> &parameters) override {
    call(static_cast<int>(CallbackAction::SET_PARAMETERS), 0);

    /* Copy parameters into a non-const buffer, needed by boost::mpi */
    std::map<std::string, Variant> p(parameters);
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);

    m_p.set_parameters(p);
  }

  std::map<std::string, Parameter> all_parameters() const override {
    return m_p.all_parameters();
  }

  std::map<std::string, Variant> get_parameters() const override {
    return m_p.get_parameters();
  }
};

} /* namespace ScriptInterface */

#endif
