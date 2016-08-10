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
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include "ScriptInterface.hpp"

#include "core/errorhandling.hpp"
#include "core/utils/parallel/InstanceCallback.hpp"
#include "core/utils/parallel/ParallelObject.hpp"

#include <iostream>

namespace ScriptInterface {

template <typename T>
class ParallelScriptInterfaceSlave : public Communication::InstanceCallback {
protected:
  ParallelScriptInterfaceSlave() : m_p(new T()) {}

  enum class CallbackAction {
    SET_ID,
    SET_PARAMETER,
    SET_PARAMETERS,
    CALL_METHOD,
    DELETE
  };

  std::shared_ptr<T> m_p;

private:
  static std::map<int, int> &get_translation_table() {
    static std::map<int, int> m_translation_table;

    return m_translation_table;
  }

  static int translate_id(int id) { return get_translation_table().at(id); }

  void mpi_slave(int action, int id) override {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << std::endl;

    switch (CallbackAction(action)) {
    case CallbackAction::SET_ID:
      std::cout << __PRETTY_FUNCTION__ << " mapping " << id << " -> "
                << m_p->id() << std::endl;
      get_translation_table()[id] = m_p->id();
      break;
    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

      if (m_p->all_parameters()[d.first].type() == ParameterType::OBJECT) {
        std::cout << Communication::mpiCallbacks().comm().rank() << ": "
                  << __PRETTY_FUNCTION__ << " mapping parameter " << d.first
                  << ": " << boost::get<int>(d.second) << " -> ";
        std::cout << translate_id(boost::get<int>(d.second)) << std::endl;
      }

      m_p->set_parameter(d.first, d.second);

      break;
    }
    case CallbackAction::SET_PARAMETERS: {
      std::map<std::string, Variant> parameters;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), parameters,
                            0);

      m_p->set_parameters(parameters);

      break;
    }
    case CallbackAction::CALL_METHOD: {
      /* Name of the method to call */
      std::string name;
      /* Parameters for the method */
      VariantMap parameters;

      /* Broadcast method name and parameters */
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), name, 0);
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), parameters,
                            0);

      /* Forward to the local instance. */
      m_p->call_method(name, parameters);

      break;
    }
    case CallbackAction::DELETE: {
      delete this;
      break;
    }
    }
  }
};

class ParallelScriptInterfaceBase : public ScriptInterfaceBase {
public:
  virtual std::shared_ptr<ScriptInterfaceBase> const &
  get_underlying_object() const = 0;
};

template <typename T>
class ParallelScriptInterface
    : public ParallelScriptInterfaceBase,
      public ParallelScriptInterfaceSlave<T>,
      public Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave<T>> {
public:
  using ParallelScriptInterfaceSlave<T>::m_p;
  using ParallelScriptInterfaceSlave<T>::call;
  using typename ParallelScriptInterfaceSlave<T>::CallbackAction;

  ParallelScriptInterface() {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << std::endl;

    call(static_cast<int>(CallbackAction::SET_ID), m_p->id());
  }

  std::shared_ptr<ScriptInterfaceBase> const &
  get_underlying_object() const override {
    return m_p;
  };

  const std::string name() const override { return m_p->name(); }

  void set_parameter(const std::string &name, const Variant &value) override {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << std::endl;

    auto d = std::make_pair(name, value);

    call(static_cast<int>(CallbackAction::SET_PARAMETER), 0);

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);
    m_p->set_parameter(name, value);
  }

  void
  set_parameters(const std::map<std::string, Variant> &parameters) override {
    call(static_cast<int>(CallbackAction::SET_PARAMETERS), 0);

    /* Copy parameters into a non-const buffer, needed by boost::mpi */
    std::map<std::string, Variant> p(parameters);
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);

    m_p->set_parameters(p);
  }

  std::map<std::string, Parameter> all_parameters() const override {
    return m_p->all_parameters();
  }

  std::map<std::string, Variant> get_parameters() const override {
    return m_p->get_parameters();
  }

  void call_method(const std::string &name,
                   const VariantMap &parameters) override {
    call(static_cast<int>(CallbackAction::CALL_METHOD), 0);
  }

public:
  static void register_new(std::string const &name) {
    /* Register the object creation callback */
    Utils::Parallel::ParallelObject<
        ParallelScriptInterfaceSlave<T>>::register_callback();

    /* Register with the factory */
    Utils::Factory<ScriptInterfaceBase>::register_new<
        ParallelScriptInterface<T>>(name);
  }
};

} /* namespace ScriptInterface */

#endif
