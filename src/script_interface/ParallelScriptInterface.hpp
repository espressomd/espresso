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
#include "utils/Factory.hpp"

namespace ScriptInterface {

class ParallelScriptInterfaceSlaveBase {
protected:
  static std::map<int, int> &get_translation_table() {
    static std::map<int, int> m_translation_table;

    return m_translation_table;
  }

  static int translate_id(int id) { return get_translation_table().at(id); }
};

template <typename T>
class ParallelScriptInterfaceSlave : public Communication::InstanceCallback,
                                     private ParallelScriptInterfaceSlaveBase {
protected:
  friend Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave<T>>;
  ParallelScriptInterfaceSlave() : m_p(ScriptInterfaceBase::make_shared<T>()) {
    std::cout << __PRETTY_FUNCTION__ << " m_p->id() =  " << m_p->id()
              << ", m_p = " << m_p.get() << std::endl;
  }

  enum class CallbackAction {
    SET_ID,
    SET_PARAMETER,
    SET_PARAMETERS,
    CALL_METHOD,
    DELETE
  };

public:
  std::shared_ptr<T> m_p;

private:
  void mpi_slave(int action, int id) override {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << std::endl;

    switch (CallbackAction(action)) {
    case CallbackAction::SET_ID:
      get_translation_table()[id] = m_p->id();
      std::cout << __PRETTY_FUNCTION__ << " mapping " << id << " -> "
                << translate_id(id) << std::endl;
      break;
    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

      /* If the parameter is a object we have to tranlate it first to a
         local id.
      */
      if (m_p->all_parameters()[d.first].type() == ParameterType::OBJECT) {
        const int global_id = boost::get<int>(d.second);

        std::cout << Communication::mpiCallbacks().comm().rank() << ": "
                  << __PRETTY_FUNCTION__ << " mapping parameter " << d.first
                  << ": " << global_id << std::endl;

        const int local_id = translate_id(global_id);
        m_p->set_parameter(d.first, local_id);
      } else {
        m_p->set_parameter(d.first, d.second);
      }

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
  virtual std::shared_ptr<ScriptInterfaceBase>
  get_underlying_object() const = 0;

  virtual bool operator==(ParallelScriptInterfaceBase const &rhs) {
    return this->get_underlying_object() == rhs.get_underlying_object();
  }

  virtual bool operator!=(ParallelScriptInterfaceBase const &rhs) {
    return !(*this == rhs);
  }
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
    call(static_cast<int>(CallbackAction::SET_ID), m_p->id());
  }

  std::shared_ptr<ScriptInterfaceBase> get_underlying_object() const override {
    return std::static_pointer_cast<ScriptInterfaceBase>(m_p);
  };

  const std::string name() const override { return m_p->name(); }

  void set_parameter(const std::string &name, const Variant &value) override {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << std::endl;

    auto d = std::make_pair(name, transform_object_parameter(name, value));
    call(static_cast<int>(CallbackAction::SET_PARAMETER), 0);

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

    m_p->set_parameter(d.first, d.second);
  }

  Variant transform_object_parameter(std::string const &name,
                                     Variant const &value) const {
    const int outer_id = boost::get<int>(value);

    auto so_ptr = get_instance(outer_id).lock();

    auto po_ptr =
        std::dynamic_pointer_cast<ParallelScriptInterfaceBase>(so_ptr);

    std::cout << __PRETTY_FUNCTION__ << " po_ptr = " << po_ptr.get()
              << std::endl;

    if (po_ptr != nullptr) {
      auto underlying_object = po_ptr->get_underlying_object();
      std::cout << __PRETTY_FUNCTION__ << " outer_id = " << outer_id
                << ", inner_id = " << po_ptr->get_underlying_object()->id()
                << std::endl;
      return underlying_object->id();
    } else {
      throw std::runtime_error(
          "Parameters passed to Parallel entities must also be parallel.");
    }
  }

  void
  set_parameters(const std::map<std::string, Variant> &parameters) override {
    call(static_cast<int>(CallbackAction::SET_PARAMETERS), 0);

    /* Copy parameters into a non-const buffer, needed by boost::mpi */
    std::map<std::string, Variant> p(parameters);
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);

    /* Unwrap the object ids */
    for (auto &it : p) {
      if (all_parameters()[it.first].type() == ParameterType::OBJECT) {
        it.second = transform_object_parameter(it.first, it.second);
      }
    }

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
    std::cout << __PRETTY_FUNCTION__ << " name = " << name << std::endl;

    /* Register with the factory */
    Utils::Factory<ScriptInterfaceBase>::register_new<
        ParallelScriptInterface<T>>(name);

    /* Register the object creation callback */
    Utils::Parallel::ParallelObject<
        ParallelScriptInterfaceSlave<T>>::register_callback();
  }
};

} /* namespace ScriptInterface */

#endif
