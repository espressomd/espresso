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

#include "ScriptInterfaceBase.hpp"

#include "ParallelScriptInterfaceSlave.hpp"
#include "core/errorhandling.hpp"

#include "utils/Factory.hpp"

namespace ScriptInterface {
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
      private Communication::InstanceCallback,
      public Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave<T>> {
public:
  using CallbackAction =
      typename ParallelScriptInterfaceSlave<T>::CallbackAction;

  std::shared_ptr<T> m_p;

  ParallelScriptInterface() : m_p(ScriptInterfaceBase::make_shared<T>()) {

    call(static_cast<int>(CallbackAction::SET_ID));

    ObjectId global_id{m_p->id()};

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), global_id, 0);
  }

  std::shared_ptr<ScriptInterfaceBase> get_underlying_object() const override {
    return std::static_pointer_cast<ScriptInterfaceBase>(m_p);
  };

  const std::string name() const override { return m_p->name(); }

  void set_parameter(const std::string &name, const Variant &value) override {
    std::pair<std::string, Variant> d(name, Variant());

    if (valid_parameters()[name].type() == ParameterType::OBJECTID) {
      d = std::make_pair(name, map_parallel_to_local_id(name, value));
    } else {
      d = std::make_pair(name, value);
    }

    call(static_cast<int>(CallbackAction::SET_PARAMETER), 0);

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

    m_p->set_parameter(d.first, d.second);
  }

  Variant map_local_to_parallel_id(std::string const &name,
                                   Variant const &value) const {
    if (boost::get<ObjectId>(value) != ObjectId()) {
      return obj_map.at(name)->id();
    } else {
      return value;
    }
  }

  Variant map_parallel_to_local_id(std::string const &name,
                                   Variant const &value) {
    const auto outer_id = boost::get<ObjectId>(value);

    auto so_ptr = get_instance(outer_id).lock();

    auto po_ptr =
        std::dynamic_pointer_cast<ParallelScriptInterfaceBase>(so_ptr);

    if (po_ptr != nullptr) {
      /* Store a pointer to the object to keep it alive */
      obj_map[name] = po_ptr;

      /* and return the id of the underlying object */
      auto underlying_object = po_ptr->get_underlying_object();

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

    /* Unwrap the object ids */
    for (auto &it : p) {
      if (valid_parameters()[it.first].type() == ParameterType::OBJECTID) {
        it.second = map_parallel_to_local_id(it.first, it.second);
      }
    }

    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);

    m_p->set_parameters(p);
  }

  std::map<std::string, Parameter> valid_parameters() const override {
    return m_p->valid_parameters();
  }

  std::map<std::string, Variant> get_parameters() const override {
    auto p = m_p->get_parameters();

    /* Wrap the object ids */
    for (auto &it : p) {
      if (valid_parameters()[it.first].type() == ParameterType::OBJECTID) {
        it.second = map_local_to_parallel_id(it.first, it.second);
      }
    }

    return p;
  }

  Variant call_method(const std::string &name,
                      const VariantMap &parameters) override {
    InstanceCallback::call(static_cast<int>(CallbackAction::CALL_METHOD), 0);
    VariantMap p = parameters;

    /* Unwrap the object ids */
    for (auto &it : p) {
      if (it.second.which() == static_cast<int>(ParameterType::OBJECTID)) {
        it.second = map_parallel_to_local_id(it.first, it.second);
      }
    }

    auto d = std::make_pair(name, p);
    /* Broadcast method name and parameters */
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

    return m_p->call_method(name, p);
  }

private:
  std::map<std::string, std::shared_ptr<ParallelScriptInterfaceBase>> obj_map;

public:
  static void register_new(std::string const &name) {
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
