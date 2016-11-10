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

#include "ScriptInterface.hpp"

#include "ParallelScriptInterfaceSlave.hpp"
#include "core/errorhandling.hpp"

namespace ScriptInterface {

class ParallelScriptInterface : public ScriptInterfaceBase,
                                private Communication::InstanceCallback {
public:
  using CallbackAction = ParallelScriptInterfaceSlave::CallbackAction;

  ParallelScriptInterface(std::string const &name);

  static void initialize();

  bool operator==(ParallelScriptInterface const &rhs);
  bool operator!=(ParallelScriptInterface const &rhs);

  std::shared_ptr<ScriptInterfaceBase> get_underlying_object() const {
    return std::static_pointer_cast<ScriptInterfaceBase>(m_p);
  };

  /* Script interface */
  const std::string name() const override { return m_p->name(); }
  void set_parameter(const std::string &name, const Variant &value) override;
  void
  set_parameters(const std::map<std::string, Variant> &parameters) override;
  std::map<std::string, Parameter> valid_parameters() const override {
    return m_p->valid_parameters();
  }

  std::map<std::string, Variant> get_parameters() const override {
    auto p = m_p->get_parameters();

    /* Wrap the object ids */
    for (auto &it : p) {
      if (is_objectid(it.second)) {
        it.second = map_local_to_parallel_id(it.first, it.second);
      }
    }

    return p;
  }

  /* Id mapping */
  Variant map_local_to_parallel_id(std::string const &name,
                                   Variant const &value) const {
    return obj_map.at(name);
  }

  Variant map_parallel_to_local_id(std::string const &name,
                                   Variant const &value) {
    const auto outer_id = boost::get<ObjectId>(value);

    auto so_ptr = get_instance(outer_id).lock();

    auto po_ptr = std::dynamic_pointer_cast<ParallelScriptInterface>(so_ptr);

    if (po_ptr != nullptr) {
      /* Store a pointer to the object */
      obj_map[name] = po_ptr->id();

      /* and return the id of the underlying object */
      auto underlying_object = po_ptr->get_underlying_object();

      return underlying_object->id();
    } else {
      throw std::runtime_error(
          "Parameters passed to Parallel entities must also be parallel.");
    }
  }

  Variant call_method(const std::string &name,
                      const VariantMap &parameters) override {
    InstanceCallback::call(static_cast<int>(CallbackAction::CALL_METHOD), 0);
    VariantMap p = parameters;

    /* Unwrap the object ids */
    for (auto &it : p) {
      if (it.second.which() == static_cast<int>(ParameterType::OBJECT)) {
        it.second = map_parallel_to_local_id(it.first, it.second);
      }
    }

    auto d = std::make_pair(name, p);
    /* Broadcast method name and parameters */
    boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

    return m_p->call_method(name, p);
  }

private:
  /* Data members */
  std::shared_ptr<ScriptInterfaceBase> m_p;
  std::map<std::string, ObjectId> obj_map;
};

} /* namespace ScriptInterface */

#endif
