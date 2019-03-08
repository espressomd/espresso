/*
Copyright (C) 2010-2018 The ESPResSo project

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
#include "ParallelScriptInterfaceSlave.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

namespace ScriptInterface {

VariantMap ParallelScriptInterfaceSlave::bcast_variant_map() const {
  VariantMap ret;
  boost::mpi::broadcast(m_cb->comm(), ret, 0);

  /* If the parameter is a object we have to translate it first to a
     local id.
  */
  translate_id(ret);

  return ret;
}

ParallelScriptInterfaceSlave::ParallelScriptInterfaceSlave()
    : m_callback_id(m_cb, [this](CallbackAction a) { mpi_slave(a); }) {}

std::map<ObjectId, ObjectId> &
ParallelScriptInterfaceSlave::get_translation_table() {
  static std::map<ObjectId, ObjectId> m_translation_table;

  return m_translation_table;
}

void ParallelScriptInterfaceSlave::mpi_slave(CallbackAction action) {
  switch (action) {
  case CallbackAction::NEW: {
    std::pair<ObjectId, std::string> what;
    boost::mpi::broadcast(m_cb->comm(), what, 0);

    m_p = ScriptInterfaceBase::make_shared(
        what.second, ScriptInterfaceBase::CreationPolicy::LOCAL);

    get_translation_table()[what.first] = m_p->id();

    break;
  }
  case CallbackAction::CONSTRUCT: {
    auto const parameters = bcast_variant_map();

    m_p->construct(parameters);
    break;
  }
  case CallbackAction::SET_PARAMETER: {
    std::pair<std::string, Variant> d;
    boost::mpi::broadcast(m_cb->comm(), d, 0);

    /* If the parameter is a object we have to translate it first to a
       local id.
    */
    translate_id(d.second);
    m_p->set_parameter(d.first, d.second);
    break;
  }
  case CallbackAction::SET_PARAMETERS: {
    auto parameters = bcast_variant_map();

    m_p->set_parameters(parameters);

    break;
  }
  case CallbackAction::CALL_METHOD: {
    /* Name of the method and parameters */
    std::pair<std::string, VariantMap> d;

    /* Broadcast method name and parameters */
    boost::mpi::broadcast(m_cb->comm(), d, 0);

    translate_id(d.second);

    /* Forward to the local instance. */
    m_p->call_method(d.first, d.second);

    break;
  }
  case CallbackAction::DELETE: {
    delete this;
    break;
  }
  }
}

Communication::MpiCallbacks *ParallelScriptInterfaceSlave::m_cb = nullptr;

} /* namespace ScriptInterface */
