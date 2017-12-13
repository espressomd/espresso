#include "ParallelScriptInterfaceSlave.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include "core/utils/Factory.hpp"

namespace ScriptInterface {

ParallelScriptInterfaceSlave::ParallelScriptInterfaceSlave() {
  m_cb->add(Communication::MpiCallbacks::function_type(
      [this](int a, int) { mpi_slave(a, 0); }));
}

std::map<ObjectId, ObjectId> &
ParallelScriptInterfaceSlave::get_translation_table() {
  static std::map<ObjectId, ObjectId> m_translation_table;

  return m_translation_table;
}

void ParallelScriptInterfaceSlave::mpi_slave(int action, int = 0) {
  switch (CallbackAction(action)) {
  case CallbackAction::CREATE: {
    std::pair<ObjectId, std::string> what;
    boost::mpi::broadcast(m_cb->comm(), what, 0);

    m_p = ScriptInterfaceBase::make_shared(what.second);

    get_translation_table()[what.first] = m_p->id();

    break;
  }
  case CallbackAction::SET_PARAMETER: {
    std::pair<std::string, Variant> d;
    boost::mpi::broadcast(m_cb->comm(), d, 0);

    /* If the parameter is a object we have to tranlate it first to a
       local id.
    */
    translate_id(d.second);
    m_p->set_parameter(d.first, d.second);
    break;
  }
  case CallbackAction::SET_PARAMETERS: {
    std::map<std::string, Variant> parameters;
    boost::mpi::broadcast(m_cb->comm(), parameters, 0);

    /* If the parameter is a object we have to tranlate it first to a
       local id.
    */
    for (auto &p : parameters) {
      translate_id(p.second);
    }

    m_p->set_parameters(parameters);

    break;
  }
  case CallbackAction::CALL_METHOD: {
    /* Name of the method and parameters */
    std::pair<std::string, VariantMap> d;

    /* Broadcast method name and parameters */
    boost::mpi::broadcast(m_cb->comm(), d, 0);

    for (auto &p : d.second) {
      translate_id(p.second);
    }

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
