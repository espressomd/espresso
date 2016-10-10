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

#ifndef SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_SLAVE_HPP
#define SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_SLAVE_HPP

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include "core/utils/parallel/InstanceCallback.hpp"
#include "core/utils/parallel/ParallelObject.hpp"

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
public:
  enum class CallbackAction {
    SET_ID,
    SET_PARAMETER,
    SET_PARAMETERS,
    CALL_METHOD,
    DELETE
  };

protected:
  friend Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave<T>>;
  ParallelScriptInterfaceSlave() : m_p(ScriptInterfaceBase::make_shared<T>()) {}

public:
  std::shared_ptr<T> m_p;

private:
  void mpi_slave(int action, int id) override {
    switch (CallbackAction(action)) {
    case CallbackAction::SET_ID:
      get_translation_table()[id] = m_p->id();
      break;

    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

      /* If the parameter is a object we have to tranlate it first to a
         local id.
      */
      if (m_p->valid_parameters()[d.first].type() == ParameterType::OBJECT) {
        const int global_id = boost::get<OId>(d.second).id;
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

      auto const valid_parameters = m_p->valid_parameters();

      for (auto &it : parameters) {
        if (valid_parameters.at(it.first).type() == ParameterType::OBJECT) {
          const int global_id = boost::get<OId>(it.second).id;
          it.second = translate_id(global_id);
        }
      }

      m_p->set_parameters(parameters);

      break;
    }
    case CallbackAction::CALL_METHOD: {
      /* Name of the method and para// meters */
      // std::pair<std::string, VariantMap> d;

      // /* Broadcast method name and parameters */
      // boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

      // /* Forward to the local instance. */
      // m_p->call_method(d.first, d.second);

      break;
    }
    case CallbackAction::DELETE: {
      delete this;
      break;
    }
    }
  }
};

} /* namespace ScriptInterface */

#endif
