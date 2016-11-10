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

#include "core/utils/Factory.hpp"
#include "core/utils/parallel/InstanceCallback.hpp"
#include "core/utils/parallel/ParallelObject.hpp"

namespace ScriptInterface {

class ParallelScriptInterfaceSlaveBase {};

class ParallelScriptInterfaceSlave : public Communication::InstanceCallback,
                                     private ParallelScriptInterfaceSlaveBase {
public:
  enum class CallbackAction {
    CREATE,
    SET_PARAMETER,
    SET_PARAMETERS,
    CALL_METHOD,
    DELETE
  };

private:
  friend Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave>;
  ParallelScriptInterfaceSlave() {}

  std::shared_ptr<ScriptInterfaceBase> m_p;

  static std::map<ObjectId, ObjectId> &get_translation_table() {
    static std::map<ObjectId, ObjectId> m_translation_table;

    return m_translation_table;
  }

  /* If the variant encapsulates an object id we translate the
     master id to a local one */
  static void translate_id(Variant &v) {
    try {
      const ObjectId global_id = boost::get<ObjectId>(v);

      v = get_translation_table().at(global_id);
      /* We catch only the bad_get exception, if the id does
         not exsits .at throws out_of_range, which is a real
         error and should be propagated. */
    } catch (boost::bad_get &) {
      ;
    }
  }

private:
  void mpi_slave(int action, int) override {
    switch (CallbackAction(action)) {
    case CallbackAction::CREATE: {
      std::pair<ObjectId, std::string> what;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), what, 0);

      m_p = ScriptInterfaceBase::make_shared(what.second);

      get_translation_table()[what.first] = m_p->id();

      break;
    }
    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

      /* If the parameter is a object we have to tranlate it first to a
         local id.
      */
      translate_id(d.second);
      m_p->set_parameter(d.first, d.second);
      break;
    }
    case CallbackAction::SET_PARAMETERS: {
      std::map<std::string, Variant> parameters;
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), parameters,
                            0);

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
      /* Name of the method and para// meters */
      std::pair<std::string, VariantMap> d;

      /* Broadcast method name and parameters */
      boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

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
};

} /* namespace ScriptInterface */

#endif
