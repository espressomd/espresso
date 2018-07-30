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

#include "ScriptInterfaceBase.hpp"
#include "core/utils/parallel/ParallelObject.hpp"

namespace ScriptInterface {

class ParallelScriptInterfaceSlaveBase {};

class ParallelScriptInterfaceSlave : private ParallelScriptInterfaceSlaveBase {
public:
  enum class CallbackAction {
    NEW,
    CONSTRUCT,
    SET_PARAMETER,
    SET_PARAMETERS,
    CALL_METHOD,
    DELETE
  };

private:
  friend class ParallelScriptInterface;
  static Communication::MpiCallbacks *m_cb;
  friend Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave>;
  ParallelScriptInterfaceSlave();

  std::shared_ptr<ScriptInterfaceBase> m_p;

  static std::map<ObjectId, ObjectId> &get_translation_table();

  /* If the variant encapsulates an object id we translate the
     master id to a local one */
  static void translate_id(Variant &v) {
    if (is_objectid(v)) {
      v = get_translation_table().at(boost::get<ObjectId>(v));
    }
  }

  VariantMap bcast_variant_map() const;

private:
  void mpi_slave(int action, int);
};
} /* namespace ScriptInterface */

#endif
