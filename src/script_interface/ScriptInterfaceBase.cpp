/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "ScriptInterfaceBase.hpp"
#include "ParallelScriptInterface.hpp"
#include "ScriptInterface.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <sstream>

namespace ScriptInterface {
Utils::Factory<ScriptInterfaceBase> factory;

std::shared_ptr<ScriptInterfaceBase>
ScriptInterfaceBase::make_shared(std::string const &name,
                                 CreationPolicy policy) {
  std::shared_ptr<ScriptInterfaceBase> sp;

  switch (policy) {
  case CreationPolicy::LOCAL:
    sp = factory.make(name);
    break;
  case CreationPolicy::GLOBAL:
    sp =
        std::shared_ptr<ScriptInterfaceBase>(new ParallelScriptInterface(name));
    break;
  }

  /* Set the policy and the name */
  sp->set_policy(policy);
  sp->set_name(name);

  /* Id of the newly created instance */
  const auto id = sp->id();

  /* Now get a reference to the corresponding weak_ptr in ObjectId and update
     it with our shared ptr, so that everybody uses the same ref count.
  */
  sp->get_instance(id) = sp;

  return sp;
}

std::weak_ptr<ScriptInterfaceBase> &
ScriptInterfaceBase::get_instance(ObjectId id) {
  return Utils::AutoObjectId<ScriptInterfaceBase>::get_instance(id);
}

/* Checkpointing functions. */

/**
 * @brief Returns a binary representation of the state often
 *        the instance, as returned by get_state().
 */
std::string ScriptInterfaceBase::serialize() const {
  return {};
}

/**
 * @brief Creates a new instance from a binary state,
 *        as returned by serialize().
 */
std::shared_ptr<ScriptInterfaceBase>
ScriptInterfaceBase::unserialize(std::string const &state) {
  return {};
}

} /* namespace ScriptInterface */
