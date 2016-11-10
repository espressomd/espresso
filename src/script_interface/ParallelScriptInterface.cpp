/*
  Copyright (C) 2015,2016 The ESPResSo project

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

#include "ParallelScriptInterface.hpp"
#include "utils/parallel/ParallelObject.hpp"

namespace ScriptInterface {

bool ParallelScriptInterface::operator==(ParallelScriptInterface const &rhs) {
  return this->get_underlying_object() == rhs.get_underlying_object();
}

bool ParallelScriptInterface::operator!=(ParallelScriptInterface const &rhs) {
  return !(*this == rhs);
}

ParallelScriptInterface::ParallelScriptInterface(std::string const &name) {
  /* Create the slaves */
  Utils::Parallel::ParallelObject<ParallelScriptInterfaceSlave>::make();

  /* Create local object */
  m_p = ScriptInterfaceBase::make_shared(
      name, ScriptInterfaceBase::CreationPolicy::LOCAL);

  /* Bcast class name and global id to the slaves */
  call(static_cast<int>(CallbackAction::CREATE));
  std::pair<ObjectId, std::string> what = std::make_pair(m_p->id(), name);
  boost::mpi::broadcast(Communication::mpiCallbacks().comm(), what, 0);
}

void ParallelScriptInterface::initialize() {
  Utils::Parallel::ParallelObject<
      ParallelScriptInterfaceSlave>::register_callback();
}

void ParallelScriptInterface::set_parameter(const std::string &name,
                                            const Variant &value) {
  std::pair<std::string, Variant> d(name, Variant());

  if (is_objectid(value)) {
    d.second = map_parallel_to_local_id(name, value);
  } else {
    d.second = value;
  }

  call(static_cast<int>(CallbackAction::SET_PARAMETER), 0);

  boost::mpi::broadcast(Communication::mpiCallbacks().comm(), d, 0);

  m_p->set_parameter(d.first, d.second);
}

void ParallelScriptInterface::set_parameters(const VariantMap &parameters) {
  call(static_cast<int>(CallbackAction::SET_PARAMETERS), 0);

  /* Copy parameters into a non-const buffer, needed by boost::mpi */
  auto p = parameters;

  /* Unwrap the object ids */
  for (auto &it : p) {
    if (is_objectid(it.second)) {
      it.second = map_parallel_to_local_id(it.first, it.second);
    }
  }

  boost::mpi::broadcast(Communication::mpiCallbacks().comm(), p, 0);

  m_p->set_parameters(p);
}

} /* namespace ScriptInterface */
