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

#include "ScriptInterfaceBase.hpp"
#include "ParallelScriptInterface.hpp"
#include "Serializer.hpp"
#include "utils/Factory.hpp"

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <sstream>

namespace iostreams = boost::iostreams;

namespace ScriptInterface {
std::shared_ptr<ScriptInterfaceBase>
ScriptInterfaceBase::make_shared(std::string const &name,
                                 CreationPolicy policy) {
  std::shared_ptr<ScriptInterfaceBase> sp;

  switch (policy) {
  case CreationPolicy::LOCAL:
    sp = Utils::Factory<ScriptInterfaceBase>::make(name);
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
 * @brief  Return a Variant representation of the state of the object.
 *
 * This should returne the internal state of the instance, so that
 * the instance can be restored from this information.  The default
 * implementation stores all the public parameters, including object
 * parameters that are captured by calling get_state on them.
 */
Variant ScriptInterfaceBase::get_state() const {
  std::vector<Variant> state;

  auto params = this->get_parameters();
  state.reserve(params.size());

  for (auto const &p : params) {
    state.push_back(std::vector<Variant>{
        {p.first, boost::apply_visitor(Serializer{}, p.second)}});
  }

  return state;
}

void ScriptInterfaceBase::set_state(Variant const &state) {
  using boost::get;
  using std::vector;

  VariantMap params;
  UnSerializer u;

  for (auto const &v : get<vector<Variant>>(state)) {
    auto const &p = get<vector<Variant>>(v);
    params[get<std::string>(p.at(0))] = boost::apply_visitor(u, p.at(1));
  }

  this->construct(params);
}

/**
 * @brief Returns a binary representation of the state often
 *        the instance, as returned by @f get_state().
 */
std::string ScriptInterfaceBase::serialize() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);
  auto v = Serializer{}(this->id());

  oa << v;
  return ss.str();
}

/**
 * @brief Creates a new instance from a binary state,
 *        as returned by @f serialize().
 */
std::shared_ptr<ScriptInterfaceBase>
ScriptInterfaceBase::unserialize(std::string const &state) {
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  Variant v;
  ia >> v;

  UnSerializer u;
  auto oid = boost::get<ObjectId>(boost::apply_visitor(u, v));

  return get_instance(oid).lock();
}

} /* namespace ScriptInterface */
