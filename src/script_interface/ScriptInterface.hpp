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

#ifndef SCRIPT_INTERFACE_SCRIPT_INTERFACE_HPP
#define SCRIPT_INTERFACE_SCRIPT_INTERFACE_HPP

#include <type_traits>

#include "Variant.hpp"

#include "ScriptInterfaceBase.hpp"
#include "auto_parameters/AutoParameters.hpp"
#include "get_value.hpp"
#include "initialize.hpp"
#include "utils/Factory.hpp"

namespace ScriptInterface {
template <typename T> static void register_new(std::string const &name) {
  static_assert(std::is_base_of<ScriptInterfaceBase, T>::value, "");

  /* Register with the factory */
  Utils::Factory<ScriptInterfaceBase>::register_new<T>(name);
}

template <typename T> static void register_new() {
  register_new<T>(T{}.name());
}

inline std::shared_ptr<ScriptInterfaceBase> get_instance(Variant value) {
  const auto id = boost::get<ObjectId>(value);

  return ScriptInterfaceBase::get_instance(id).lock();
}
} /* namespace ScriptInterface */

#endif
