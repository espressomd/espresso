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

#include "ObjectHandle.hpp"
#include "GlobalContext.hpp"

namespace ScriptInterface {
void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  manager()->notify_set_parameter(this, name, value);

  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  manager()->nofity_call_method(this, name, params);

  return this->do_call_method(name, params);
}

void ObjectHandle::delete_remote() {
  if (m_manager)
    manager()->nofity_delete_handle(this);
}

ObjectHandle::~ObjectHandle() { this->do_destroy(); }
} /* namespace ScriptInterface */
