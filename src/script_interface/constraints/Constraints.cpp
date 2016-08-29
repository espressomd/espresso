/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include <stdexcept>

#include "Constraint.hpp"
#include "Constraints.hpp"

#include "core/constraints.hpp"

#include <iostream>

namespace ScriptInterface {
namespace Constraints {

Variant Constraints::call_method(std::string const &method,
                                 VariantMap const &parameters) {
  Variant par = parameters.at("constraint");

  auto so_ptr = ScriptInterface::get_instance(par);

  auto constraint_ptr =
      std::dynamic_pointer_cast<ScriptInterface::Constraints::Constraint>(
          so_ptr);

  if (constraint_ptr == nullptr)
    throw std::runtime_error("Wrong type");

  if (method == "add")
    ::Constraints::constraints.add(constraint_ptr->constraint());
  if (method == "remove")
    ::Constraints::constraints.remove(constraint_ptr->constraint());

  return {};
}
}
}
