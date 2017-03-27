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

#ifndef SCRIPT_INTERFACE_CONSTRAINTS_CONSTRAINT_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_CONSTRAINT_HPP

#include "core/constraints/Constraint.hpp"
#include "core/utils/Factory.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface {
namespace Constraints {

class Constraint : public ScriptInterfaceBase {
public:
  Constraint()
      : m_constraint(new ::Constraints::Constraint()) {}

  const std::string name() const override { return "Constraints::Constraint"; }


  std::shared_ptr<::Constraints::Constraint> constraint() {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::Constraint> m_constraint;

};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
