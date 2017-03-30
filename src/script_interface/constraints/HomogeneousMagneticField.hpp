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

#ifndef SCRIPT_INTERFACE_CONSTRAINTS_HOMOGENEOUSMAGNETICFIELD_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_HOMOGENEOUSMAGNETICFIELD_HPP

#include "core/constraints/Constraint.hpp"
#include "core/constraints/HomogeneousMagneticField.hpp"

namespace ScriptInterface {
namespace Constraints {

class HomogeneousMagneticField : public Constraint {
public:
  HomogeneousMagneticField()
    : m_constraint(new ::Constraints::HomogeneousMagneticField())
  {}

  const std::string name() const override { return "Constraints::HomogeneousMagneticField"; }

  VariantMap get_parameters() const override {
    return {{"H", homogeneous_magnetic_field()->H()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"H", {ParameterType::DOUBLE_VECTOR, 3, true}}
    };
  }

  void set_parameter(std::string const &name, ScriptInterface::Variant const &value) override {
    if (name == "H") {
      homogeneous_magnetic_field()->set_H(get_value<Vector3d>(value));
    }
  }

  std::shared_ptr<::Constraints::Constraint> constraint() {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<const ::Constraints::Constraint> constraint() const {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<::Constraints::HomogeneousMagneticField> homogeneous_magnetic_field() const {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::HomogeneousMagneticField> m_constraint;

};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
