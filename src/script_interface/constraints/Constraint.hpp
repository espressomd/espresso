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
#include "script_interface/shapes/Shape.hpp"

namespace ScriptInterface {
namespace Constraints {

class Constraint : public ScriptInterfaceBase {
public:
  Constraint()
      : m_constraint(new ::Constraints::Constraint()), m_shape(nullptr) {}

  const std::string name() const override { return "Constraints::Constraint"; }

  VariantMap get_parameters() const override {
    return {{"only_positive", m_constraint->only_positive()},
            {"penetrable", m_constraint->penetrable()},
            {"particle_type", m_constraint->type()},
            {"shape", (m_shape != nullptr) ? m_shape->id() : ObjectId()},
            {"ext_electric_field", m_constraint->ext_electric_field()},
            {"ext_magn_field", m_constraint->ext_magn_field()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"only_positive", {ParameterType::INT, true}},
            {"penetrable", {ParameterType::INT, true}},
            {"particle_type", {ParameterType::INT, true}},
            {"shape", {ParameterType::OBJECTID, true}},
            {"ext_electric_field", {ParameterType::DOUBLE, true}},
            {"ext_magn_field", {ParameterType::DOUBLE, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if (name == "shape") {
      m_shape = get_value<std::shared_ptr<Shapes::Shape>>(value);

      if (m_shape) {
        m_constraint->set_shape(m_shape->shape());
      }
    }

    SET_PARAMETER_HELPER("only_positive", m_constraint->only_positive());
    SET_PARAMETER_HELPER("penetrable", m_constraint->penetrable());
    SET_PARAMETER_HELPER("particle_type", m_constraint->type());
    SET_PARAMETER_HELPER("ext_electric_field", m_constraint->ext_electric_field());
    SET_PARAMETER_HELPER("ext_magn_field", m_constraint->ext_magn_field());
  }

  Variant call_method(std::string const &name, VariantMap const &) override {
    if (name == "total_force") {
      return m_constraint->total_force();
    }

    return false;
  }

  std::shared_ptr<::Constraints::Constraint> constraint() {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::Constraint> m_constraint;

  /* Keep a reference to the shape */
  std::shared_ptr<Shapes::Shape> m_shape;
};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
