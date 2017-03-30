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

#ifndef SCRIPT_INTERFACE_CONSTRAINTS_SHAPEBASEDCONSTRAINT_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_SHAPEBASEDCONSTRAINT_HPP

#include "core/constraints/Constraint.hpp"
#include "core/constraints/ShapeBasedConstraint.hpp"
#include "core/utils/Factory.hpp"
#include "script_interface/shapes/Shape.hpp"

namespace ScriptInterface {
namespace Constraints {

class ShapeBasedConstraint : public Constraint {
public:
  ShapeBasedConstraint()
      : m_constraint(new ::Constraints::ShapeBasedConstraint()), m_shape(nullptr) {}

  const std::string name() const override { return "Constraints::ShapeBasedConstraint"; }

  VariantMap get_parameters() const override {
    return {{"only_positive", shape_based_constraint()->only_positive()},
            {"penetrable", shape_based_constraint()->penetrable()},
            {"particle_type", shape_based_constraint()->type()},
            {"shape", (m_shape != nullptr) ? m_shape->id() : ObjectId()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"only_positive", {ParameterType::INT, true}},
            {"penetrable", {ParameterType::INT, true}},
            {"particle_type", {ParameterType::INT, true}},
            {"shape", {ParameterType::OBJECTID, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if ((name == "shape") && (boost::get<ObjectId>(value) != ObjectId())) {
      std::shared_ptr<ScriptInterfaceBase> so_ptr =
          ScriptInterface::get_instance(value);

      auto shape_ptr =
          std::dynamic_pointer_cast<ScriptInterface::Shapes::Shape>(so_ptr);

      /* We are expecting a ScriptInterface::Shapes::Shape here,
         throw if not. That means the assigned object had the wrong type. */
      if (shape_ptr != nullptr) {
        shape_based_constraint()->set_shape(shape_ptr->shape());
        /* Store a reference */
        m_shape = shape_ptr;
      } else {
        throw std::runtime_error("shape parameter expects a Shapes::Shape");
      }
    }

    SET_PARAMETER_HELPER("only_positive", shape_based_constraint()->only_positive());
    SET_PARAMETER_HELPER("penetrable", shape_based_constraint()->penetrable());
    SET_PARAMETER_HELPER("particle_type", shape_based_constraint()->type());
  }

  Variant call_method(std::string const &name, VariantMap const &) override {
    if (name == "total_force") {
      return shape_based_constraint()->total_force();
    }

    return false;
  }


  std::shared_ptr<::Constraints::Constraint> constraint() {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<const ::Constraints::Constraint> constraint() const {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<::Constraints::ShapeBasedConstraint> shape_based_constraint() const {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::ShapeBasedConstraint> m_constraint;

  /* Keep a reference to the shape */
  std::shared_ptr<ScriptInterfaceBase> m_shape;
};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
