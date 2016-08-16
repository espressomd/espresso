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

#include "ScriptInterface.hpp"
#include "core/constraints/Constraint.hpp"
#include "core/utils/Factory.hpp"
#include "shapes/Shape.hpp"

#include <memory>

namespace ScriptInterface {
namespace Constraints {

class Constraint : public ScriptInterfaceBase {
public:
  Constraint() : m_constraint(new ::Constraints::Constraint()) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }
  const std::string name() const override { return "Constraints::Constraint"; }

  VariantMap get_parameters() const override {
    return {{"only_positive", m_constraint->only_positive()},
            {"penetrable", m_constraint->penetrable()},
            {"type", m_constraint->type()},
            {"shape",
             (m_shape != nullptr) ? m_shape->id() : ScriptInterface::NOT_SET}};
  }

  ParameterMap all_parameters() const override {
    return {{"only_positive", {ParameterType::INT, true}},
            {"penetrable", {ParameterType::INT, true}},
            {"type", {ParameterType::INT, true}},
            {"shape", {ParameterType::OBJECT, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if (name == "shape") {
      std::shared_ptr<ScriptInterfaceBase> so_ptr =
          ScriptInterface::get_instance(value);
      std::cout << __PRETTY_FUNCTION__ << " value = " << boost::get<int>(value)
                << std::endl;

      std::cout << __PRETTY_FUNCTION__ << " so_ptr = " << so_ptr.get()
                << std::endl;

      auto shape_ptr =
          std::dynamic_pointer_cast<ScriptInterface::Shapes::Shape>(so_ptr);

      /* We are expecting a ScriptInterface::Shapes::Shape here,
         throw if not. That means the assigned object had the wrong type. */
      if (shape_ptr != nullptr) {
        m_constraint->set_shape(shape_ptr->shape());
        /* Store a reference */
        m_shape = shape_ptr;
      } else {
        throw std::runtime_error("shape parameter expects a Shapes::Shape");
      }
    }

    SET_PARAMETER_HELPER("only_positive", m_constraint->only_positive());
    SET_PARAMETER_HELPER("penetrable", m_constraint->penetrable());
    SET_PARAMETER_HELPER("type", m_constraint->type());
  }

  std::shared_ptr<::Constraints::Constraint> constraint() {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::Constraint> m_constraint;

  /* Keep a reference to the shape */
  std::shared_ptr<ScriptInterfaceBase> m_shape;
};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
