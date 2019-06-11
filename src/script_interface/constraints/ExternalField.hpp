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

#ifndef SCRIPT_INTERFACE_CONSTRAINTS_EXTERNAL_FIELD_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_EXTERNAL_FIELD_HPP

#include "core/constraints/ExternalField.hpp"
#include "script_interface/ScriptInterface.hpp"

#include "couplings.hpp"
#include "fields.hpp"

#include "core/grid.hpp"

namespace ScriptInterface {
namespace Constraints {

template <typename Coupling, typename Field>
class ExternalField : public Constraint {
  using CoreField = ::Constraints::ExternalField<Coupling, Field>;

public:
  ExternalField() {
    add_parameters(detail::coupling_parameters<Coupling>(
        [this]() { return m_constraint->coupling(); }));
    add_parameters(detail::field_parameters<Field>(
        [this]() { return m_constraint->field(); }));
  }

  void construct(VariantMap const &args) override {
    m_constraint = std::make_shared<CoreField>(
        detail::make_coupling<Coupling>(args), detail::make_field<Field>(args));
  }

  Variant call_method(const std::string &name,
                      VariantMap const &args) override {
    if (name == "_eval_field") {
      return m_constraint->field()(get_value<Utils::Vector3d>(args, "x"),
                                   get_value_or<double>(args, "t", 0.));
    }
    return none;
  }

  std::shared_ptr<::Constraints::Constraint> constraint() override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }

  std::shared_ptr<const ::Constraints::Constraint> constraint() const override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }

private:
  /* The actual constraint */
  std::shared_ptr<CoreField> m_constraint;
};
} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
