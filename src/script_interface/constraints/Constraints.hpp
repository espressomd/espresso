/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "Constraint.hpp"

#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/system/Leaf.hpp"
#include "script_interface/system/System.hpp"

#include "core/constraints/Constraints.hpp"
#include "core/system/System.hpp"

#include <memory>

namespace ScriptInterface {
namespace Constraints {

using Constraints_t =
    ObjectList<Constraint, AutoParameters<ObjectList<Constraint, System::Leaf>,
                                          System::Leaf>>;

class Constraints : public Constraints_t {
  using Base = Constraints_t;
  std::shared_ptr<::Constraints::Constraints> m_handle;
  std::unique_ptr<VariantMap> m_params;

  bool has_in_core(std::shared_ptr<Constraint> const &obj_ptr) const override {
    return m_handle->contains(obj_ptr->constraint());
  }
  void add_in_core(std::shared_ptr<Constraint> const &obj_ptr) override {
    m_handle->add(obj_ptr->constraint());
    obj_ptr->bind_system(m_system.lock());
  }
  void remove_in_core(std::shared_ptr<Constraint> const &obj_ptr) override {
    m_handle->remove(obj_ptr->constraint());
  }

  void do_construct(VariantMap const &params) override {
    m_handle = std::make_shared<::Constraints::Constraints>();
    m_handle->bind_system(::System::get_system().shared_from_this());
    m_params = std::make_unique<VariantMap>(params);
  }

  void on_bind_system(::System::System &system) override {
    m_handle = system.constraints;
    m_handle->bind_system(m_system.lock());
    Base::do_construct(*m_params);
    m_params.reset();
  }
};

} // namespace Constraints
} // namespace ScriptInterface
