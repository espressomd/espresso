/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#include "AccumulatorBase.hpp"

#include "core/accumulators/AutoUpdateAccumulators.hpp"
#include "core/system/System.hpp"

#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/system/Leaf.hpp"
#include "script_interface/system/System.hpp"

#include <memory>

namespace ScriptInterface {
namespace Accumulators {

using AutoUpdateAccumulators_t = ObjectList<
    AccumulatorBase,
    AutoParameters<ObjectList<AccumulatorBase, System::Leaf>, System::Leaf>>;

class AutoUpdateAccumulators : public AutoUpdateAccumulators_t {
  using Base = AutoUpdateAccumulators_t;
  std::shared_ptr<::Accumulators::AutoUpdateAccumulators> m_handle;
  std::unique_ptr<VariantMap> m_params;

  bool
  has_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) const override {
    return m_handle->contains(obj_ptr->accumulator().get());
  }

  void add_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) override {
    m_handle->add(obj_ptr->accumulator().get());
  }

  void
  remove_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) override {
    m_handle->remove(obj_ptr->accumulator().get());
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

  void on_bind_system(::System::System &system) override {
    m_handle = system.auto_update_accumulators;
    m_handle->bind_system(m_system.lock());
    Base::do_construct(*m_params);
    m_params.reset();
  }
};
} /* namespace Accumulators */
} /* namespace ScriptInterface */
