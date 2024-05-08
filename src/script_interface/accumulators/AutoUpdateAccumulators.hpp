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

#include "core/accumulators.hpp"
#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface {
namespace Accumulators {
class AutoUpdateAccumulators : public ObjectList<AccumulatorBase> {
  bool
  has_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) const override {
    return ::Accumulators::auto_update_contains(obj_ptr->accumulator().get());
  }

  void add_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) override {
    ::Accumulators::auto_update_add(obj_ptr->accumulator().get());
  }

  void
  remove_in_core(std::shared_ptr<AccumulatorBase> const &obj_ptr) override {
    ::Accumulators::auto_update_remove(obj_ptr->accumulator().get());
  }

private:
  // disable serialization: pickling done by the python interface
  std::string get_internal_state() const override { return {}; }
  void set_internal_state(std::string const &) override {}
};
} /* namespace Accumulators */
} /* namespace ScriptInterface */
