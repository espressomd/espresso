/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/system/System.hpp"

#include "core/accumulators/AccumulatorBase.hpp"
#include "core/system/System.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Accumulators {

class AccumulatorBase : public AutoParameters<AccumulatorBase> {
public:
  AccumulatorBase() {
    add_parameters({{"delta_N",
                     [this](const Variant &v) {
                       accumulator()->delta_N() = get_value<int>(v);
                     },
                     [this]() { return accumulator()->delta_N(); }}});
  }
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "shape") {
      auto const shape = accumulator()->shape();
      return std::vector<int>{shape.begin(), shape.end()};
    }
    if (method == "reload_from_checkpoint") {
      accumulator()->set_internal_state(
          get_value<std::string>(parameters, "state"));
      return {};
    }
    return {};
  }
  virtual std::shared_ptr<const ::Accumulators::AccumulatorBase>
  accumulator() const = 0;
  virtual std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() = 0;

protected:
  auto get_core_system_pointer(VariantMap const &params) const {
    auto const *system = &::System::get_system();
    if (params.contains("system")) {
      system =
          &(get_value<std::shared_ptr<System::System const>>(params, "system")
                ->get_system());
    }
    return system;
  }

private:
  std::string get_internal_state() const override {
    return accumulator()->get_internal_state();
  }

  void set_internal_state(std::string const &state) override {
    call_method("reload_from_checkpoint", {{"state", state}});
  }
};

} // namespace Accumulators
} // namespace ScriptInterface
