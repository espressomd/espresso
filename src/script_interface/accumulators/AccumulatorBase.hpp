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

#include "core/accumulators/AccumulatorBase.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

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
  Variant call_method(std::string const &method, VariantMap const &parameters) {
    if (method == "shape") {
      auto const shape = accumulator()->shape();
      return std::vector<int>{shape.begin(), shape.end()};
    }
    return {};
  }
  virtual std::shared_ptr<const ::Accumulators::AccumulatorBase>
  accumulator() const = 0;
  virtual std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() = 0;
};

} // namespace Accumulators
} // namespace ScriptInterface
