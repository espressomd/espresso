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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_OBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_OBSERVABLE_HPP

#include "script_interface/ScriptInterface.hpp"

#include "core/observables/Observable.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Observables {

/** Base class for script interfaces to core observables classes */
class Observable : public ObjectHandle {
public:
  virtual std::shared_ptr<::Observables::Observable> observable() const = 0;
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "calculate") {
      std::vector<double> out{};
      context()->parallel_try_catch([this, &out]() {
        out = observable()->operator()(context()->get_comm());
      });
      return out;
    }
    if (method == "shape") {
      auto const shape = observable()->shape();
      return std::vector<int>{shape.begin(), shape.end()};
    }
    return {};
  }
};
} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
