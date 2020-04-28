/*
 * Copyright (C) 2015-2019 The ESPResSo project
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

#include "script_interface/initialize.hpp"
#include "script_interface/ScriptInterface.hpp"

#include "AutoUpdateAccumulators.hpp"
#include "Correlator.hpp"
#include "MeanVarianceCalculator.hpp"
#include "TimeSeries.hpp"

namespace ScriptInterface {
namespace Accumulators {

void initialize() {
  ScriptInterface::register_new<
      ScriptInterface::Accumulators::AutoUpdateAccumulators>(
      "Accumulators::AutoUpdateAccumulators");

  ScriptInterface::register_new<
      ScriptInterface::Accumulators::MeanVarianceCalculator>(
      "Accumulators::MeanVarianceCalculator");

  ScriptInterface::register_new<ScriptInterface::Accumulators::TimeSeries>(
      "Accumulators::TimeSeries");

  ScriptInterface::register_new<ScriptInterface::Accumulators::Correlator>(
      "Accumulators::Correlator");
}
} /* namespace Accumulators */
} /* namespace ScriptInterface */
