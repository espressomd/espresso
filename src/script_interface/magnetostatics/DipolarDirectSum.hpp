/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ESPRESSO_DIPOLES

#include "Actor.hpp"

#include "core/magnetostatics/dipolar_direct_sum.hpp"

#include "script_interface/code_info/CodeInfo.hpp"
#include "script_interface/get_value.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Dipoles {

class DipolarDirectSum : public Actor<DipolarDirectSum, ::DipolarDirectSum> {
public:
  DipolarDirectSum() {
    add_parameters({
        {"n_replicas", AutoParameter::read_only,
         [this]() { return actor()->n_replicas; }},
        {"gpu", AutoParameter::read_only,
         [this]() { return actor()->is_gpu(); }},
    });
  }

  void do_construct(VariantMap const &params) override {
    auto const gpu = get_value_or(params, "gpu", false);
    context()->parallel_try_catch([this, &params, gpu]() {
      if (gpu) {
        std::vector<std::string> required_features;
        required_features.emplace_back("CUDA");
        CodeInfo::check_features(required_features);
      }
      m_actor = std::make_shared<CoreActorClass>(
          get_value<double>(params, "prefactor"),
          get_value<int>(params, "n_replicas"), gpu);
    });
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // ESPRESSO_DIPOLES
