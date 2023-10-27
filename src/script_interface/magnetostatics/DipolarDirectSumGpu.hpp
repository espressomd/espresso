/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifdef DIPOLAR_DIRECT_SUM

#include "Actor.hpp"

#include "core/magnetostatics/dipolar_direct_sum_gpu.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Dipoles {

class DipolarDirectSumGpu
    : public Actor<DipolarDirectSumGpu, ::DipolarDirectSumGpu> {
public:
  DipolarDirectSumGpu() = default;

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([this, &params]() {
      m_actor = std::make_shared<CoreActorClass>(
          get_value<double>(params, "prefactor"));
    });
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DIPOLAR_DIRECT_SUM
