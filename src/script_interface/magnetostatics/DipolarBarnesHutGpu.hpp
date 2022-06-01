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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_MAGNETOSTATICS_DIPOLAR_BH_GPU_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_MAGNETOSTATICS_DIPOLAR_BH_GPU_HPP

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "Actor.hpp"

#include "core/magnetostatics/barnes_hut_gpu.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Dipoles {

class DipolarBarnesHutGpu
    : public Actor<DipolarBarnesHutGpu, ::DipolarBarnesHutGpu> {
public:
  DipolarBarnesHutGpu() {
    add_parameters({
        {"epssq", AutoParameter::read_only,
         [this]() { return actor()->m_epssq; }},
        {"itolsq", AutoParameter::read_only,
         [this]() { return actor()->m_itolsq; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_actor =
        std::make_shared<CoreActorClass>(get_value<double>(params, "prefactor"),
                                         get_value<double>(params, "epssq"),
                                         get_value<double>(params, "itolsq"));
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DIPOLAR_BARNES_HUT
#endif
