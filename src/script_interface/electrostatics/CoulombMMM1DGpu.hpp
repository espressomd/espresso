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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_MMM1D_GPU_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_MMM1D_GPU_HPP

#include "config.hpp"

#ifdef MMM1D_GPU

#include "Actor.hpp"

#include "core/electrostatics/mmm1d_gpu.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Coulomb {

class CoulombMMM1DGpu : public Actor<CoulombMMM1DGpu, ::CoulombMMM1DGpu> {

public:
  CoulombMMM1DGpu() {
    add_parameters({
        {"is_tuned", AutoParameter::read_only,
         [this]() { return actor()->is_tuned(); }},
        {"far_switch_radius", AutoParameter::read_only,
         [this]() { return actor()->far_switch_radius; }},
        {"maxPWerror", AutoParameter::read_only,
         [this]() { return actor()->maxPWerror; }},
        {"bessel_cutoff", AutoParameter::read_only,
         [this]() { return actor()->bessel_cutoff; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([this, &params]() {
      m_actor = std::make_shared<CoreActorClass>(
          get_value<double>(params, "prefactor"),
          get_value<double>(params, "maxPWerror"),
          get_value<double>(params, "far_switch_radius"),
          get_value<int>(params, "bessel_cutoff"));
    });
    set_charge_neutrality_tolerance(params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // MMM1D_GPU
#endif
