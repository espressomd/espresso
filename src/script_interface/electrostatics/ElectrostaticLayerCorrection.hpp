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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_ELC_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_ELC_HPP

#include "config/config.hpp"

#ifdef P3M

#include "Actor.hpp"

#include "CoulombP3M.hpp"

#include "core/electrostatics/elc.hpp"

#include "script_interface/get_value.hpp"

#include "boost/variant.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Coulomb {

class ElectrostaticLayerCorrection
    : public Actor<ElectrostaticLayerCorrection,
                   ::ElectrostaticLayerCorrection> {

  using BaseSolver = boost::variant<
#ifdef CUDA
      std::shared_ptr<CoulombP3MGPU>,
#endif // CUDA
      std::shared_ptr<CoulombP3M>>;
  BaseSolver m_solver;

public:
  ElectrostaticLayerCorrection() {
    add_parameters({
        {"maxPWerror", AutoParameter::read_only,
         [this]() { return actor()->elc.maxPWerror; }},
        {"gap_size", AutoParameter::read_only,
         [this]() { return actor()->elc.gap_size; }},
        {"far_cut", AutoParameter::read_only,
         [this]() { return actor()->elc.far_cut; }},
        {"neutralize", AutoParameter::read_only,
         [this]() { return actor()->elc.neutralize; }},
        {"delta_mid_top", AutoParameter::read_only,
         [this]() { return actor()->elc.delta_mid_top; }},
        {"delta_mid_bot", AutoParameter::read_only,
         [this]() { return actor()->elc.delta_mid_bot; }},
        {"const_pot", AutoParameter::read_only,
         [this]() { return actor()->elc.const_pot; }},
        {"pot_diff", AutoParameter::read_only,
         [this]() { return actor()->elc.pot_diff; }},
        {"actor", AutoParameter::read_only,
         [this]() {
           return boost::apply_visitor(
               [](auto &solver) { return Variant{solver}; }, m_solver);
         }},
    });
  }

  void do_construct(VariantMap const &params) override {
    ::ElectrostaticLayerCorrection::BaseSolver solver;
    auto so_ptr = get_value<ObjectRef>(params, "actor");
    context()->parallel_try_catch([&]() {
#ifdef CUDA
      if (auto so_solver = std::dynamic_pointer_cast<CoulombP3MGPU>(so_ptr)) {
        solver = so_solver->actor();
        m_solver = so_solver;
      } else
#endif // CUDA
        if (auto so_solver = std::dynamic_pointer_cast<CoulombP3M>(so_ptr)) {
          solver = so_solver->actor();
          m_solver = so_solver;
        } else {
          throw std::invalid_argument("Parameter 'actor' of type " +
                                      so_ptr->name().to_string() +
                                      " isn't supported by ELC");
        }
    });
    context()->parallel_try_catch([&]() {
      auto elc = elc_data{get_value<double>(params, "maxPWerror"),
                          get_value<double>(params, "gap_size"),
                          get_value<double>(params, "far_cut"),
                          get_value<bool>(params, "neutralize"),
                          get_value<double>(params, "delta_mid_top"),
                          get_value<double>(params, "delta_mid_bot"),
                          get_value<bool>(params, "const_pot"),
                          get_value<double>(params, "pot_diff")};
      m_actor =
          std::make_shared<CoreActorClass>(std::move(elc), std::move(solver));
    });
    set_charge_neutrality_tolerance(params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // P3M
#endif
