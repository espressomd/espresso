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
#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_DIPOLAR_LAYER_CORRECTION_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_DIPOLAR_LAYER_CORRECTION_HPP

#include "config/config.hpp"

#ifdef DIPOLES

#include "Actor.hpp"

#include "DipolarDirectSum.hpp"
#include "DipolarP3M.hpp"

#include "core/magnetostatics/dlc.hpp"

#include "script_interface/get_value.hpp"

#include "boost/variant.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Dipoles {

class DipolarLayerCorrection
    : public Actor<DipolarLayerCorrection, ::DipolarLayerCorrection> {
  using DipolarDSR = DipolarDirectSum;
  using BaseSolver = boost::variant<
#ifdef DP3M
      std::shared_ptr<DipolarP3M>,
#endif
      std::shared_ptr<DipolarDSR>>;
  BaseSolver m_solver;

public:
  DipolarLayerCorrection() {
    add_parameters({
        {"maxPWerror", AutoParameter::read_only,
         [this]() { return actor()->dlc.maxPWerror; }},
        {"gap_size", AutoParameter::read_only,
         [this]() { return actor()->dlc.gap_size; }},
        {"far_cut", AutoParameter::read_only,
         [this]() { return actor()->dlc.far_cut; }},
        {"actor", AutoParameter::read_only,
         [this]() {
           return boost::apply_visitor(
               [](auto &solver) { return Variant{solver}; }, m_solver);
         }},
    });
  }

  void do_construct(VariantMap const &params) override {
    ::DipolarLayerCorrection::BaseSolver solver;
    auto so_ptr = get_value<ObjectRef>(params, "actor");
    context()->parallel_try_catch([&]() {
#ifdef DP3M
      if (auto so_solver = std::dynamic_pointer_cast<DipolarP3M>(so_ptr)) {
        solver = so_solver->actor();
        m_solver = so_solver;
      } else
#endif // DP3M
        if (auto so_solver = std::dynamic_pointer_cast<DipolarDSR>(so_ptr)) {
          solver = so_solver->actor();
          m_solver = so_solver;
        } else {
          throw std::invalid_argument("Parameter 'actor' of type " +
                                      so_ptr->name().to_string() +
                                      " isn't supported by DLC");
        }
    });
    context()->parallel_try_catch([&]() {
      auto dlc = dlc_data(get_value<double>(params, "maxPWerror"),
                          get_value<double>(params, "gap_size"),
                          get_value<double>(params, "far_cut"));
      m_actor =
          std::make_shared<CoreActorClass>(std::move(dlc), std::move(solver));
    });
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DIPOLES
#endif
