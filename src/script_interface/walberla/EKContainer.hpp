/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#ifdef WALBERLA

#include "EKPoissonSolver.hpp"
#include "EKSpecies.hpp"

#include "core/grid_based_algorithms/ek_container.hpp"

#include <script_interface/ObjectList.hpp>
#include <script_interface/ScriptInterface.hpp>

#include <memory>

namespace ScriptInterface::walberla {

class EKContainer : public ObjectList<EKSpecies> {
  void add_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    EK::ek_container.add(obj_ptr->get_ekinstance());
  }
  void remove_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    EK::ek_container.remove(obj_ptr->get_ekinstance());
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "set_tau") {
      EK::ek_container.set_tau(get_value<double>(parameters, "tau"));
      return none;
    }
    if (method == "get_tau") {
      return EK::ek_container.get_tau();
    }
    if (method == "set_poisson_solver") {
      auto obj_ptr =
          get_value<std::shared_ptr<EKPoissonSolver>>(parameters.at("object"));
      EK::ek_container.set_poisson_solver(obj_ptr->get_instance());
      return none;
    }
    if (method == "reset_poisson_solver") {
      EK::ek_container.set_poisson_solver(nullptr);
      return none;
    }
    if (method == "is_poisson_solver_set") {
      return EK::ek_container.is_poisson_solver_set();
    }

    return ObjectList<EKSpecies>::do_call_method(method, parameters);
  }
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
