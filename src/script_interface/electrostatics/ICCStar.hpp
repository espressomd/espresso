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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_ICC_STAR_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_ICC_STAR_HPP

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "core/actor/registration.hpp"
#include "core/electrostatics/icc.hpp"

#include <utils/Vector.hpp>

#include "script_interface/Context.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace Coulomb {

class ICCStar : public AutoParameters<ICCStar> {
  using CoreActorClass = ::ICCStar;
  std::shared_ptr<CoreActorClass> m_actor;

public:
  ICCStar() {
    add_parameters({
        {"n_icc", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.n_icc; }},
        {"max_iterations", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.max_iterations; }},
        {"eps_out", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.eps_out; }},
        {"areas", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.areas; }},
        {"epsilons", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.epsilons; }},
        {"sigmas", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.sigmas; }},
        {"convergence", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.convergence; }},
        {"normals", AutoParameter::read_only,
         [this]() {
           return make_vector_of_variants(actor()->icc_cfg.normals);
         }},
        {"ext_field", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.ext_field; }},
        {"relaxation", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.relaxation; }},
        {"citeration", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.citeration; }},
        {"first_id", AutoParameter::read_only,
         [this]() { return actor()->icc_cfg.first_id; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    auto icc_parameters = ::icc_data{
        get_value<int>(params, "n_icc"),
        get_value<int>(params, "max_iterations"),
        get_value<double>(params, "eps_out"),
        get_value<std::vector<double>>(params, "areas"),
        get_value<std::vector<double>>(params, "epsilons"),
        get_value<std::vector<double>>(params, "sigmas"),
        get_value<double>(params, "convergence"),
        get_value<std::vector<Utils::Vector3d>>(params, "normals"),
        get_value<Utils::Vector3d>(params, "ext_field"),
        get_value<double>(params, "relaxation"),
        0,
        get_value<int>(params, "first_id"),
    };
    context()->parallel_try_catch([&]() {
      m_actor = std::make_shared<CoreActorClass>(std::move(icc_parameters));
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "activate") {
      context()->parallel_try_catch([&]() {
        auto &system = System::get_system();
        add_actor(context()->get_comm(), system.coulomb.impl->extension,
                  m_actor, [&system]() { system.on_coulomb_change(); });
      });
      return {};
    }
    return {};
  }

  std::shared_ptr<CoreActorClass> actor() { return m_actor; }
  std::shared_ptr<CoreActorClass const> actor() const { return m_actor; }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ELECTROSTATICS
#endif
