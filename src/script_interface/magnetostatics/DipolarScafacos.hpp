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

#ifdef SCAFACOS_DIPOLES

#include "Actor.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/communication.hpp"
#include "core/magnetostatics/scafacos.hpp"
#include "core/scafacos/ScafacosContextBase.hpp"

#include "script_interface/scafacos/scafacos.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Dipoles {

class DipolarScafacos : public Actor<DipolarScafacos, ::DipolarScafacos> {
  std::shared_ptr<boost::mpi::environment> m_mpi_env_lock;

public:
  DipolarScafacos() {
    add_parameters({
        {"method_name", AutoParameter::read_only,
         [this]() { return actor()->get_method(); }},
        {"method_params", AutoParameter::read_only,
         [this]() {
           return Scafacos::deserialize_parameters(actor()->get_parameters());
         }},
    });
  }

  ~DipolarScafacos() override {
    m_actor.reset();
    m_mpi_env_lock.reset();
  }

  void do_construct(VariantMap const &params) override {
    auto const method_name = get_value<std::string>(params, "method_name");
    auto const param_list = params.at("method_params");
    auto const prefactor = get_value<double>(params, "prefactor");

    context()->parallel_try_catch([&]() {
      if (prefactor <= 0.) {
        throw std::domain_error("Parameter 'prefactor' must be > 0");
      }
      ScafacosContextBase::sanity_check_method(method_name);
      auto const method_params = Scafacos::serialize_parameters(param_list);
      m_actor = make_dipolar_scafacos(method_name, method_params);
      actor()->set_prefactor(prefactor);
    });
    // MPI communicator is needed to destroy the FFT plans
    m_mpi_env_lock = ::Communication::mpiCallbacksHandle()->share_mpi_env();
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_available_methods") {
      return make_vector_of_variants(Scafacos::available_methods());
    }
    return Actor<SIActorClass, CoreActorClass>::do_call_method(name, params);
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // SCAFACOS_DIPOLES
