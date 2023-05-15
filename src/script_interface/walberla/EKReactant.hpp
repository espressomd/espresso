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

#include "EKSpecies.hpp"

#include <walberla_bridge/electrokinetics/reactions/EKReactant.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <memory>

namespace ScriptInterface::walberla {
class EKReactant : public AutoParameters<::walberla::EKReactant> {
public:
  void do_construct(VariantMap const &args) override {
    m_ekreactant = std::make_shared<::walberla::EKReactant>(
        get_value<std::shared_ptr<EKSpecies>>(args, "ekspecies")
            ->get_ekinstance(),
        get_value<double>(args, "stoech_coeff"),
        get_value<double>(args, "order"));

    add_parameters({{"order", AutoParameter::read_only,
                     [this]() { return m_ekreactant->get_order(); }},
                    {"stoech_coeff",
                     [this](Variant const &v) {
                       m_ekreactant->set_stoech_coefficient(
                           get_value<double>(v));
                     },
                     [this]() { return m_ekreactant->get_stoech_coeff(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactant> get_instance() {
    return m_ekreactant;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::EKReactant> m_ekreactant;
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
