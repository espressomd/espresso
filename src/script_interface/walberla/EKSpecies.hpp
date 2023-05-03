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

#include "LatticeWalberla.hpp"

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <utils/math/int_pow.hpp>

#include <memory>
#include <string>

namespace ScriptInterface::walberla {

class EKSpecies : public AutoParameters<EKSpecies> {
private:
  double m_conv_diffusion;
  double m_conv_ext_efield;
  double m_conv_density;
  double m_conv_flux;
  double m_tau;
  double m_density;

public:
  EKSpecies() {
    add_parameters(
        {{"diffusion",
          [this](Variant const &v) {
            m_ekinstance->set_diffusion(get_value<double>(v) *
                                        m_conv_diffusion);
          },
          [this]() {
            return m_ekinstance->get_diffusion() / m_conv_diffusion;
          }},
         {"kT",
          [this](Variant const &v) {
            m_ekinstance->set_kT(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_kT(); }},
         {"valency",
          [this](Variant const &v) {
            m_ekinstance->set_valency(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_valency(); }},
         {"ext_efield",
          [this](Variant const &v) {
            m_ekinstance->set_ext_efield(get_value<Utils::Vector3d>(v) *
                                         m_conv_ext_efield);
          },
          [this]() {
            return m_ekinstance->get_ext_efield() / m_conv_ext_efield;
          }},
         {"advection",
          [this](Variant const &v) {
            m_ekinstance->set_advection(get_value<bool>(v));
          },
          [this]() { return m_ekinstance->get_advection(); }},
         {"friction_coupling",
          [this](Variant const &v) {
            m_ekinstance->set_friction_coupling(get_value<bool>(v));
          },
          [this]() { return m_ekinstance->get_friction_coupling(); }},
         {"is_single_precision", AutoParameter::read_only,
          [this]() { return not m_ekinstance->is_double_precision(); }},
         {"tau", AutoParameter::read_only, [this]() { return m_tau; }},
         {"density", AutoParameter::read_only,
          [this]() { return m_density / m_conv_density; }},
         {"shape", AutoParameter::read_only,
          [this]() {
            return m_ekinstance->get_lattice().get_grid_dimensions();
          }},
         {"lattice", AutoParameter::read_only,
          [this]() { return m_lattice; }}});
  }

  void do_construct(VariantMap const &args) override;

  [[nodiscard]] std::shared_ptr<::EKinWalberlaBase> const
  get_ekinstance() const {
    return m_ekinstance;
  }

  [[nodiscard]] std::shared_ptr<LatticeWalberla> get_lattice() {
    return m_lattice;
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override;

  [[nodiscard]] double get_conversion_factor_density() const noexcept {
    return m_conv_density;
  }
  [[nodiscard]] double get_conversion_factor_flux() const noexcept {
    return m_conv_flux;
  }

  LatticeModel::units_map get_latice_to_md_units_conversion() const {
    return {
        {"density", 1. / m_conv_density},
        {"flux", 1. / m_conv_flux},
    };
  }

private:
  std::shared_ptr<::EKinWalberlaBase> m_ekinstance;

  std::shared_ptr<LatticeWalberla> m_lattice;

  void load_checkpoint(std::string const &filename, int mode);
  void save_checkpoint(std::string const &filename, int mode);
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
