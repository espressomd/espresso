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

#include "LatticeModel.hpp"
#include "LatticeWalberla.hpp"
#include "VTKHandle.hpp"

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include <script_interface/ScriptInterface.hpp>

#include <utils/math/int_pow.hpp>

#include <memory>
#include <string>

namespace ScriptInterface::walberla {

class EKVTKHandle;

class EKSpecies : public LatticeModel<::EKinWalberlaBase, EKVTKHandle> {
  using Base = LatticeModel<::EKinWalberlaBase, EKVTKHandle>;
  double m_conv_diffusion;
  double m_conv_ext_efield;
  double m_conv_energy;
  double m_conv_density;
  double m_conv_flux;
  double m_tau;
  double m_density;

public:
  EKSpecies() {
    add_parameters({
        {"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
        {"diffusion",
         [this](Variant const &v) {
           m_instance->set_diffusion(get_value<double>(v) * m_conv_diffusion);
         },
         [this]() { return m_instance->get_diffusion() / m_conv_diffusion; }},
        {"kT",
         [this](Variant const &v) {
           m_instance->set_kT(get_value<double>(v) * m_conv_energy);
         },
         [this]() { return m_instance->get_kT() / m_conv_energy; }},
        {"valency",
         [this](Variant const &v) {
           m_instance->set_valency(get_value<double>(v));
         },
         [this]() { return m_instance->get_valency(); }},
        {"ext_efield",
         [this](Variant const &v) {
           m_instance->set_ext_efield(get_value<Utils::Vector3d>(v) *
                                      m_conv_ext_efield);
         },
         [this]() { return m_instance->get_ext_efield() / m_conv_ext_efield; }},
        {"advection",
         [this](Variant const &v) {
           m_instance->set_advection(get_value<bool>(v));
         },
         [this]() { return m_instance->get_advection(); }},
        {"friction_coupling",
         [this](Variant const &v) {
           m_instance->set_friction_coupling(get_value<bool>(v));
         },
         [this]() { return m_instance->get_friction_coupling(); }},
        {"single_precision", AutoParameter::read_only,
         [this]() { return not m_instance->is_double_precision(); }},
        {"tau", AutoParameter::read_only, [this]() { return m_tau; }},
        {"density", AutoParameter::read_only,
         [this]() { return m_density / m_conv_density; }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_instance->get_lattice().get_grid_dimensions(); }},
        {"vtk_writers", AutoParameter::read_only,
         [this]() { return serialize_vtk_writers(); }},
    });
  }

  void do_construct(VariantMap const &args) override;

  [[nodiscard]] auto get_ekinstance() const { return m_instance; }
  [[nodiscard]] auto get_lattice() const { return m_lattice; }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override;

  [[nodiscard]] auto get_conversion_factor_density() const noexcept {
    return m_conv_density;
  }
  [[nodiscard]] auto get_conversion_factor_flux() const noexcept {
    return m_conv_flux;
  }

  ::LatticeModel::units_map get_latice_to_md_units_conversion() const override {
    return {
        {"density", 1. / m_conv_density},
        {"flux", 1. / m_conv_flux},
    };
  }

private:
  void load_checkpoint(std::string const &filename, int mode);
  void save_checkpoint(std::string const &filename, int mode);
};

class EKVTKHandle : public VTKHandleBase<::EKinWalberlaBase> {
  static std::unordered_map<std::string, int> const obs_map;

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
