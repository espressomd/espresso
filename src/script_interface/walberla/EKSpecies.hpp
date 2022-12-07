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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_SPECIES_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_SPECIES_HPP

#include "config/config.hpp"

#ifdef WALBERLA

#include "walberla_bridge/LatticeWalberla.hpp"

#include "walberla_bridge/electrokinetics/ek_walberla_init.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

namespace ScriptInterface::walberla {

class EKSpecies : public AutoParameters<EKSpecies> {
public:
  void do_construct(VariantMap const &args) override {
    m_lattice = get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice");
    const auto single_precision =
        get_value_or<bool>(args, "single_precision", false);
    m_ekinstance = new_ek_walberla(
        m_lattice->lattice(), get_value<double>(args, "diffusion"),
        get_value<double>(args, "kT"), get_value<double>(args, "valency"),
        get_value<Utils::Vector3d>(args, "ext_efield"),
        get_value<double>(args, "density"), get_value<bool>(args, "advection"),
        get_value<bool>(args, "friction_coupling"), single_precision);

    add_parameters(
        {{"diffusion",
          [this](Variant const &v) {
            m_ekinstance->set_diffusion(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_diffusion(); }},
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
            m_ekinstance->set_ext_efield(get_value<Utils::Vector3d>(v));
          },
          [this]() { return m_ekinstance->get_ext_efield(); }},
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
         {"shape", AutoParameter::read_only,
          [this]() {
            return m_ekinstance->get_lattice().get_grid_dimensions();
          }},
         {"lattice", AutoParameter::read_only,
          [this]() { return m_lattice; }}});
  }

  [[nodiscard]] std::shared_ptr<::EKinWalberlaBase> get_ekinstance() {
    return m_ekinstance;
  }

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "update_flux_boundary_from_shape") {
      m_ekinstance->update_flux_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    if (method == "update_density_boundary_from_shape") {
      m_ekinstance->update_density_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    if (method == "clear_flux_boundaries") {
      m_ekinstance->clear_flux_boundaries();
      return none;
    }
    if (method == "clear_density_boundaries") {
      m_ekinstance->clear_density_boundaries();
      return none;
    }
    return none;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::EKinWalberlaBase> m_ekinstance;

  std::shared_ptr<LatticeWalberla> m_lattice;
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
