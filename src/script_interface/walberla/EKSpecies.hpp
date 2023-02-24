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

#include "LatticeWalberla.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include "walberla_bridge/electrokinetics/ek_walberla_init.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <walberla_bridge/electrokinetics/EKWalberlaNodeState.hpp>

#include <utils/math/int_pow.hpp>

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

  void do_construct(VariantMap const &args) override {
    m_lattice = get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice");
    const auto single_precision =
        get_value_or<bool>(args, "single_precision", false);

    // unit conversions
    auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
    m_tau = get_value<double>(args, "tau");

    m_conv_diffusion = m_tau / Utils::int_pow<2>(agrid);
    m_conv_ext_efield = Utils::int_pow<2>(m_tau) / agrid;
    m_conv_density = Utils::int_pow<3>(agrid);
    m_conv_flux = m_tau * Utils::int_pow<2>(agrid);

    auto const ek_diffusion =
        get_value<double>(args, "diffusion") * m_conv_diffusion;
    auto const ek_ext_efield =
        get_value<Utils::Vector3d>(args, "ext_efield") * m_conv_ext_efield;
    auto const ek_density = m_density =
        get_value<double>(args, "density") * m_conv_density;

    m_ekinstance = new_ek_walberla(
        m_lattice->lattice(), ek_diffusion, get_value<double>(args, "kT"),
        get_value<double>(args, "valency"), ek_ext_efield, ek_density,
        get_value<bool>(args, "advection"),
        get_value<bool>(args, "friction_coupling"), single_precision);
  }

  [[nodiscard]] std::shared_ptr<::EKinWalberlaBase> const get_ekinstance() const {
    return m_ekinstance;
  }

  [[nodiscard]] std::shared_ptr<LatticeWalberla> get_lattice() {
    return m_lattice;
  }

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "update_flux_boundary_from_shape") {
      auto value_view =
          get_value<std::vector<double>>(parameters, "value_view");
      std::transform(value_view.begin(), value_view.end(), value_view.begin(),
                     [this](double v) { return v * m_conv_flux; });

      m_ekinstance->update_flux_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"), value_view);
      return none;
    }
    if (method == "update_density_boundary_from_shape") {
      auto value_view =
          get_value<std::vector<double>>(parameters, "value_view");
      std::transform(value_view.begin(), value_view.end(), value_view.begin(),
                     [this](double v) { return v * m_conv_density; });
      m_ekinstance->update_density_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"), value_view);
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
    if (method == "save_checkpoint") {
      auto const path = get_value<std::string>(parameters, "path");
      auto const mode = get_value<int>(parameters, "mode");
      save_checkpoint(path, mode);
    }
    if (method == "load_checkpoint") {
      auto const path = get_value<std::string>(parameters, "path");
      auto const mode = get_value<int>(parameters, "mode");
      load_checkpoint(path, mode);
    }
    return none;
  }

  [[nodiscard]] double get_conversion_factor_density() const noexcept {
    return m_conv_density;
  }
  [[nodiscard]] double get_conversion_factor_flux() const noexcept {
    return m_conv_flux;
  }

private:
  std::shared_ptr<::EKinWalberlaBase> m_ekinstance;

  std::shared_ptr<LatticeWalberla> m_lattice;

  void load_checkpoint(std::string const &filename, int mode);
  void save_checkpoint(std::string const &filename, int mode);
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
