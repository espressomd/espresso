/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include "core/grid_based_algorithms/lb_walberla_instance.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaNodeState.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

class Fluid : public AutoParameters<Fluid> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  std::shared_ptr<::LBWalberlaParams> m_lb_params;
  bool m_is_single_precision;
  bool m_is_active;
  double m_conv_dist;
  double m_conv_visc;
  double m_conv_temp;
  double m_conv_dens;
  double m_conv_speed;
  double m_conv_press;
  double m_conv_force;
  double m_conv_force_dens;
  int m_seed;

public:
  Fluid() {
    add_parameters({
        {"is_single_precision", AutoParameter::read_only,
         [this]() { return m_is_single_precision; }},
        {"is_active", AutoParameter::read_only,
         [this]() { return m_is_active; }},
        {"agrid", AutoParameter::read_only,
         [this]() { return m_lb_params->get_agrid(); }},
        {"tau", AutoParameter::read_only,
         [this]() { return m_lb_params->get_tau(); }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_lb_fluid->get_lattice().get_grid_dimensions(); }},
        {"kT", AutoParameter::read_only,
         [this]() { return m_lb_fluid->get_kT() / m_conv_temp; }},
        {"seed", AutoParameter::read_only, [this]() { return m_seed; }},
        {"rng_state",
         [this](const Variant &v) {
           context()->parallel_try_catch([&]() {
             m_lb_fluid->set_rng_state(
                 static_cast<uint64_t>(get_value<int>(v)));
           });
         },
         [this]() { return static_cast<int>(m_lb_fluid->get_rng_state()); }},
        {"density", AutoParameter::read_only,
         [this]() { return m_lb_fluid->get_density() / m_conv_dens; }},
        {"kinematic_viscosity",
         [this](const Variant &v) {
           auto const visc = m_conv_visc * get_value<double>(v);
           m_lb_fluid->set_viscosity(visc);
         },
         [this]() { return m_lb_fluid->get_viscosity() / m_conv_visc; }},
        {"ext_force_density",
         [this](const Variant &v) {
           auto const ext_f = m_conv_force_dens * get_value<Utils::Vector3d>(v);
           m_lb_fluid->set_external_force(ext_f);
         },
         [this]() {
           return m_lb_fluid->get_external_force() / m_conv_force_dens;
         }},
    });
  }

  void do_construct(VariantMap const &params) override {
    auto const lb_lattice_si =
        get_value<std::shared_ptr<LatticeWalberla>>(params, "lattice");
    auto const tau = get_value<double>(params, "tau");
    auto const agrid = get_value<double>(lb_lattice_si->get_parameter("agrid"));
    m_conv_dist = 1. / agrid;
    m_conv_visc = Utils::int_pow<1>(tau) / Utils::int_pow<2>(agrid);
    m_conv_temp = Utils::int_pow<2>(tau) / Utils::int_pow<2>(agrid);
    m_conv_dens = Utils::int_pow<3>(agrid);
    m_conv_speed = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    m_conv_press = Utils::int_pow<2>(tau) * Utils::int_pow<1>(agrid);
    m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
    m_conv_force_dens = Utils::int_pow<2>(tau) * Utils::int_pow<2>(agrid);
    m_lb_params = std::make_shared<::LBWalberlaParams>(agrid, tau);
    m_is_active = false;
    m_seed = get_value<int>(params, "seed");
    auto const lb_lattice = lb_lattice_si->lattice();
    auto const lb_visc =
        m_conv_visc * get_value<double>(params, "kinematic_viscosity");
    auto const lb_dens = m_conv_dens * get_value<double>(params, "density");
    auto const lb_temp = m_conv_temp * get_value<double>(params, "kT");
    auto const ext_f = m_conv_force_dens *
                       get_value<Utils::Vector3d>(params, "ext_force_density");
    m_is_single_precision = get_value<bool>(params, "single_precision");
    m_lb_fluid = init_lb_walberla(lb_lattice, *m_lb_params, lb_visc, lb_dens,
                                  lb_temp, m_seed, m_is_single_precision);
    if (m_lb_fluid) {
      m_lb_fluid->set_external_force(ext_f);
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  /** Non-owning pointer to the LB fluid. */
  std::weak_ptr<::LBWalberlaBase> lb_fluid() { return m_lb_fluid; }

  /** Non-owning pointer to the LB parameters. */
  std::weak_ptr<::LBWalberlaParams> lb_params() { return m_lb_params; }

  LatticeModel::units_map get_latice_to_md_units_conversion() const {
    return {
        {"density", 1. / m_conv_dens},
        {"velocity", 1. / m_conv_speed},
        {"pressure", 1. / m_conv_press},
    };
  }

private:
  void load_checkpoint(std::string const &filename, int mode);
  void save_checkpoint(std::string const &filename, int mode);
  std::vector<Variant> get_average_pressure_tensor() const;
  Variant get_interpolated_velocity(Utils::Vector3d const &pos) const;
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
