/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_FLUIDWALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_FLUIDWALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaBase.hpp>

#include "LatticeWalberla.hpp"

#include "core/communication.hpp"
#include "core/grid_based_algorithms/lb_interface.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "core/integrate.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <memory>
#include <string>
#include <tuple>

namespace ScriptInterface::walberla {

class FluidWalberla : public AutoParameters<FluidWalberla> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  std::shared_ptr<::LBWalberlaParams> m_lb_params;
  bool m_is_single_precision;
  bool m_is_active;
  double m_conv_visc;
  double m_conv_temp;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  std::tuple<double, double, double, int, Utils::Vector3d> m_ctor_params;

public:
  FluidWalberla() {
    add_parameters({
        {"is_single_precision", AutoParameter::read_only,
         [this]() { return m_is_single_precision; }},
        {"is_active", AutoParameter::read_only,
         [this]() { return m_is_active; }},
        {"is_initialized", AutoParameter::read_only,
         [this]() { return m_lb_fluid != nullptr; }},
        {"agrid", AutoParameter::read_only,
         [this]() { return m_lb_params->get_agrid(); }},
        {"tau", AutoParameter::read_only,
         [this]() { return m_lb_params->get_tau(); }},
        {"shape", AutoParameter::read_only,
         [this]() {
           return (m_lb_fluid) ? m_lb_fluid->get_grid_dimensions()
                               : Utils::Vector3i::broadcast(-1);
         }},
        {"kT", AutoParameter::read_only,
         [this]() { return std::get<2>(m_ctor_params) / m_conv_temp; }},
        {"seed", AutoParameter::read_only,
         [this]() { return std::get<3>(m_ctor_params); }},
        {"rng_state",
         [this](const Variant &v) {
           try {
             if (m_lb_fluid)
               m_lb_fluid->set_rng_state(
                   static_cast<uint64_t>(get_value<int>(v)));
             else
               throw std::runtime_error(
                   "Cannot set 'rng_state' before walberla is initialized");
           } catch (const std::exception &e) {
             if (this_node == 0) {
               throw;
             }
           }
         },
         [this]() {
           try {
             if (m_lb_fluid)
               return static_cast<int>(m_lb_fluid->get_rng_state());
             else
               throw std::runtime_error(
                   "Cannot get 'rng_state' before walberla is initialized");
           } catch (const std::exception &e) {
             if (this_node == 0) {
               throw;
             }
           }
           return -1;
         }},
        {"density", AutoParameter::read_only,
         [this]() {
           return ((m_lb_fluid) ? m_lb_fluid->get_density()
                                : std::get<1>(m_ctor_params)) /
                  m_conv_dens;
         }},
        {"viscosity",
         [this](const Variant &v) {
           auto const visc = m_conv_visc * get_value<double>(v);
           if (m_lb_fluid)
             m_lb_fluid->set_viscosity(visc);
           else
             std::get<0>(m_ctor_params) = visc;
         },
         [this]() {
           return ((m_lb_fluid) ? m_lb_fluid->get_viscosity()
                                : std::get<0>(m_ctor_params)) /
                  m_conv_visc;
         }},
        {"ext_force_density",
         [this](const Variant &v) {
           auto const ext_f = m_conv_force * get_value<Utils::Vector3d>(v);
           if (m_lb_fluid)
             m_lb_fluid->set_external_force(ext_f);
           else
             std::get<4>(m_ctor_params) = ext_f;
         },
         [this]() {
           return ((m_lb_fluid) ? m_lb_fluid->get_external_force()
                                : std::get<4>(m_ctor_params)) /
                  m_conv_force;
         }},
    });
  }

  void do_construct(VariantMap const &params) override {
    // construction of the LB object is deferred until the first
    // activation, because the box length and MD time step can change
    // freely when the LB object is not active yet
    auto const tau = get_value<double>(params, "tau");
    auto const agrid = get_value<double>(params, "agrid");
    m_conv_visc = Utils::int_pow<1>(tau) / Utils::int_pow<2>(agrid);
    m_conv_temp = Utils::int_pow<2>(tau) / Utils::int_pow<2>(agrid);
    m_conv_dens = Utils::int_pow<3>(agrid);
    m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
    m_conv_force = Utils::int_pow<2>(agrid) * Utils::int_pow<2>(tau);
    m_lb_params = std::make_shared<::LBWalberlaParams>(agrid, tau);
    m_ctor_params = {
        m_conv_visc * get_value<double>(params, "viscosity"),
        m_conv_dens * get_value<double>(params, "density"),
        m_conv_temp * get_value<double>(params, "kT"),
        get_value<int>(params, "seed"),
        m_conv_force * get_value<Utils::Vector3d>(params, "ext_force_density")};
    m_is_active = false;
    m_is_single_precision = get_value<bool>(params, "single_precision");
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "instantiate") {
      auto const [lb_visc, lb_dens, lb_kT, lb_seed, ext_f] = m_ctor_params;
      auto const lb_lattice =
          get_value<std::shared_ptr<LatticeWalberla>>(params, "lattice")
              ->lattice();
      m_lb_fluid = init_lb_walberla(*lb_lattice, *m_lb_params, lb_visc, lb_dens,
                                    lb_kT, lb_seed, m_is_single_precision);
      m_lb_fluid->set_external_force(ext_f);
    }
    if (name == "activate") {
      auto const fail = activate_lb_walberla(m_lb_fluid, m_lb_params);
      if (not fail) {
        m_is_active = true;
      }
    }
    if (name == "deactivate") {
      deactivate_lb_walberla();
      m_is_active = false;
    }

    return {};
  }

  /** Non-owning pointer to the LB fluid. */
  std::weak_ptr<::LBWalberlaBase> lb_fluid() { return m_lb_fluid; }

  /** Non-owning pointer to the LB parameters. */
  std::weak_ptr<::LBWalberlaParams> lb_params() { return m_lb_params; }
};

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
