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

#include "Protocol.hpp"

#include "core/BoxGeometry.hpp"
#include "core/lees_edwards/LeesEdwardsBC.hpp"
#include "core/lees_edwards/lees_edwards.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/system/Leaf.hpp"

#include <memory>
#include <stdexcept>

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards, System::Leaf> {
  std::shared_ptr<::LeesEdwards::LeesEdwards> m_lees_edwards;
  std::shared_ptr<Protocol> m_protocol;
  std::unique_ptr<VariantMap> m_params;

public:
  LeesEdwards() {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            auto &system = get_system();
            if (is_none(value)) {
              context()->parallel_try_catch([&]() {
                auto constexpr invalid_dir = LeesEdwardsBC::invalid_dir;
                system.lb.lebc_sanity_checks(invalid_dir, invalid_dir);
              });
              m_protocol = nullptr;
              system.box_geo->set_lees_edwards_bc(LeesEdwardsBC{});
              m_lees_edwards->unset_protocol();
              system.on_lees_edwards_change();
              return;
            }
            context()->parallel_try_catch([this]() {
              auto constexpr invalid_dir = LeesEdwardsBC::invalid_dir;
              auto const &lebc = get_lebc();
              if (lebc.shear_direction == invalid_dir or
                  lebc.shear_plane_normal == invalid_dir) {
                throw std::runtime_error(
                    "Parameters 'shear_plane_normal' and 'shear_direction' "
                    "must be initialized together with 'protocol' on first "
                    "activation via set_boundary_conditions()");
              }
            });
            auto const protocol = get_value<std::shared_ptr<Protocol>>(value);
            context()->parallel_try_catch([&]() {
              try {
                set_protocol(protocol);
              } catch (...) {
                set_protocol(m_protocol);
                throw;
              }
            });
          },
          [this]() {
            if (m_protocol)
              return make_variant(m_protocol);
            return make_variant(none);
          }},
         {"shear_velocity", AutoParameter::read_only,
          [this]() { return get_lebc().shear_velocity; }},
         {"pos_offset", AutoParameter::read_only,
          [this]() { return get_lebc().pos_offset; }},
         {"shear_direction", AutoParameter::read_only,
          [this]() { return get_shear_name(get_lebc().shear_direction); }},
         {"shear_plane_normal", AutoParameter::read_only,
          [this]() { return get_shear_name(get_lebc().shear_plane_normal); }}});
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_boundary_conditions") {
      context()->parallel_try_catch([this, &params]() {
        auto const variant = params.at("protocol");
        if (is_none(variant)) {
          do_set_parameter("protocol", variant);
          return;
        }
        // check input arguments
        auto const protocol = get_value<std::shared_ptr<Protocol>>(variant);
        auto const shear_direction = get_shear_axis(params, "shear_direction");
        auto const shear_plane_normal =
            get_shear_axis(params, "shear_plane_normal");
        if (shear_plane_normal == shear_direction) {
          throw std::invalid_argument("Parameters 'shear_direction' and "
                                      "'shear_plane_normal' must differ");
        }
        auto &system = get_system();
        system.lb.lebc_sanity_checks(shear_direction, shear_plane_normal);
        // update box geometry and cell structure
        auto const old_shear_direction = get_lebc().shear_direction;
        auto const old_shear_plane_normal = get_lebc().shear_plane_normal;
        try {
          system.box_geo->set_lees_edwards_bc(
              LeesEdwardsBC{0., 0., shear_direction, shear_plane_normal});
          set_protocol(protocol);
        } catch (...) {
          system.box_geo->set_lees_edwards_bc(LeesEdwardsBC{
              0., 0., old_shear_direction, old_shear_plane_normal});
          set_protocol(m_protocol);
          throw;
        }
      });
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

private:
  unsigned int get_shear_axis(VariantMap const &params, std::string name) {
    auto const value = get_value<std::string>(params, name);
    if (value == "x") {
      return 0u;
    }
    if (value == "y") {
      return 1u;
    }
    if (value == "z") {
      return 2u;
    }
    throw std::invalid_argument("Parameter '" + name + "' is invalid");
  }

  Variant get_shear_name(unsigned int axis) {
    if (axis == 0u) {
      return {std::string("x")};
    }
    if (axis == 1u) {
      return {std::string("y")};
    }
    if (axis == 2u) {
      return {std::string("z")};
    }
    return {none};
  }

  LeesEdwardsBC const &get_lebc() const {
    return get_system().box_geo->lees_edwards_bc();
  }

  void on_bind_system(::System::System &system) override {
    m_lees_edwards = system.lees_edwards;
    auto const &params = *m_params;
    if (not params.empty()) {
      do_call_method("set_boundary_conditions", params);
    }
    m_params.reset();
  }

  void set_protocol(std::shared_ptr<Protocol> const &protocol) {
    auto &system = get_system();
    if (protocol) {
      m_lees_edwards->set_protocol(protocol->protocol());
    } else {
      system.box_geo->set_lees_edwards_bc(LeesEdwardsBC{});
      m_lees_edwards->unset_protocol();
    }
    system.on_lees_edwards_change();
    m_protocol = protocol;
  }
};

} // namespace LeesEdwards
} // namespace ScriptInterface
