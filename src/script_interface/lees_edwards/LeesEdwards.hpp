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
#ifndef SCRIPT_INTERFACE_LEES_EDWARDS_LEES_EDWARDS_HPP
#define SCRIPT_INTERFACE_LEES_EDWARDS_LEES_EDWARDS_HPP

#include "Protocol.hpp"

#include "core/grid.hpp"
#include "core/grid_based_algorithms/lb_interface.hpp"
#include "core/lees_edwards/LeesEdwardsBC.hpp"
#include "core/lees_edwards/lees_edwards.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <stdexcept>

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards> {
  std::shared_ptr<Protocol> m_protocol;
  LeesEdwardsBC const &m_lebc = ::box_geo.lees_edwards_bc();

public:
  LeesEdwards() : m_protocol{nullptr} {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            if (is_none(value)) {
              context()->parallel_try_catch([]() {
                auto constexpr invalid_dir = LeesEdwardsBC::invalid_dir;
                LB::lebc_sanity_checks(invalid_dir, invalid_dir);
              });
              m_protocol = nullptr;
              ::box_geo.set_lees_edwards_bc(LeesEdwardsBC{});
              ::LeesEdwards::unset_protocol();
              return;
            }
            context()->parallel_try_catch([this]() {
              auto constexpr invalid_dir = LeesEdwardsBC::invalid_dir;
              if (m_lebc.shear_direction == invalid_dir or
                  m_lebc.shear_plane_normal == invalid_dir) {
                throw std::runtime_error(
                    "Parameters 'shear_plane_normal' and 'shear_direction' "
                    "must be initialized together with 'protocol' on first "
                    "activation via set_boundary_conditions()");
              }
            });
            m_protocol = get_value<std::shared_ptr<Protocol>>(value);
            ::LeesEdwards::set_protocol(m_protocol->protocol());
          },
          [this]() {
            if (m_protocol)
              return make_variant(m_protocol);
            return make_variant(none);
          }},
         {"shear_velocity", AutoParameter::read_only,
          [this]() { return m_lebc.shear_velocity; }},
         {"pos_offset", AutoParameter::read_only,
          [this]() { return m_lebc.pos_offset; }},
         {"shear_direction", AutoParameter::read_only,
          [this]() { return get_shear_name(m_lebc.shear_direction); }},
         {"shear_plane_normal", AutoParameter::read_only,
          [this]() { return get_shear_name(m_lebc.shear_plane_normal); }}});
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_boundary_conditions") {
      context()->parallel_try_catch([this, &params]() {
        auto const protocol = params.at("protocol");
        if (is_none(protocol)) {
          do_set_parameter("protocol", protocol);
          return;
        }
        // check input arguments
        m_protocol = get_value<std::shared_ptr<Protocol>>(protocol);
        auto const shear_direction = get_shear_axis(params, "shear_direction");
        auto const shear_plane_normal =
            get_shear_axis(params, "shear_plane_normal");
        if (shear_plane_normal == shear_direction) {
          throw std::invalid_argument("Parameters 'shear_direction' and "
                                      "'shear_plane_normal' must differ");
        }
        LB::lebc_sanity_checks(shear_direction, shear_plane_normal);
        // update box geometry and cell structure
        ::box_geo.set_lees_edwards_bc(
            LeesEdwardsBC{0., 0., shear_direction, shear_plane_normal});
        ::LeesEdwards::set_protocol(m_protocol->protocol());
      });
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    if (not params.empty()) {
      do_call_method("set_boundary_conditions", params);
    }
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
};

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
