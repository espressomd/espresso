/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "LBFluid.hpp"

#include "LatticeSlice.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <cassert>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

using VelocityBounceBackType = boost::optional<Utils::Vector3d>;

struct LBFieldSerializer {
  template <typename T> static Variant serialize(std::vector<T> const &values) {
    if constexpr (std::is_same_v<T, VelocityBounceBackType>) {
      std::vector<Variant> vec;
      vec.reserve(values.size());
      for (auto const &opt : values) {
        if (opt) {
          vec.emplace_back(Variant{*opt});
        } else {
          vec.emplace_back(Variant{None{}});
        }
      }
      return {vec};
    } else if constexpr (std::is_same_v<T, int> or std::is_same_v<T, double>) {
      return {values};
    } else {
      return make_vector_of_variants(values);
    }
  }

  template <typename T>
  static std::vector<T> deserialize(Variant const &variant) {
    std::vector<T> values;
    if constexpr (std::is_same_v<T, VelocityBounceBackType>) {
      auto const vector_variants = get_value<std::vector<Variant>>(variant);
      for (auto const &value : vector_variants) {
        if (is_none(value)) {
          values.emplace_back(boost::none);
        } else {
          values.emplace_back(get_value<Utils::Vector3d>(value));
        }
      }
    } else if constexpr (std::is_same_v<T, double>) {
      if (is_type<std::vector<int>>(variant)) {
        auto const values_int = get_value<std::vector<int>>(variant);
        values.reserve(values_int.size());
        for (auto const val : values_int) {
          values.emplace_back(static_cast<double>(val));
        }
      } else {
        values = get_value<std::vector<T>>(variant);
      }
    } else {
      values = get_value<std::vector<T>>(variant);
    }
    return values;
  }
};

class LBFluidSlice : public LatticeSlice<LBFieldSerializer> {
  using LatticeModel = ::LBWalberlaBase;
  std::shared_ptr<LatticeModel> m_lb_fluid;
  std::shared_ptr<LBFluid> m_lb_sip;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;
  std::unordered_map<std::string, std::vector<int>> m_shape_val;

public:
  void do_construct(VariantMap const &params) override {
    m_lb_sip = get_value<std::shared_ptr<LBFluid>>(params, "parent_sip");
    m_lb_fluid = m_lb_sip->get_lb_fluid();
    auto const lb_params = m_lb_sip->get_lb_params();
    auto const tau = lb_params->get_tau();
    auto const agrid = lb_params->get_agrid();
    m_conv_dens = Utils::int_pow<3>(agrid);
    m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
    m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
    m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    m_shape = get_value<std::vector<int>>(params, "shape");
    m_slice_lower_corner =
        get_value<Utils::Vector3i>(params, "slice_lower_corner");
    m_slice_upper_corner =
        get_value<Utils::Vector3i>(params, "slice_upper_corner");
    m_shape_val["density"] = {1};
    m_shape_val["population"] = {static_cast<int>(m_lb_fluid->stencil_size())};
    m_shape_val["velocity"] = {3};
    m_shape_val["velocity_at_boundary"] = {1};
    m_shape_val["is_boundary"] = {1};
    m_shape_val["last_applied_force"] = {3};
    m_shape_val["pressure_tensor"] = {3, 3};
    m_shape_val["pressure_tensor_neq"] = {3, 3};
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  ::LatticeWalberla const &get_lattice() const override {
    return m_lb_fluid->get_lattice();
  }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
