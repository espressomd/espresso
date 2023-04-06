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

#include "FluidWalberla.hpp"

#include "LatticeIndices.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

class FluidSlice : public AutoParameters<FluidSlice, LatticeIndices> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  std::shared_ptr<FluidWalberla> m_lb_sip;
  Utils::Vector3i m_slice_lower_corner;
  Utils::Vector3i m_slice_upper_corner;
  std::vector<int> m_shape;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;
  std::unordered_map<std::string, std::vector<int>> m_shape_val;

public:
  FluidSlice() { add_parameters({}); }

  void do_construct(VariantMap const &params) override {
    m_lb_sip = get_value<std::shared_ptr<FluidWalberla>>(params, "parent_sip");
    m_lb_fluid = m_lb_sip->lb_fluid().lock();
    assert(m_lb_fluid);
    auto const &lb_params = m_lb_sip->lb_params().lock();
    assert(lb_params);
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
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

private:
  auto get_sentinel_index(::LatticeWalberla const &lattice) const {
    return -(static_cast<int>(lattice.get_ghost_layers()) + 1);
  }

  auto get_slices_bounding_boxes() const {
    auto const &lattice = m_lb_fluid->get_lattice();
    auto const &slice_lower_corner = m_slice_lower_corner;
    auto const &slice_upper_corner = m_slice_upper_corner;
    assert(slice_upper_corner <= lattice.get_grid_dimensions());
    assert(slice_lower_corner >= Utils::Vector3i::broadcast(0));
    auto const sentinel = get_sentinel_index(lattice);
    auto [local_lower_corner, local_upper_corner] =
        lattice.get_local_grid_range();
    for (auto const i : {0, 1, 2}) {
      if (local_lower_corner[i] >= slice_upper_corner[i] or
          slice_lower_corner[i] >= local_upper_corner[i]) {
        local_lower_corner[i] = sentinel;
        local_upper_corner[i] = sentinel;
      } else {
        if (slice_lower_corner[i] > local_lower_corner[i]) {
          local_lower_corner[i] = slice_lower_corner[i];
        }
        if (slice_upper_corner[i] < local_upper_corner[i]) {
          local_upper_corner[i] = slice_upper_corner[i];
        }
      }
    }
    return std::make_tuple(slice_lower_corner, slice_upper_corner,
                           local_lower_corner, local_upper_corner);
  }

private:
  template <typename T>
  Variant gather_3d(VariantMap const &params, std::vector<int> const &data_dims,
                    std::vector<T> (::LBWalberlaBase::*getter)(
                        Utils::Vector3i const &, Utils::Vector3i const &) const,
                    double units_conversion = 1.) const;

  template <typename T>
  void scatter_3d(VariantMap const &params, std::vector<int> const &data_dims,
                  void (::LBWalberlaBase::*setter)(Utils::Vector3i const &,
                                                   Utils::Vector3i const &,
                                                   std::vector<T> const &),
                  double units_conversion = 1.);
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
