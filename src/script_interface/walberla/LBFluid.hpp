/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "LatticeModel.hpp"
#include "LatticeWalberla.hpp"
#include "VTKHandle.hpp"

#include "core/lb/LBWalberla.hpp"
#include "core/lb/Solver.hpp"
#include "core/system/System.hpp"

#include <script_interface/ScriptInterface.hpp>

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

class LBVTKHandle : public VTKHandleBase<::LBWalberlaBase> {
  static std::unordered_map<std::string, int> const obs_map;

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

inline void
lb_throw_if_expired(std::optional<ResourceObserver> const &mpi_obs) {
  if (not(mpi_obs and mpi_obs->is_valid())) {
    throw std::runtime_error(
        "the MPI Cartesian communicator of this LB object has expired");
  }
}

class LBFluid : public LatticeModel<::LBWalberlaBase, LBVTKHandle> {
protected:
  using Base = LatticeModel<::LBWalberlaBase, LBVTKHandle>;
  std::shared_ptr<::LB::LBWalberlaParams> m_lb_params;
  std::optional<ResourceObserver> m_mpi_cart_comm_observer;
  bool m_is_active;
  double m_conv_dist;
  double m_conv_visc;
  double m_conv_dens;
  double m_conv_speed;
  double m_conv_press;
  double m_conv_force;
  double m_conv_force_dens;
  double m_conv_energy;

public:
  LBFluid() {
    add_parameters({
        {"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
        {"single_precision", AutoParameter::read_only,
         [this]() { return not m_instance->is_double_precision(); }},
        {"gpu", AutoParameter::read_only,
         [this]() { return m_instance->is_gpu(); }},
        {"is_active", AutoParameter::read_only,
         [this]() { return m_is_active; }},
        {"agrid", AutoParameter::read_only,
         [this]() { return m_lb_params->get_agrid(); }},
        {"tau", AutoParameter::read_only,
         [this]() { return m_lb_params->get_tau(); }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_instance->get_lattice().get_grid_dimensions(); }},
        {"kT", AutoParameter::read_only,
         [this]() { return m_instance->get_kT() / m_conv_energy; }},
        {"seed", AutoParameter::read_only,
         [this]() { return static_cast<int>(m_instance->get_seed()); }},
        {"rng_state",
         [this](Variant const &v) {
           auto const rng_state = get_value<int>(v);
           context()->parallel_try_catch([&]() {
             if (rng_state < 0) {
               throw std::domain_error("Parameter 'rng_state' must be >= 0");
             }
             m_instance->set_rng_state(static_cast<uint64_t>(rng_state));
           });
         },
         [this]() {
           auto const opt = m_instance->get_rng_state();
           return (opt) ? Variant{static_cast<int>(*opt)} : Variant{None{}};
         }},
        {"density", AutoParameter::read_only,
         [this]() { return m_instance->get_density() / m_conv_dens; }},
        {"kinematic_viscosity",
         [this](Variant const &v) {
           auto const visc = m_conv_visc * get_value<double>(v);
           m_instance->set_viscosity(visc);
         },
         [this]() { return m_instance->get_viscosity() / m_conv_visc; }},
        {"ext_force_density",
         [this](Variant const &v) {
           auto const ext_f = m_conv_force_dens * get_value<Utils::Vector3d>(v);
           m_instance->set_external_force(ext_f);
         },
         [this]() {
           return m_instance->get_external_force() / m_conv_force_dens;
         }},
        {"vtk_writers", AutoParameter::read_only,
         [this]() { return serialize_vtk_writers(); }},
    });
  }

  void do_construct(VariantMap const &params) override;

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  [[nodiscard]] auto get_lb_fluid() const { return m_instance; }
  [[nodiscard]] auto get_lb_params() const { return m_lb_params; }
  [[nodiscard]] auto get_mpi_cart_comm_observer() const {
    return m_mpi_cart_comm_observer;
  }

  ::LatticeModel::units_map
  get_lattice_to_md_units_conversion() const override {
    return {
        {"density", 1. / m_conv_dens},
        {"velocity", 1. / m_conv_speed},
        {"pressure", 1. / m_conv_press},
    };
  }

protected:
  void make_instance(VariantMap const &params) override;

private:
  void load_checkpoint(std::filesystem::path const &path, int mode);
  void save_checkpoint(std::filesystem::path const &path, int mode);
  std::vector<Variant> get_average_pressure_tensor() const;
  Variant get_boundary_force_from_shape(std::vector<int> const &raster) const;
  Variant get_boundary_force() const;
  Variant get_interpolated_velocity(Utils::Vector3d const &pos) const;
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
