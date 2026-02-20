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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "LatticeModel.hpp"
#include "LatticeWalberla.hpp"
#include "VTKHandle.hpp"

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>

#include <script_interface/ScriptInterface.hpp>

#include <utils/math/int_pow.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace ScriptInterface::walberla {

class EKVTKHandle : public VTKHandleBase<::EKinWalberlaBase> {
  static std::unordered_map<std::string, int> const obs_map;

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

inline void
ek_throw_if_expired(std::optional<ResourceObserver> const &mpi_obs) {
  if (not(mpi_obs and mpi_obs->is_valid())) {
    throw std::runtime_error(
        "the MPI Cartesian communicator of this EK object has expired");
  }
}

class EKSpecies : public LatticeModel<::EKinWalberlaBase, EKVTKHandle> {
protected:
  using Base = LatticeModel<::EKinWalberlaBase, EKVTKHandle>;
  std::optional<ResourceObserver> m_mpi_cart_comm_observer;
  double m_conv_diffusion;
  double m_conv_ext_efield;
  double m_conv_energy;
  double m_conv_density;
  double m_conv_flux;
  double m_tau;
  double m_density;

public:
  EKSpecies() {
    add_parameters(
        {{"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
         {"single_precision", AutoParameter::read_only,
          [this]() { return not m_instance->is_double_precision(); }},
         {"gpu", AutoParameter::read_only,
          [this]() { return m_instance->is_gpu(); }},
         {"diffusion",
          [this](Variant const &v) {
            m_instance->set_diffusion(get_value<double>(v) * m_conv_diffusion);
          },
          [this]() { return m_instance->get_diffusion() / m_conv_diffusion; }},
         {"kT",
          [this](Variant const &v) {
            context()->parallel_try_catch([&]() {
              auto const kT = get_value<double>(v);
              if (kT < 0.) {
                throw std::domain_error("Parameter 'kT' must be >= 0");
              }
              m_instance->set_kT(kT * m_conv_energy);
            });
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
          [this]() {
            return m_instance->get_ext_efield() / m_conv_ext_efield;
          }},
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
         {"tau", AutoParameter::read_only, [this]() { return m_tau; }},
         {"density", AutoParameter::read_only,
          [this]() { return m_density / m_conv_density; }},
         {"thermalized", AutoParameter::read_only,
          [this]() { return m_instance->is_thermalized(); }},
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
         {"shape", AutoParameter::read_only,
          [this]() { return m_instance->get_lattice().get_grid_dimensions(); }},
         {"vtk_writers", AutoParameter::read_only,
          [this]() { return serialize_vtk_writers(); }}});
  }

  void do_construct(VariantMap const &params) override;

  [[nodiscard]] auto get_ekinstance() const { return m_instance; }
  [[nodiscard]] auto get_lattice() const { return m_lattice; }
  [[nodiscard]] auto get_mpi_cart_comm_observer() const {
    return m_mpi_cart_comm_observer;
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override;

  [[nodiscard]] auto get_conversion_factor_density() const noexcept {
    return m_conv_density;
  }
  [[nodiscard]] auto get_conversion_factor_flux() const noexcept {
    return m_conv_flux;
  }

  ::LatticeModel::units_map
  get_lattice_to_md_units_conversion() const override {
    return {
        {"density", 1. / m_conv_density},
        {"flux", 1. / m_conv_flux},
    };
  }

  void flux_boundary_ghost_layer_size_sanity_check() const {
    context()->parallel_try_catch([&]() {
      if (get_lattice()->lattice()->get_ghost_layers() < 2 and
          context()->get_comm().size() > 1) {
        throw std::runtime_error("The number of ghostlayers should be > 1 "
                                 "when using flux boundaries and MPI");
      }
    });
  }

protected:
  void make_instance(VariantMap const &params) override;

private:
  void load_checkpoint(std::filesystem::path const &path, int mode);
  void save_checkpoint(std::filesystem::path const &path, int mode);
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
