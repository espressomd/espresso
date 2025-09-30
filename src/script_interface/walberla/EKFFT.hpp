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

#ifdef ESPRESSO_WALBERLA_FFT

#include "EKPoissonSolver.hpp"

#include "LatticeWalberla.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/communication.hpp"

#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <utils/math/int_pow.hpp>

#include <memory>

namespace ScriptInterface::walberla {

class EKFFT : public EKPoissonSolver {
protected:
  std::unique_ptr<ResourceManager> m_resources_lock;
  std::shared_ptr<LatticeWalberla> m_lattice;
  double m_conv_permittivity;
  bool m_single_precision;

public:
  void make_instance(VariantMap const &args) override {
    // unit conversions
    auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
    m_conv_permittivity = Utils::int_pow<2>(agrid);
    auto const permittivity =
        get_value<double>(args, "permittivity") * m_conv_permittivity;

    m_instance = ::walberla::new_ek_poisson_fft(
        m_lattice->lattice(), permittivity, m_single_precision);
    m_instance->setup_fft(false);
  }

  void do_construct(VariantMap const &args) override {
    m_single_precision = get_value_or<bool>(args, "single_precision", false);
    m_lattice = get_value<decltype(m_lattice)>(args, "lattice");
    m_vtk_writers =
        get_value_or<decltype(m_vtk_writers)>(args, "vtk_writers", {});

    make_instance(args), m_resources_lock = std::make_unique<ResourceManager>();
    // MPI communicator is needed to destroy the FFT plans
    m_resources_lock->acquire_lock(::communication_environment->get_mpi_env());
    for (auto &vtk : m_vtk_writers) {
      vtk->attach_to_lattice(m_instance, get_lattice_to_md_units_conversion());
    }
  }

  EKFFT() {
    add_parameters({
        {"permittivity",
         [this](Variant const &v) {
           m_instance->set_permittivity(get_value<double>(v) *
                                        m_conv_permittivity);
         },
         [this]() {
           return m_instance->get_permittivity() / m_conv_permittivity;
         }},
        {"single_precision", AutoParameter::read_only,
         [this]() { return m_single_precision; }},
        {"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_instance->get_lattice().get_grid_dimensions(); }},
    });
  }

  ~EKFFT() override {
    m_lattice.reset();
    m_instance.reset();
    m_resources_lock.reset();
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver>
  get_instance() const noexcept override {
    return m_instance;
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA_FFT
