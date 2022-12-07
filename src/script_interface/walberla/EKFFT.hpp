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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_FFT_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_FFT_HPP

#include "config/config.hpp"

#ifdef WALBERLA
#ifdef WALBERLA_FFT

#include "EKPoissonSolver.hpp"

#include "walberla_bridge/electrokinetics/ek_walberla_init.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <memory>

namespace ScriptInterface::walberla {

class EKFFT : public EKPoissonSolver {
public:
  void do_construct(VariantMap const &args) override {
    m_single_precision = get_value_or<bool>(args, "single_precision", false);
    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();
    auto permittivity = get_value<double>(args, "permittivity");

    m_instance = new_ek_poisson_fft(lattice, permittivity, m_single_precision);

    add_parameters({{"permittivity",
                     [this](Variant const &v) {
                       m_instance->set_permittivity(get_value<double>(v));
                     },
                     [this]() { return m_instance->get_permittivity(); }},
                    {"is_single_precision", AutoParameter::read_only,
                     [this]() { return m_single_precision; }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver>
  get_instance() const noexcept override {
    return m_instance;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::PoissonSolver> m_instance;

  bool m_single_precision;
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA_FFT
#endif // WALBERLA
#endif
