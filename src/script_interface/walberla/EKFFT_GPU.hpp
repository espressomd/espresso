/*
 * Copyright (C) 2025 The ESPResSo project
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
#ifdef ESPRESSO_WALBERLA_FFT

#include "EKPoissonSolver.hpp"
#include "LatticeWalberla.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/communication.hpp"

#include <walberla_bridge/electrokinetics/ek_poisson_fft_gpu_init.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>
#include <script_interface/walberla/EKFFT.hpp>

#include <utils/math/int_pow.hpp>

namespace ScriptInterface::walberla {

class EKFFTGPU : public EKFFT {

public:
  void make_instance(VariantMap const &args) override {
    // unit conversions
    auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
    m_conv_permittivity = Utils::int_pow<2>(agrid);
    auto const permittivity =
        get_value<double>(args, "permittivity") * m_conv_permittivity;

    m_instance = ::walberla::new_ek_poisson_fft_cuda(
        m_lattice->lattice(), permittivity, m_single_precision);
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA_FFT
#endif // ESPRESSO_WALBERLA
