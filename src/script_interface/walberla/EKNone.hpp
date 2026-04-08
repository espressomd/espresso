/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "EKPoissonSolver.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>
#include <script_interface/code_info/CodeInfo.hpp>

#include <memory>

namespace ScriptInterface::walberla {

class EKNone : public EKPoissonSolver {
  std::shared_ptr<::walberla::PoissonSolver> m_instance;
  std::shared_ptr<LatticeWalberla> m_lattice;
  bool m_gpu;
  bool m_single_precision;

protected:
  void make_instance(VariantMap const &args) override {

    auto *make_new_instance = &::walberla::new_ek_poisson_none;
    if (m_gpu) {
      std::vector<std::string> required_features;
      required_features.emplace_back("CUDA");
      CodeInfo::check_features(required_features);
#ifdef ESPRESSO_CUDA
      make_new_instance = &::walberla::new_ek_poisson_none_cuda;
#endif
    }
    m_instance = make_new_instance(m_lattice->lattice(), m_single_precision);
  }

public:
  void do_construct(VariantMap const &args) override {
    m_gpu = get_value_or<bool>(args, "gpu", false);
    m_single_precision = get_value_or<bool>(args, "single_precision", m_gpu);
    m_lattice = get_value<decltype(m_lattice)>(args, "lattice");

    make_instance(args);
    add_parameters({
        {"single_precision", AutoParameter::read_only,
         [this]() { return m_single_precision; }},
        {"gpu", AutoParameter::read_only, [this]() { return m_gpu; }},
        {"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
    });
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver>
  get_instance() const noexcept override {
    return m_instance;
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
