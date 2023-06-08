/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_P3M_GPU_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_P3M_GPU_HPP

/**
 * @file
 * P3M electrostatics on GPU.
 */

#include "config/config.hpp"

#ifdef P3M
#ifdef CUDA

#include "electrostatics/p3m.hpp"

#include "ParticleRange.hpp"

struct CoulombP3MGPU : public CoulombP3M {
  using CoulombP3M::CoulombP3M;

  void init();
  void on_activation() {
    request_gpu();
    CoulombP3M::on_activation();
    if (is_tuned()) {
      init_cpu_kernels();
    }
  }

  void add_long_range_forces(ParticleRange const &);
  void request_gpu() const;

private:
  /**
   * @brief Initialize the CPU kernels.
   * This operation is time-consuming and sets up data members
   * that are only relevant for ELC force corrections.
   */
  void init_cpu_kernels();
};

#endif // CUDA
#endif // P3M

#endif
