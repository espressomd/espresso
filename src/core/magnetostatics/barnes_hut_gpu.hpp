/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_BARNES_HUT_GPU_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_BARNES_HUT_GPU_HPP

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "magnetostatics/barnes_hut_gpu_cuda.cuh"

struct DipolarBarnesHutGpu {
  double prefactor;
  double m_epssq;
  double m_itolsq;
  DipolarBarnesHutGpu(double prefactor, double epssq, double itolsq);
  ~DipolarBarnesHutGpu() { deallocBH(&m_bh_data); }

  void on_activation() const {}
  void on_boxl_change() const {}
  void on_node_grid_change() const {}
  void on_periodicity_change() const {}
  void on_cell_structure_change() const {}
  void init() const {}
  void sanity_checks() const {}

  void add_long_range_forces();
  void long_range_energy();

private:
  /// Container for pointers to device memory.
  BHData m_bh_data = {0,       0,       0,       nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr};
  int initialize_data_structure();
};

#endif // DIPOLAR_BARNES_HUT

#endif
