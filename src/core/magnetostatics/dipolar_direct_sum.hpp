/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifdef ESPRESSO_DIPOLES

#include "magnetostatics/actor.hpp"

#include <utils/Vector.hpp>

/**
 * @brief Dipolar all with all and no replica.
 * Handling of a system of dipoles where no replicas exist.
 * Assumes minimum image convention for those axis in which the
 * system is periodic.
 */
struct DipolarDirectSum : public Dipoles::Actor<DipolarDirectSum> {
  int n_replicas;
  bool m_is_gpu;
  DipolarDirectSum(double prefactor, int n_replicas, bool gpu);

  bool is_gpu() const { return m_is_gpu; }
  void on_activation() const {
#ifdef ESPRESSO_CUDA
    if (m_is_gpu) {
      on_activation_gpu();
    }
#endif
  }
#ifdef ESPRESSO_CUDA
  void on_activation_gpu() const;
#endif
  void on_boxl_change() const {}
  void on_node_grid_change() const {}
  void on_periodicity_change() const {}
  void on_cell_structure_change() const {}
  void init() const {}
  void sanity_checks() const {}

  /**
   * @brief Calculate long-range dipolar energy.
   * The GPU implementation stores the energy on a GPU accumulator
   * and zero is returned from this method.
   */
  double long_range_energy() const {
#ifdef ESPRESSO_CUDA
    if (m_is_gpu) {
      long_range_energy_gpu();
      return 0.;
    }
#endif
    return long_range_energy_cpu();
  }
  void add_long_range_forces() const {
#ifdef ESPRESSO_CUDA
    if (m_is_gpu) {
      return add_long_range_forces_gpu();
    }
#endif
    add_long_range_forces_cpu();
  }
  double long_range_energy_cpu() const;
  void add_long_range_forces_cpu() const;
#ifdef ESPRESSO_CUDA
  void long_range_energy_gpu() const;
  void add_long_range_forces_gpu() const;
#endif

  /**
   * @brief Calculate the dipolar pressure tensor.
   * Not implemented for the GPU variant: emits a runtime warning
   * and returns a zero tensor in that case.
   */
  Utils::Vector9d long_range_pressure() const;
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  void dipole_field_at_part_cpu() const;
#endif
};

#endif // ESPRESSO_DIPOLES
