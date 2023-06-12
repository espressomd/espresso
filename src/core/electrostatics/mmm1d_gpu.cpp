/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef MMM1D_GPU

#include "electrostatics/mmm1d_gpu.hpp"

#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"

#include <stdexcept>

CoulombMMM1DGpu::CoulombMMM1DGpu(double prefactor, double maxPWerror,
                                 double far_switch_radius, int bessel_cutoff)
    : maxPWerror{maxPWerror}, far_switch_radius{far_switch_radius},
      far_switch_radius_sq{-1.}, bessel_cutoff{bessel_cutoff}, m_is_tuned{
                                                                   false} {
  set_prefactor(prefactor);
  if (maxPWerror <= 0.) {
    throw std::domain_error("Parameter 'maxPWerror' must be > 0");
  }
  if (far_switch_radius <= 0. and far_switch_radius != -1.) {
    throw std::domain_error("Parameter 'far_switch_radius' must be > 0");
  }
  if (far_switch_radius > 0. and far_switch_radius > box_geo.length()[2]) {
    throw std::domain_error(
        "Parameter 'far_switch_radius' must not be larger than box length");
  }
  if (bessel_cutoff < 0 and bessel_cutoff != -1) {
    throw std::domain_error("Parameter 'bessel_cutoff' must be > 0");
  }

  auto &gpu_particle_data = System::get_system().gpu;
  gpu_particle_data.enable_property(GpuParticleData::prop::force);
  gpu_particle_data.enable_property(GpuParticleData::prop::pos);
  gpu_particle_data.enable_property(GpuParticleData::prop::q);
  if (this_node == 0) {
    modpsi_init();
  }
}

void CoulombMMM1DGpu::sanity_checks_periodicity() const {
  if (box_geo.periodic(0) || box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("MMM1D requires periodicity (False, False, True)");
  }
}

void CoulombMMM1DGpu::sanity_checks_cell_structure() const {
  if (local_geo.cell_structure_type() !=
      CellStructureType::CELL_STRUCTURE_NSQUARE) {
    throw std::runtime_error("MMM1D requires the N-square cellsystem");
  }
}

void CoulombMMM1DGpu::tune() {
  System::get_system().gpu.update();
  if (this_node == 0) {
    setup();
    tune(maxPWerror, far_switch_radius, bessel_cutoff);
  }
  m_is_tuned = true;
}

#endif
