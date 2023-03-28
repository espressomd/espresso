/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#ifdef DIPOLAR_BARNES_HUT

#include "magnetostatics/barnes_hut_gpu.hpp"
#include "magnetostatics/barnes_hut_gpu_cuda.cuh"

#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"

DipolarBarnesHutGpu::DipolarBarnesHutGpu(double prefactor, double epssq,
                                         double itolsq)
    : prefactor{prefactor}, m_epssq{epssq}, m_itolsq{itolsq} {
  if (prefactor <= 0.) {
    throw std::domain_error("Parameter 'prefactor' must be > 0");
  }
  if (m_itolsq <= 0.) {
    throw std::domain_error("Parameter 'itolsq' must be > 0");
  }
  if (m_epssq <= 0.) {
    throw std::domain_error("Parameter 'epssq' must be > 0");
  }
  auto &system = EspressoSystemInterface::Instance();
  system.requestFGpu();
  system.requestTorqueGpu();
  system.requestRGpu();
  system.requestDipGpu();
  if (this_node == 0) {
    setBHPrecision(static_cast<float>(m_epssq), static_cast<float>(m_itolsq));
  }
}

template <class... Args, class... ArgRef>
int call_kernel(void (*fp)(Args...), ArgRef &&...args) {
  int error_code = ES_ERROR;
  try {
    fp(args...);
    error_code = ES_OK;
  } catch (std::runtime_error const &err) {
    runtimeErrorMsg() << "DipolarBarnesHutGpu: " << err.what();
  }
  return error_code;
}

int DipolarBarnesHutGpu::initialize_data_structure() {
  auto &system = EspressoSystemInterface::Instance();
  auto const n_part = static_cast<int>(system.npart_gpu());
  auto const error_code = call_kernel(allocBHmemCopy, n_part, &m_bh_data);

  if (error_code == ES_OK) {
    fill_bh_data(system.rGpuBegin(), system.dipGpuBegin(), &m_bh_data);
    initBHgpu(m_bh_data.blocks);
    buildBoxBH(m_bh_data.blocks);
    buildTreeBH(m_bh_data.blocks);
    summarizeBH(m_bh_data.blocks);
    sortBH(m_bh_data.blocks);
  }

  return error_code;
}

void DipolarBarnesHutGpu::add_long_range_forces() {
  auto &system = EspressoSystemInterface::Instance();
  system.update();
  if (this_node == 0) {
    if (initialize_data_structure() == ES_OK) {
      call_kernel(forceBH, &m_bh_data, static_cast<float>(prefactor),
                  system.fGpuBegin(), system.torqueGpuBegin());
    }
  }
}

void DipolarBarnesHutGpu::long_range_energy() {
  auto &system = EspressoSystemInterface::Instance();
  system.update();
  if (this_node == 0) {
    if (initialize_data_structure() == ES_OK) {
      auto energy = &(reinterpret_cast<CUDA_energy *>(system.eGpu())->dipolar);
      call_kernel(energyBH, &m_bh_data, static_cast<float>(prefactor), energy);
    }
  }
}

#endif // DIPOLAR_BARNES_HUT
