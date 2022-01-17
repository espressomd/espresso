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

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "DipolarBarnesHut.hpp"
#include "DipolarBarnesHut_cuda.cuh"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include "SystemInterface.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "energy.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"

DipolarBarnesHut::DipolarBarnesHut(SystemInterface &s) {
  s.requestFGpu();
  s.requestTorqueGpu();
  s.requestRGpu();
  s.requestDipGpu();
}

void DipolarBarnesHut::activate() {
  dipole.method = DIPOLAR_BH_GPU;
  mpi_bcast_coulomb_params();

  forceActors.push_back(this);
  energyActors.push_back(this);
}

void DipolarBarnesHut::deactivate() {
  dipole.method = DIPOLAR_NONE;
  mpi_bcast_coulomb_params();

  forceActors.remove(this);
  energyActors.remove(this);
}

void DipolarBarnesHut::set_params(float epssq, float itolsq) {
  m_prefactor = static_cast<float>(dipole.prefactor);
  m_epssq = epssq;
  m_itolsq = itolsq;
  setBHPrecision(m_epssq, m_itolsq);
}

int DipolarBarnesHut::initialize_data_structure(SystemInterface &s) {
  try {
    allocBHmemCopy(static_cast<int>(s.npart_gpu()), &m_bh_data);
  } catch (cuda_runtime_error const &err) {
    runtimeErrorMsg() << "DipolarBarnesHut: " << err.what();
    return ES_ERROR;
  }

  fill_bh_data(s.rGpuBegin(), s.dipGpuBegin(), &m_bh_data);
  initBHgpu(m_bh_data.blocks);
  buildBoxBH(m_bh_data.blocks);
  buildTreeBH(m_bh_data.blocks);
  summarizeBH(m_bh_data.blocks);
  sortBH(m_bh_data.blocks);

  return ES_OK;
}

void DipolarBarnesHut::computeForces(SystemInterface &s) {
  if (initialize_data_structure(s) == ES_OK) {
    if (forceBH(&m_bh_data, m_prefactor, s.fGpuBegin(), s.torqueGpuBegin())) {
      runtimeErrorMsg() << "kernels encountered a functional error";
    }
  }
}

void DipolarBarnesHut::computeEnergy(SystemInterface &s) {
  if (initialize_data_structure(s) == ES_OK) {
    auto energy = reinterpret_cast<CUDA_energy *>(s.eGpu());
    if (energyBH(&m_bh_data, m_prefactor, &(energy->dipolar))) {
      runtimeErrorMsg() << "kernels encountered a functional error";
    }
  }
}

#endif
