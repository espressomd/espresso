/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifdef DIPOLAR_DIRECT_SUM

#include "DipolarDirectSum.hpp"
#include "DipolarDirectSum_cuda.cuh"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include "SystemInterface.hpp"
#include "cuda_interface.hpp"
#include "energy.hpp"
#include "forces.hpp"
#include "grid.hpp"

DipolarDirectSum::DipolarDirectSum(SystemInterface &s) {
  s.requestFGpu();
  s.requestTorqueGpu();
  s.requestRGpu();
  s.requestDipGpu();
}

void DipolarDirectSum::activate() {
  // also necessary on 1 CPU or GPU, does more than just broadcasting
  dipole.method = DIPOLAR_DS_GPU;
  mpi_bcast_coulomb_params();

  forceActors.push_back(this);
  energyActors.push_back(this);
}

void DipolarDirectSum::deactivate() {
  dipole.method = DIPOLAR_NONE;
  mpi_bcast_coulomb_params();

  forceActors.remove(this);
  energyActors.remove(this);
}

void DipolarDirectSum::set_params() {
  m_prefactor = static_cast<float>(dipole.prefactor);
}

void DipolarDirectSum::computeForces(SystemInterface &s) {
  float box[3];
  int per[3];
  for (int i = 0; i < 3; i++) {
    box[i] = static_cast<float>(s.box()[i]);
    per[i] = box_geo.periodic(i);
  }
  DipolarDirectSum_kernel_wrapper_force(
      m_prefactor, static_cast<unsigned>(s.npart_gpu()), s.rGpuBegin(),
      s.dipGpuBegin(), s.fGpuBegin(), s.torqueGpuBegin(), box, per);
}

void DipolarDirectSum::computeEnergy(SystemInterface &s) {
  float box[3];
  int per[3];
  for (int i = 0; i < 3; i++) {
    box[i] = static_cast<float>(s.box()[i]);
    per[i] = box_geo.periodic(i);
  }
  auto energy = reinterpret_cast<CUDA_energy *>(s.eGpu());
  DipolarDirectSum_kernel_wrapper_energy(
      m_prefactor, static_cast<unsigned>(s.npart_gpu()), s.rGpuBegin(),
      s.dipGpuBegin(), box, per, &(energy->dipolar));
}

#endif
