/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// *******
// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main Espresso routines
// Functions to be exported for Espresso are in ibm_main.hpp

#ifndef IBM_CUDA_INTERFACE_HPP
#define IBM_CUDA_INTERFACE_HPP

#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "ParticleRange.hpp"

// *********** Communication functions ********
// Implemented in real C++, but called from the ibm_cuda.cu
void IBM_cuda_mpi_send_velocities(ParticleRange particles);
void IBM_cuda_mpi_get_particles(ParticleRange particles);

void ParticleVelocitiesFromLB_GPU(ParticleRange particles);

// ******** data types for CUDA and MPI communication ******
typedef struct {
  float pos[3];
  float f[3];
  bool is_virtual;
} IBM_CUDA_ParticleDataInput;

typedef struct {
  float v[3];
} IBM_CUDA_ParticleDataOutput;

// ******** global variables for CUDA and MPI communication ******
extern IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_host;
extern IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_host;

#endif

#endif
