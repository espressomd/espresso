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
/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation

#ifndef VIRTUAL_SITES_LB_INERTIALESS_TRACERS_HPP
#define VIRTUAL_SITES_LB_INERTIALESS_TRACERS_HPP

#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "ParticleRange.hpp"

// Main functions for CPU & GPU
void IBM_UpdateParticlePositions(ParticleRange particles);

// Main functions for CPU implementation - called from integrate.cpp
void IBM_ForcesIntoFluid_CPU();

// Main functions for GPU implementation - called from integrate.cpp
// These are in ibm_cuda.cu
void IBM_ForcesIntoFluid_GPU(ParticleRange particles);
void IBM_ResetLBForces_GPU();

#endif

#endif
