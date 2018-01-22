/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation


#ifndef IBM_MAIN_HPP
#define IBM_MAIN_HPP

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "ParticleRange.hpp"

// Main functions for CPU & GPU
void IBM_UpdateParticlePositions(ParticleRange particles);

// Main functions for CPU implementation - called from integrate.cpp
void IBM_ForcesIntoFluid_CPU();
void IBM_ResetLBForces_CPU();

// Main functions for GPU implementation - called from integrate.cpp
// These are in ibm_cuda.cu
void IBM_ForcesIntoFluid_GPU(ParticleRange particles);
void IBM_ResetLBForces_GPU();

#endif

#endif
