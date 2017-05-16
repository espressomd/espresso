/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation


#ifndef IBM_MAIN_HPP
#define IBM_MAIN_HPP

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "ParticleRange.hpp"

// Dummy functions. They are required by Espresso, but they don't do anything here.
// We have our own update functions.
void update_mol_pos_particle(Particle *);
void update_mol_vel_particle(Particle *);
void distribute_mol_force();

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
