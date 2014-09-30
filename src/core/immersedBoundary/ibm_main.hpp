/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation


#ifndef IBM_MAIN_HPP
#define IBM_MAIN_HPP

#ifdef IMMERSED_BOUNDARY

// Dummy functions. They are required by Espresso, but they don't do anything here.
// We have our own update functions.
void update_mol_pos_particle(Particle *);
void update_mol_vel_particle(Particle *);
void distribute_mol_force();

// Main functions for CPU implementation
void IBM_ForcesIntoFluid();
void IBM_UpdateParticlePositions();
void IBM_ResetLBForces();

#endif

#endif