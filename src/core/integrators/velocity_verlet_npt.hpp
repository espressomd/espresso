#ifndef INTEGRATORS_VELOCITY_VERLET_NPT_HPP
#define INTEGRATORS_VELOCITY_VERLET_NPT_HPP

#include "config.hpp"

#ifdef NPT
#include "ParticleRange.hpp"
#include "particle_data.hpp"

void velocity_verlet_npt_step_1(const ParticleRange &particles);

void velocity_verlet_npt_step_2(const ParticleRange &particles);
#endif

#endif
