#ifndef _LBTRACERS_H
#define _LBTRACERS_H

#include "config.hpp"
#include "particle_data.hpp"

#ifdef LBTRACERS

//Update Position ~ Euler/Runge-Kutta
void update_mol_pos_particle(Particle *);
//Update Velocity ~ Get interpolated velocity of LB
void update_mol_vel_particle(Particle *);
//Since no 'real' particles are involved, this function will stay empty
void distribute_mol_force();

#endif

#endif
