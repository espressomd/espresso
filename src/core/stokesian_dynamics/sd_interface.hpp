#ifndef STOKESIAN_DYNAMICS_INTERFACE_H
#define STOKESIAN_DYNAMICS_INTERFACE_H

#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS

void set_sd_viscosity(double eta);
double get_sd_viscosity();
void set_sd_device(std::string const &dev);
std::string get_sd_device();
void propagate_vel_pos_sd();

#endif // STOKESIAN_DYNAMICS

#endif // STOKESIAN_DYNAMICS_INTERFACE_H
