#ifndef STOKESIAN_DYNAMICS_INTERFACE_H
#define STOKESIAN_DYNAMICS_INTERFACE_H

#include "config.hpp"

#include <string>
#include <unordered_map>

#ifdef STOKESIAN_DYNAMICS

void set_sd_viscosity(double eta);
double get_sd_viscosity();

void set_sd_device(std::string const &dev);
std::string get_sd_device();

void set_sd_radius_dict(std::unordered_map<int, double> const &x);
std::unordered_map<int, double> get_sd_radius_dict();

void set_sd_kT(double kT);
double get_sd_kT();

void set_sd_flags(int flg);
int get_sd_flags();

void set_sd_seed(std::size_t seed);
std::size_t get_sd_seed();

void propagate_vel_pos_sd();

#endif // STOKESIAN_DYNAMICS

#endif // STOKESIAN_DYNAMICS_INTERFACE_H
