/*
 * Copyright (C) 2010-2020 The ESPResSo project
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

/** @file
 *  See @cite durlofsky87a for the Stokesian dynamics method used here.
 *  See @cite banchio03a and @cite brady88a for the thermalization method.
 */

#ifndef STOKESIAN_DYNAMICS_INTERFACE_H
#define STOKESIAN_DYNAMICS_INTERFACE_H

#include "ParticleRange.hpp"
#include "config.hpp"

#include <string>
#include <unordered_map>
#include <utils/Vector.hpp>

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

/** Takes the forces and torques on all particles and computes their
 *  velocities. Acts globally on particles on all nodes; i.e. particle data
 *  is gathered from all nodes and their velocities and angular velocities are
 *  set according to the Stokesian Dynamics method.
 */
void propagate_vel_pos_sd(const ParticleRange &particles);

#endif // STOKESIAN_DYNAMICS

#endif // STOKESIAN_DYNAMICS_INTERFACE_H
