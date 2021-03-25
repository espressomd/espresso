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

#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS
#include "ParticleRange.hpp"

#include <boost/mpi/communicator.hpp>

#include <unordered_map>

void set_sd_viscosity(double eta);

void set_sd_radius_dict(std::unordered_map<int, double> const &x);

void set_sd_kT(double kT);
double get_sd_kT();

void set_sd_flags(int flg);

/** Takes the forces and torques on all particles and computes their
 *  velocities. Acts globally on particles on all nodes; i.e. particle data
 *  is gathered from all nodes and their velocities and angular velocities are
 *  set according to the Stokesian Dynamics method.
 */
void propagate_vel_pos_sd(const ParticleRange &particles,
                          const boost::mpi::communicator &comm,
                          double time_step);

#endif // STOKESIAN_DYNAMICS

#endif // STOKESIAN_DYNAMICS_INTERFACE_H
