/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef _STATISTICS_H
#define _STATISTICS_H
/** \file
 *  Statistical tools to analyze simulations.
 *
 *  Implementation in statistics.cpp.
 */

#include "PartCfg.hpp"

#include <vector>

/** Calculate the minimal distance of two particles with types in set1 resp.
 *  set2.
 *  @param partCfg particle collection.
 *  @param set1 types of particles
 *  @param set2 types of particles
 *  @return the minimal distance of two particles
 */
double mindist(PartCfg &partCfg, const std::vector<int> &set1,
               const std::vector<int> &set2);

/** Find all particles within a given radius @p r_catch around a position.
 *  @param partCfg    @copybrief PartCfg
 *  @param pos        position of sphere center
 *  @param r_catch    the sphere radius
 *  @param planedims  orientation of coordinate system
 *
 *  @return List of ids close to @p pos.
 */
std::vector<int> nbhood(PartCfg &partCfg, const Utils::Vector3d &pos,
                        double r_catch, const Utils::Vector3i &planedims);

/** Calculate minimal distance to point.
 *  @param partCfg particle selection
 *  @param pos  point
 *  @param pid  if a valid particle id, this particle is omitted from
 *              minimization (this is a good idea if @p pos is the
 *              position of a particle).
 *  @return the minimal distance of a particle to coordinates @p pos
 */
double distto(PartCfg &partCfg, const Utils::Vector3d &pos, int pid = -1);

/** Calculate the distribution of particles around others.
 *
 *  Calculates the distance distribution of particles with types given
 *  in the @p p1_types list around particles with types given in the
 *  @p p2_types list. The distances range from @p r_min to @p r_max, binned
 *  into @p r_bins bins which are either equidistant (@p log_flag==false) or
 *  logarithmically equidistant (@p log_flag==true). The result is stored
 *  in the @p array dist.
 *  @param partCfg  particle collection.
 *  @param p1_types list with types of particles to find the distribution for.
 *  @param p2_types list with types of particles the others are distributed
 *                  around.
 *  @param r_min    Minimal distance for the distribution.
 *  @param r_max    Maximal distance for the distribution.
 *  @param r_bins   Number of bins.
 *  @param log_flag Whether the bins are (logarithmically) equidistant.
 *  @param low      particles closer than @p r_min
 *  @param dist     Array to store the result (size: @p r_bins).
 */
void calc_part_distribution(PartCfg &partCfg, std::vector<int> const &p1_types,
                            std::vector<int> const &p2_types, double r_min,
                            double r_max, int r_bins, bool log_flag,
                            double *low, double *dist);

/** Calculate the spherically averaged structure factor.
 *
 *  Calculates the spherically averaged structure factor of particles of a
 *  given type. The possible wave vectors are given by q = 2PI/L sqrt(nx^2 +
 *  ny^2 + nz^2).
 *  The S(q) is calculated up to a given length measured in 2PI/L (the
 *  recommended order of the wave vector is less than 20).
 *  The data is stored starting with q=1, and contains alternatingly S(q-1) and
 *  the number of wave vectors l with l^2=q. Only if the second number is
 *  nonzero, the first is meaningful. This means the q=1 entries are sf[0]=S(1)
 *  and sf[1]=1. For q=7, there are no possible wave vectors, so
 *  sf[2*(7-1)]=sf[2*(7-1)+1]=0.
 *
 *  @param partCfg   particle collection
 *  @param p_types   list with types of particles to be analyzed
 *  @param order     the maximum wave vector length in 2PI/L
 */
std::vector<double> calc_structurefactor(PartCfg &partCfg,
                                         std::vector<int> const &p_types,
                                         int order);

std::vector<std::vector<double>> modify_stucturefactor(int order,
                                                       double const *sf);

/** Calculate the center of mass of a special type of the current configuration.
 *  \param part_type  type of the particle
 */
Utils::Vector3d centerofmass(PartCfg &, int part_type);

/** Calculate the angular momentum of a special type of the current
 *  configuration.
 *  \param type  type of the particle
 */
Utils::Vector3d angularmomentum(PartCfg &, int type);

/** Calculate the center of mass of a special type of a saved configuration.
 *  \param partCfg     @copybrief PartCfg
 *  \param type        type of the particle, -1 for all
 *  \param MofImatrix  Center of mass
 */
void momentofinertiamatrix(PartCfg &partCfg, int type, double *MofImatrix);

/** Calculate total momentum of the system (particles & LB fluid).
 *  Inputs are bools to include particles and fluid in the linear momentum
 *  calculation
 *  @return Result for this processor
 */
Utils::Vector3d calc_linear_momentum(int include_particles,
                                     int include_lbfluid);

#endif
