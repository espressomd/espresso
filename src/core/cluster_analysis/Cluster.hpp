/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CLUSTER_ANALYSIS_CLUSTER_HPP
#define CLUSTER_ANALYSIS_CLUSTER_HPP

#include <utils/Vector.hpp>
#include <vector>

#include "Particle.hpp"

namespace ClusterAnalysis {

/** @brief Represents a single cluster of particles */
class Cluster {
public:
  /** @brief Ids of the particles in the cluster */
  std::vector<int> particles;
  /** @brief add a particle to the cluster */
  void add_particle(const Particle &p) { particles.push_back(p.p.identity); }
  /** @brief Calculate the center of mass of the cluster */
  Utils::Vector3d
  center_of_mass_subcluster(std::vector<int> &subcl_partcicle_ids);
  Utils::Vector3d center_of_mass();
  /** @brief Longest distance between any combination of two particles */
  double longest_distance();
  /** @brief Calculate radius of gyration of the cluster */
  double radius_of_gyration();
  double radius_of_gyration_subcluster(std::vector<int> &subcl_particle_ids);
  /** @brief Calculate the fractal dimension
   *  N(r) via r^d, where N(r) counts the number of particles in a sphere
   *  of radius n, and d denotes the fractal dimension.
   *  The fitting is done by the Gnu Scientific Library.
   *  @param dr   increment for when constructing the discrete version of N(r)
   *
   *  @return fractal dimension, rms error of the fit */
  std::pair<double, double> fractal_dimension(double dr);
};

} // namespace ClusterAnalysis

#endif
