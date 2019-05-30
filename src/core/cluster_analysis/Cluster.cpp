/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#ifdef GSL
#include "gsl/gsl_fit.h"
#endif
#include <utils/Vector.hpp>
#include <vector>

#include "Cluster.hpp"

namespace ClusterAnalysis {

// Center of mass of an aggregate
Utils::Vector3d Cluster::center_of_mass() {
  return center_of_mass_subcluster(particles);
}

// Center of mass of an aggregate
Utils::Vector3d
Cluster::center_of_mass_subcluster(std::vector<int> &subcl_partcicle_ids) {
  Utils::Vector3d com{};

  // The distances between the particles "folded", such that all distances
  // are smaller than box_l/2 in a periodic system. The 1st particle
  // of the cluster is arbitrarily chosen as reference.

  Utils::Vector3d reference_position = folded_position(partCfg()[particles[0]]);
  Utils::Vector3d dist_to_reference;
  double total_mass = 0.;
  for (int pid :
       subcl_partcicle_ids) // iterate over all particle ids within a cluster
  {
    const Utils::Vector3d folded_pos = folded_position(partCfg()[pid]);
    get_mi_vector(dist_to_reference, folded_pos,
                  reference_position); // add current particle positions
    com = com + dist_to_reference * partCfg()[pid].p.mass;
    total_mass += partCfg()[pid].p.mass;
  }

  // Normalize by number of particles
  com = com * 1. / total_mass;

  // Re-add reference position
  com = com + reference_position;

  // Fold into simulation box

  for (int i = 0; i < 3; i++) {
    com[i] = fmod(com[i], box_l[i]);
  }
  return com;
}

double Cluster::longest_distance() {
  double ld = 0.;
  for (auto a = particles.begin(); a != particles.end(); a++) {
    for (auto b = a; ++b != particles.end();) {
      auto const dist =
          get_mi_vector(partCfg()[*a].r.p, partCfg()[*b].r.p).norm();

      // Larger than previous largest distance?
      ld = std::max(ld, dist);
    }
  }
  return ld;
}

// Radius of gyration
double Cluster::radius_of_gyration() {
  return radius_of_gyration_subcluster(particles);
}

double
Cluster::radius_of_gyration_subcluster(std::vector<int> &subcl_particle_ids) {
  // Center of mass
  Utils::Vector3d com = center_of_mass_subcluster(subcl_particle_ids);
  double sum_sq_dist = 0.;
  for (auto const pid : subcl_particle_ids) {
    // calculate square length of this distance
    sum_sq_dist += get_mi_vector(com, partCfg()[pid].r.p).norm2();
  }

  return sqrt(sum_sq_dist / subcl_particle_ids.size());
}

template <typename T>
std::vector<std::size_t> sort_indices(const std::vector<T> &v) {

  // Unsorted for unsorted vector (0..n-1)
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
  return idx;
}

std::pair<double, double> Cluster::fractal_dimension(double dr) {
#ifdef GSL
  Utils::Vector3d com = center_of_mass();
  // calculate Df using linear regression on the logarithms of the radii of
  // gyration against the number of particles in sub-clusters. Particles are
  // included step by step from the center of mass outwards

  // Distances of particles from the center of mass
  std::vector<double> distances;

  for (auto const &it : particles) {
    distances.push_back(get_mi_vector(com.begin(), partCfg()[it].r.p)
                            .norm()); // add distance from the current particle
                                      // to the com in the distances vectors
  }

  // Get particle indices in the cluster which yield distances  sorted in
  // ascending order from center of mass.
  auto particle_indices = sort_indices(distances);

  // Particle ids in the current sub-cluster
  std::vector<int> subcluster_ids;

  std::vector<double> log_pcounts;   // particle count
  std::vector<double> log_diameters; // corresponding radii of gyration
  double last_dist = 0;
  for (auto const idx : particle_indices) {
    subcluster_ids.push_back(particles[idx]);
    if (distances[idx] < last_dist + dr)
      continue;

    last_dist = distances[idx];
    if (subcluster_ids.size() == 1)
      continue;
    double current_rg = radius_of_gyration_subcluster(subcluster_ids);
    log_pcounts.push_back(log(subcluster_ids.size()));
    log_diameters.push_back(log(current_rg * 2.0));
  }
  // usage: Function: int gsl_fit_linear (const double * x, const size_t
  // xstride, const double * y, const size_t ystride, size_t n, double * c0,
  // double * c1, double * cov00, double * cov01, double * cov11, double *
  // sumsq)
  const int n = log_pcounts.size();
  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear(&log_diameters.front(), 1, &log_pcounts.front(), 1, n, &c0,
                 &c1, &cov00, &cov01, &cov11, &sumsq);
  double mean_sq_residual = sumsq / log_diameters.size();
  return {c1, mean_sq_residual};
#else
  runtimeErrorMsg() << "GSL (gnu scientific library) is required for fractal "
                       "dimension calculation.";
  return {0, 0};
#endif
}

} // namespace ClusterAnalysis
