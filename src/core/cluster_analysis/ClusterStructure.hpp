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

#ifndef CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP
#define CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP

#include "pair_criteria/pair_criteria.hpp"
#include <map>

#include "Cluster.hpp"
#include "Particle.hpp"
#include "pair_criteria/pair_criteria.hpp"

namespace ClusterAnalysis {

/** @brief Holds the result and parameters of a cluster analysis */
class ClusterStructure {
public:
  ClusterStructure();
  /** @brief Map holding the individual clusters. The key is an integer cluster
   * id */
  std::map<int, std::shared_ptr<Cluster>> clusters;
  /** @brief Map between particle ids and corresponding cluster ids */
  std::map<int, int> cluster_id;
  /** @brief Clear data structures */
  void clear();
  /** @brief Run cluster analysis, consider all particle pairs */
  void run_for_all_pairs();
  /** @brief Run cluster analysis, consider pairs of particles connected by a
   * bonded interaction */
  void run_for_bonded_particles();
  /** Is particle p part of a cluster */
  bool part_of_cluster(const Particle &p);
  /** Sets the pair criterion which decides if two particles are neighbors */
  void
  set_pair_criterion(std::shared_ptr<PairCriteria::PairCriterion> const &c) {
    m_pair_criterion = c;
  }

  PairCriteria::PairCriterion const &pair_criterion() const {
    return *m_pair_criterion;
  }

private:
  /** @brief Clusters that turn out to be the same during the analysis process
   * (i.e., if two particles are neighbors that already belong to different
   * clusters
   */
  std::map<int, int> m_cluster_identities;

  /** @brief pair criterion which decides whether two particles are neighbors */
  std::shared_ptr<PairCriteria::PairCriterion> m_pair_criterion;

  /** @brief Consider an individual pair of particles during cluster analysis */
  void add_pair(const Particle &p1, const Particle &p2);
  /** Merge clusters and populate their structures */
  void merge_clusters();
  /** @brief Follow a chain of cluster identities during analysis */
  inline int find_id_for(int x);
  /** @brief Get next free cluster id */
  inline int get_next_free_cluster_id();
};

} // namespace ClusterAnalysis
#endif
