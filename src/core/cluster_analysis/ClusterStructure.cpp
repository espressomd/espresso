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
#include "ClusterStructure.hpp"
#include "Cluster.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include <algorithm>
#include <stdexcept>
#include <utils/for_each_pair.hpp>

namespace ClusterAnalysis {

ClusterStructure::ClusterStructure() { clear(); }

void ClusterStructure::clear() {
  clusters.clear();
  cluster_id.clear();
  m_cluster_identities.clear();
}

inline bool ClusterStructure::part_of_cluster(const Particle &p) {
  return cluster_id.find(p.p.identity) != cluster_id.end();
}

// Analyze the cluster structure of the given particles
void ClusterStructure::run_for_all_pairs() {
  // clear data structs
  clear();

  // Iterate over pairs
  Utils::for_each_pair(partCfg().begin(), partCfg().end(),
                       [this](const Particle &p1, const Particle &p2) {
                         this->add_pair(p1, p2);
                       });
  merge_clusters();
}

void ClusterStructure::run_for_bonded_particles() {
  clear();
  for (const auto &p : partCfg()) {
    for (auto const &bond : p.bonds()) {
      if (bond.partner_ids().size() == 1) {
        add_pair(p, get_particle_data(bond.partner_ids()[0]));
      }
    }
  }

  merge_clusters();
}

void ClusterStructure::add_pair(const Particle &p1, const Particle &p2) {
  // * check, if there's a neighbor
  //   * No: Then go on to the next particle
  // * Yes: Then if
  //   * One of them belongs to a cluster, give the other one the same cluster
  //     id.
  //   * None of them belongs to a cluster: Give them both a new cluster id
  //   * Both belong to different clusters: Mark the clusters as identical
  //   * so that they can be put together later
  if (!m_pair_criterion) {
    runtimeErrorMsg() << "No cluster criterion defined";
    return;
  }
  // If the two particles are neighbors...
  if (m_pair_criterion->decide(p1, p2)) {

    if // None belongs to a cluster
        ((!part_of_cluster(p1)) && (!part_of_cluster(p2))) {
      // Both particles belong to the same, new cluster
      const int cid = get_next_free_cluster_id();

      // assign the cluster_ids
      cluster_id[p1.p.identity] = cid;
      cluster_id[p2.p.identity] = cid;
    } else if // p2 belongs to a cluster but p1 doesn't
        (part_of_cluster(p2) && !part_of_cluster(p1)) {
      // Give p1 the same cluster id as p2
      cluster_id[p1.p.identity] = find_id_for(cluster_id.at(p2.p.identity));
    } else if // i belongs to a cluster but j doesn't
        (part_of_cluster(p1) && !part_of_cluster(p2)) {
      // give p2 the cluster id from p1
      cluster_id[p2.p.identity] = find_id_for(cluster_id.at(p1.p.identity));
    } else if // Both belong to different clusters
        (part_of_cluster(p1) && part_of_cluster(p2) &&
         cluster_id.at(p1.p.identity) != cluster_id.at(p2.p.identity)) {
      // Clusters of p1 and p2 are one and the same. Add an identity to the list
      // The higher number must be inserted as first value of the pair
      // because the substitutions later have to be done in descending order
      const int cid1 = find_id_for(cluster_id.at(p1.p.identity));
      const int cid2 = find_id_for(cluster_id.at(p2.p.identity));
      if (cid1 > cid2) {
        m_cluster_identities[cid1] = cid2;
      } else if (cid1 < cid2) {
        m_cluster_identities[cid2] = cid1;
      }
      // else do nothing. The clusters are already noted for merging.
      // Connected clusters will be merged later
    }
    // The case for both particles being in the same cluster does not need to be
    // treated.
  }
}

void ClusterStructure::merge_clusters() {
  // Relabel particles according to the cluster identities map
  // Also create empty cluster objects for the final cluster id

  // Collect needed changes in a separate map, as doing the changes on the fly
  // would screw up the iterators
  std::vector<std::pair<int, int>> to_be_changed;

  for (auto it : cluster_id) {
    // particle id is in it.first and cluster id in it.second
    // We change the cluster id according to the cluster identities
    // map
    const int cid = find_id_for(it.second);
    // We note the list of changes here, so we don't modify the map
    // while iterating
    to_be_changed.emplace_back(it.first, cid);
    // Empty cluster object
    if (clusters.find(cid) == clusters.end()) {
      clusters[cid] = std::make_shared<Cluster>();
    }
  }

  // Now act on the changes marke in above iteration
  for (auto it : to_be_changed) {
    cluster_id[it.first] = it.second;
  }

  // Now fill the cluster objects with particle ids
  // Iterate over particles, fill in the cluster map
  // to each cluster particle the corresponding cluster id
  for (auto it : cluster_id) {
    // If this is the first particle in this cluster, instance a new cluster
    // object
    if (clusters.find(it.second) == clusters.end()) {
      clusters[it.second] = std::make_shared<Cluster>();
    }
    clusters[it.second]->particles.push_back(it.first);
  }

  // Sort particles ids in the clusters
  for (const auto &c : clusters) {
    std::sort(c.second->particles.begin(), c.second->particles.end());
  }
}

int ClusterStructure::find_id_for(int x) {
  int tmp = x;
  while (m_cluster_identities.find(tmp) != m_cluster_identities.end()) {
    tmp = m_cluster_identities[tmp];
  }
  return tmp;
}

int ClusterStructure::get_next_free_cluster_id() {
  // iterate over cluster_id'
  int max_seen_cluster = 0;
  for (auto it : cluster_id) {
    int cid = it.second;
    if (max_seen_cluster < cid) {
      max_seen_cluster = cid;
    }
  }
  return max_seen_cluster + 1;
}

} // namespace ClusterAnalysis
