/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "Cluster.hpp"
#include "PartCfg.hpp"
#include "errorhandling.hpp"
#include "particle_node.hpp"

#include <utils/for_each_pair.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ClusterAnalysis {

void ClusterStructure::clear() {
  clusters.clear();
  cluster_id.clear();
  m_cluster_identities.clear();
}

// Analyze the cluster structure of the given particles
void ClusterStructure::run_for_all_pairs() {
  // clear data structs
  clear();

  // Iterate over pairs
  auto const box_geo_handle = get_box_geo();
  auto const &box_geo = *box_geo_handle;
  PartCfg partCfg{box_geo};
  Utils::for_each_pair(partCfg, [this](Particle const &p1, Particle const &p2) {
    this->add_pair(p1, p2);
  });
  merge_clusters();
}

void ClusterStructure::run_for_bonded_particles() {
  clear();
  auto const box_geo_handle = get_box_geo();
  auto const &box_geo = *box_geo_handle;
  PartCfg partCfg{box_geo};
  for (const auto &p : partCfg) {
    for (auto const bond : p.bonds()) {
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
      cluster_id[p1.id()] = cid;
      cluster_id[p2.id()] = cid;
    } else if // p2 belongs to a cluster but p1 doesn't
        (part_of_cluster(p2) && !part_of_cluster(p1)) {
      // Give p1 the same cluster id as p2
      cluster_id[p1.id()] = find_id_for(cluster_id.at(p2.id()));
    } else if // i belongs to a cluster but j doesn't
        (part_of_cluster(p1) && !part_of_cluster(p2)) {
      // give p2 the cluster id from p1
      cluster_id[p2.id()] = find_id_for(cluster_id.at(p1.id()));
    } else if // Both belong to different clusters
        (part_of_cluster(p1) && part_of_cluster(p2) &&
         cluster_id.at(p1.id()) != cluster_id.at(p2.id())) {
      // Clusters of p1 and p2 are one and the same. Add an identity to the list
      // The higher number must be inserted as first value of the pair
      // because the substitutions later have to be done in descending order
      const int cid1 = find_id_for(cluster_id.at(p1.id()));
      const int cid2 = find_id_for(cluster_id.at(p2.id()));
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
  // Relabel particles according to the cluster identities map.
  // Also create empty cluster objects for the final cluster id.

  // Collect needed changes in a separate map, as doing the changes on the fly
  // would invalidate iterators
  std::vector<std::pair<int, int>> pending_changes;

  for (auto const &[pid, cid] : cluster_id) {
    // We change the cluster id according to the cluster identities map
    auto const merged_cid = find_id_for(cid);
    // We note the list of changes here, so we don't modify the map
    // while iterating
    pending_changes.emplace_back(pid, merged_cid);
    // Initialize new cluster objects
    if (not clusters.contains(merged_cid)) {
      clusters[merged_cid] = std::make_shared<Cluster>(m_box_geo);
    }
  }

  // Now act on the changes marked in above iteration
  for (auto &[pid, merged_cid] : pending_changes) {
    cluster_id[pid] = merged_cid;
  }

  // Now fill the cluster objects with particle ids
  // Iterate over particles, fill in the cluster map
  // to each cluster particle the corresponding cluster id
  for (auto const &[pid, cid] : cluster_id) {
    clusters[cid]->particles.emplace_back(pid);
  }

  // Sort particle ids in the clusters
  for (auto const &cluster : clusters | std::views::values) {
    std::ranges::sort(cluster->particles);
  }
}

int ClusterStructure::find_id_for(int x) const {
  auto cid = x;
  while (m_cluster_identities.contains(cid)) {
    cid = m_cluster_identities.at(cid);
  }
  return cid;
}

int ClusterStructure::get_next_free_cluster_id() {
  int max_id = 0;
  if (not cluster_id.empty()) {
    max_id = *std::ranges::max_element(cluster_id | std::views::values);
  }
  return max_id + 1;
}

} // namespace ClusterAnalysis
