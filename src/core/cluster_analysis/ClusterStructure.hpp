
#ifndef CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP
#define CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP


#include <vector>
#include <map>
#include <string>
#include "pair_criteria/pair_criteria.hpp"

#include "pair_criteria/pair_criteria.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "Cluster.hpp"

namespace ClusterAnalysis {


/** @brief Represents a single cluster of particles */
/** @brief Holds the result and parameters of a cluster analysis */
class ClusterStructure {
 public:
  ClusterStructure();
  /** @brief Map holding the individual clusters. The key is an interger cluster id */
  std::map<int,std::shared_ptr<Cluster>> clusters;
  /** @brief Map between particle ids and corresponding cluster ids */
  std::map<int, int> cluster_id;
  /** @brief Clear data strucutres */
  void clear();
  /** @brief Run cluster analysis, consider all aprticle pairs */
  void run_for_all_pairs();
  /** @brief Run cluster analysis, consider pairs of particles connected by a bonded interaction */
  void run_for_bonded_particles();
  /** Is particle p part of a cluster */
  bool part_of_cluster(const Particle& p);
  /** Sets the pair criterion which decides if two particles are neighbors */
  void set_pair_criterion(std::shared_ptr<PairCriteria::PairCriterion> const &c) {
    m_pair_criterion = c;
  }

  PairCriteria::PairCriterion const &pair_criterion() const { return *m_pair_criterion; }

 private:
  /** @brief Clusters that turn out to be the same during the analysis process
  * (i.e., if two particles are neighbors that already belong to different clusters
  */
  std::map<int,int> m_cluster_identities;

  /** @brief pari criterion which decides whether two particles are neighbors */
  std::shared_ptr<PairCriteria::PairCriterion> m_pair_criterion;
  
  /** @brief Consider an individual pair of particles during cluster analysis */
  void add_pair(const Particle& p1, const Particle& p2);
  /** Merge clusters and populate their structures */
  void merge_clusters();
  /** @brief Follow a chain of cluster identities during analysis */
  inline int find_id_for(int x);
  /** @brief Get next free lucster id */
  inline int get_next_free_cluster_id();
};

}
#endif
