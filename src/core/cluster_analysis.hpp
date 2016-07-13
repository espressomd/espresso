
// Cluster structure and individual particles


#ifndef CLUSTER_H
#define CLUSTER_H


#include <vector>
#include <valarray>
#include <map>


class Cluster {
  public: 
    /** Particles in the cluster */
    std::vector<int> particles;
    // add a particle (Makes a copy of the original)
    inline void add_particle(Particle& p);
 };

// add a particle (Makes a copy of the original)
void Cluster::add_particle(Particle& p) {
{
 particles.push_back(p.p.identity);
}


class ClusterStructure {
 public:
  // Container to hold the clusters
  map<int,Cluster> clusters;
  // ID of cluster, a particle belongs to
  std::map<int, int> cluster_id;
  // Clusters that turn out to be the same (i.e., if two particles are
  // neighbors that already belong to different clusters)
  std::map<int,int> cluster_identities;
  // Clear data strucutres
  void clear();
  // Analyze the cluster structure for energy and sistance-based
  // criteria
  void analyze_pairs(NeighborCriterion* nd);
  // Analyze cluster strcutre based on the presence of a bond as criterion
  void analyze_bonds(NeighborCriterion* nd);
  // Count a pair of particles
  void count(Particle& p1, Particle& p2);
  // Merge clusters which turned out to be one and the same
  void merge_clusters();
 private:
  // Follow a chain of cluster identities
  inline int find_id_for(int x);
};

ClusterStructure& cluster_analysis();
#endif

