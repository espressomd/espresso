
// Cluster structure and individual particles


#ifndef CLUSTER_H
#define CLUSTER_H


#include <vector>
#include <valarray>
#include <map>
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "energy_inline.hpp"
#include <string>
#include "grid.hpp"

class NeighborCriterion {
  public: 
    virtual bool are_neighbors(const Particle& p1, const Particle& p2) =0;
    virtual std::string name() =0;
};



/** @brief Represents a single cluster of particles */
class Cluster {
  public: 
    /** @brief Ids of the particles in the cluster */
    std::vector<int> particles;
    /** @brief add a particle to the cluster */
    void add_particle(const Particle& p);
 };

void Cluster::add_particle(const Particle& p) {
 particles.push_back(p.p.identity);
}



/** @brief Holds the result and parameters of a cluster analysis */
class ClusterStructure {
 public:
  /** @brief Map holding the individual clusters. The key is an interger cluster id */
  std::map<int,Cluster> clusters;
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
  void set_pair_criterion(std::shared_ptr<PairCriterion> const &c) {
    m_pair_criterion = c;
  }

  PairCriterion const &pair_criterion() const { return *m_pair_criterion; }

 Private:
  /** @brief Clusters that turn out to be the same during the analysis process
  * (i.e., if two particles are neighbors that already belong to different clusters
  */
  std::map<int,int> m_cluster_identities;

  /** @brief pari criterion which decides whether two particles are neighbors */
  std::shared_ptr<PairCriterion> m_pair_criterion;
  
  /** @brief Consider an individual pair of particles during cluster analysis */
  void add_pair(Particle& p1, Particle& p2);
  /** Merge clusters and populate their structures */
  void merge_clusters();
  NeighborCriterion* nc;
  /** @brief Follow a chain of cluster identities during analysis */
  inline int find_id_for(int x);
  /** @brief Get next free lucster id */
  inline int get_next_free_cluster_id();
  
};



//take cut off value that is input data in Tcl simulation script
// 
class DistanceCriterion : public NeighborCriterion {
  public: 
    DistanceCriterion(double _cut_off) {
      cut_off=_cut_off;
    }
    virtual bool are_neighbors(const Particle& p1, const Particle& p2) {
      double vec21[3];
      get_mi_vector(vec21,p1.r.p, p2.r.p); 
      return (sqrt(sqrlen(vec21)<= cut_off));
    };
    virtual std::string name() { return "distance"; };
    double get_cut_off() {
      return cut_off;
    }
    private:
      double cut_off;
};


class EnergyCriterion : public NeighborCriterion {
  public: 
    EnergyCriterion(double _cut_off) {
      cut_off=_cut_off;
    };
    virtual bool are_neighbors(const Particle& p1, const Particle& p2)  {
      double vec21[3];
      double dist_betw_part =sqrt(distance2vec(p1.r.p, p2.r.p, vec21)) <= cut_off;
      IA_parameters *ia_params = get_ia_param(p1.p.type, p2.p.type);
      return (calc_non_bonded_pair_energy(const_cast<Particle*>(&p1), const_cast<Particle*>(&p2), ia_params,
                       vec21, dist_betw_part,dist_betw_part*dist_betw_part)) >= cut_off;
    };
    virtual std::string name() { return "energy"; };
    double get_cut_off() {
      return cut_off;
    }
    private:
      double cut_off;
};

class BondCriterion : public NeighborCriterion {
  public: 
    BondCriterion(int _bond_type) {
       bond_type=_bond_type;
    };
    virtual bool are_neighbors(const Particle& p1, const Particle& p2) {
      return bond_exists(&p1,&p2,bond_type) || bond_exists(&p2,&p1,bond_type);
    };
    virtual std::string name() { return "bond"; };
    double get_bond_type() {
      return bond_type;
    };
    private:
      double bond_type;
};



ClusterStructure& cluster_analysis();
#endif

