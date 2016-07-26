
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

class NeighborCriterion {
  public: 
    virtual bool are_neighbors(const Particle& p1, const Particle& p2) =0;
    virtual std::string name() =0;
};




class Cluster {
  public: 
    /** Particles in the cluster */
    std::vector<int> particles;
    // add a particle (Makes a copy of the original)
    void add_particle(const Particle& p);
 };

// add a particle (Makes a copy of the original)

void Cluster::add_particle(const Particle& p) {
 particles.push_back(p.p.identity);
}



class ClusterStructure {
 public:
  // Container to hold the clusters
  std::map<int,Cluster> clusters;
  // ID of cluster, a particle belongs to
  std::map<int, int> cluster_id;
  // Clusters that turn out to be the same (i.e., if two particles are
  // neighbors that already belong to different clusters)
  std::map<int,int> cluster_identities;
  // Clear data strucutres
  void clear();
  // Analyze the cluster structure for energy and sistance-based
  // criteria
  void analyze_pair();
  // Analyze cluster strcutre based on the presence of a bond as criterion
  void analyze_bonds(NeighborCriterion* nd);
  // Count a pair of particles
  void add_pair(Particle& p1, Particle& p2);
  // Merge clusters which turned out to be one and the same
  void merge_clusters();
  // Neighbor criterion
  void set_criterion(NeighborCriterion* c);
  NeighborCriterion* get_criterion() {return nc;};
 private:
  NeighborCriterion* nc;
  // Follow a chain of cluster identities
  inline int find_id_for(int x);
  //
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
      return sqrt(distance2vec(p1.r.p, p2.r.p, vec21)) <= cut_off;
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
      double vec21[3];
      return bond_exists(&p1,&p2,bond_type);
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

