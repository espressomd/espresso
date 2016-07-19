/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
  
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
#ifndef CLUSTER_ANALYSIS_H
#define CLUSTER_ANALYSIS_H

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "virtual_sites_relative.hpp"
#include "virtual_sites.hpp"
#include "collision.cpp"

#include <vector>
#include <valarray>
#include <map>


/** Particles which are clustered into */
class Cluster {
  public:
    /** All particles in a cluster */
    std::vector<int> particles;
    // adds a particle (Makes a copy of the original)
    inline void add_particle(Particle& p);
 };

// Add a particle to the array particle
void Cluster::add_particle(Particle& p) 
{
 particles.push_back(p.p.identity);
}


class ClusterStructure {
 public:
  // Container to hold the clusters : map associative array composed of a collection KEY-->VALUE pairs, such that each possible key appears at most once in the collection
  std::map<int,Cluster> clusters;
  // ID of cluster, a particle belongs to
  std::map<int, int> cluster_id;
  // Clusters that turn out to be the same (i.e., if two particles are
  // neighbors that already belong to different clusters)
  std::map<int,int> cluster_identities;
  // Clear data structures

  void clear();
  // Analyze the cluster structure for energy and distance-based
  // criteria
  void analyze_pair();

  // Neighbor criterion
  //NeighborCriterion& nc;

  // Analyze cluster structure based on the presence of a bond as criterion
  void analyze_bonds();
  // Count a pair of particles
  void count(Particle& p1, Particle& p2);
  // Merge clusters which turned out to be one and the same
  void merge_clusters();

  // Add pair to the cluster
  void add_pair();  

  // Defines criteria - distance or energy based
  void set_criterion();
 private:
 
  // Follow a chain of cluster identities: reurns the cluster ID's which should be the same
  inline int find_id_for(int x);

};

class NeighborCriterion {
   bool are_neighbors(Particle& p1, Particle& p2);
   bool are_bonded(Particle& p1, Particle& p2); 
}; 
//make here the bonding criteria

/*bool NeighborCriterion::are_neighbors(Particle& p1, Particle& p2)
{
#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((!(p1->p.isVirtual)) || (!(p2->p.isVirtual))) {
#endif
    if (!(p1.p.identity==p2.p.identity)) {
      if ((bond_exists(p1,p2, collision_params.bond_centers)) or (bond_exists(p2,p1, collision_params.bond_centers))){
        //there exist bond, which means that particles are neighbors 
        return true;
        printf("particles are neighbor!");
      } //bond exist
      else 
        return false;
    } //particles are not the same
  } //particles are virtual


}
*/

ClusterStructure& cluster_analysis();
#endif


