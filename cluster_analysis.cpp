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
#include "cluster_analysis.hpp"
#include "particle_data.hpp"
#include "collision.cpp"
#include "collision.hpp"	
#ifdef COLLISION_DETECTION


void ClusterStructure::clear() {
 clusters.clear();
 cluster_id.clear();
 cluster_identities.clear();
}


// Analyze the cluster structure of the given particles
void ClusterStructure::analyze_pair(Particle* p1, Particle* p2)
{
  // clear data structs
  clear();
  
  // Iterate over all possible pairs (n=i*j particles) - ghosts are not included
  for (int i=0;i<=max_seen_particle;i++) {
    if (! local_particles[i]) continue;
    for (int j=i+1;j<=max_seen_particle;j++) {
      if (! local_particles[j]) continue;
      count(*local_particles[i],*local_particles[j]); // maybe no * 
//      int part1, part2;
//      part1 = p1->p.identity;
//      part2 = p2->p.identity;
//      p1 = local_particle[part1];
//      p2 = local_particle[part2];      
//      //check if bond exist for the two particles
//      if (!(p1==p2)) {
//        if ((bond_exists(p1,p2, collision_params.bond_centers)) or (bond_exists(p2,p1, collision_params.bond_centers))){
          //there exist bond, which means that particles are neighbors 
//          return true;
//          printf("particles %d and %d are neighbors!\n",p1.p.identity, p2.p.identity);
        } //bond exist
        else 
          return false;
      } //particles are not the same


    }
  }
}


void ClusterStructure::add_pair(Particle& p1, Particle& p2) {
// * check, if there's a neighbor
 //   * No: Then go on to the next particle
 // * Yes: Then if
 //   * One of them belongs to a cluster, give the other one the same cluster
 //     id.
 //   * None of them belongs to a cluster: Give them both a new cluster id
 //   * Both belong to different clusters: Mark the clusters as identical 
 //   * so that they can be put together later
  if (! nc) {
  //if (! nc) {
    runtimeErrorMsg() << "No cluster criterion defined"; 
    return;
  }
  if (nc->are_neighbors(p1,p2)) {
     if // None belongs to a cluster
     ((cluster_id.find(p1->p.identity)==cluster_id.end()) && (cluster_id.find(p2->p.identity)==cluster_id.end())
    {  
      // Both particles belong to the same, new cluster
      cid=get_next_free_cluster_id();
      cluster_id[p1->p.identity]=cid;
      cluster_id[p2->p.identity]=cid;
    }
    else if // j belongs to a cluster but i doesn't
      ((cluster_id.find(p1.p.identity)==cluster_id.end()) 
      && 
      (cluster_id[p2.p.identity] !=cluster_id.end()))
    {
     // Give p1 the same cluster id as p2
     cluster_id[p1.p.identity]=find_id_for(cluster_id[p2.p.identity]);
    }
    else if // i belongs to a cluster but j doesn't: 
      ((cluster_id.find(p2.p.identity)==cluster_id.end()) 
      && 
      (cluster_id[p1.p.identity] !=cluster_id.end()))
    {
     // give p2 the cluster id from p1
     cluster_id[p2.p.identity]=find_id_for(cluster_id[p1.p.identity]);
    }
    else if // Both belong to different clusters - merge clusters
      (cluster_id[p1.p.identity] != cluster_id[p2.p.identity])
    {
     // Clusters of p1 and p2 are one and the same. Add an identity to the list
     // The lower number must be inserted as first value of tjhe pair
     // because the substituions later have to be done in ascending order
     if (cluster_id[p1.p.identity]<cluster_id[p2.p.identity])
     {
       cluster_identities[find_id_for(cluster_id[p2.p.identity])] =find_id_for(cluster_id[p1.p.identity]);
     }
     else
     {
       cluster_identities[find_id_for(cluster_id[p1.p.identity])] =find_id_for(cluster_id[p2.p.identity]);
     }
     
     // Connected clusters will be merged later
    }
    // The case for both particles being in the same cluster does not need to be
    // treated.
  }
}


void ClusterStructure::merge_clusters() {
  // Relabel particles according to the cluster identities map
  // Also create empty cluster objects for the final cluster id
  for (auto it : cluster_ids) { 
    // particle id is in it.first and cluster id in it.second
    // We change the cluster id according to the cluster identities
    // map
    cid=find_id_for(it.second);
    it.second =cid;
    // Empty cluster object
    if (clusters.find(cid)==clusters.end()) {
      clusters[cid]=Cluster();
    }
  }
  // Now fill the cluster objects with particle ids
  // Itterate over particles
  for (auto id : cluster_ids) {
    clusters[it.second].particles.push_back(it.first);
  }
}


int ClusterStructure::find_id_for(int x)
{
 while (cluster_identities.find(x)!=cluster_identities.end())
 {
  x =cluster_identities[x];
 }
 return x;
}


//ClusterStructure cluster_structure;
ClusterStructure& cluster_analysis() {
  return cluster_structure;
}


void ClusterStructure::set_criterion(NeighborCriterion c) {
  if (nc)
  {
    delete nc;
    nc=0;
  }
  nc=c;
}



#endif   

 

