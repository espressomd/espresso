#include "cluster_analysis.hpp"
#include "interaction_data.hpp"


void ClusterStructure::clear() {
 clusters.clear();
 cluster_id.clear();
 cluster_identities.clear();
}

// Analyze the cluster structure of the given particles
void ClusterStructure::analyze_pair()
{
  // clear data structs
  clear();
  
  // Iterate over pairs
  for (int i=0;i<=max_seen_particle;i++) {
    if (! local_particles[i]) continue;
    for (int j=i+1;j<=max_seen_particle;j++) {
      if (! local_particles[j]) continue;
      add_pair(*local_particles[i],*local_particles[j]); // maybe no *
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
    runtimeErrorMsg() << "No cluster criterion defined"; 
    return;
  }
  if (nc->are_neighbors(p1,p2)) {
     if // None belongs to a cluster
     ((cluster_id.find(p1.p.identity)==cluster_id.end()) && (cluster_id.find(p2.p.identity)==cluster_id.end()))
    {  
      // Both particles belong to the same, new cluster
      int cid=get_next_free_cluster_id();

      // assign the 
      cluster_id[p1.p.identity]=cid;
      cluster_id[p2.p.identity]=cid;
    }
    else if // j belongs to a cluster but i doesn't
      ((cluster_id.find(p1.p.identity)==cluster_id.end()) 
      && 
      (cluster_id.find(p2.p.identity) !=cluster_id.end()))
    {
     // Give p1 the same cluster id as p2
     cluster_id[p1.p.identity]=find_id_for(cluster_id[p2.p.identity]);
    }
    else if // i belongs to a cluster but j doesn't
      ((cluster_id.find(p2.p.identity)==cluster_id.end()) 
      && 
      (cluster_id.find(p1.p.identity) !=cluster_id.end()))
    {
     // give p2 the cluster id from p1
     cluster_id[p2.p.identity]=find_id_for(cluster_id[p1.p.identity]);
    }
    else if // Both belong to different clusters
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
  for (auto it : cluster_id) { 
    // particle id is in it.first and cluster id in it.second
    // We change the cluster id according to the cluster identities
    // map
    int cid=find_id_for(it.second);
    it.second =cid;
    // Empty cluster object
    if (clusters.find(cid)==clusters.end()) {
      clusters[cid]=Cluster();
    }
  }

  
  // Now fill the cluster objects with particle ids
  // Iterate over particles, fill in the cluster map 
  // to each cluster particle the corresponding cluster id 
  for (auto it : cluster_id) {
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

int ClusterStructure::get_next_free_cluster_id(){
  //iterate over cluster_id'
  int max_seen_cluster = 1;
  for (auto it : clusters){
    int cid=it.first;
    if (max_seen_cluster <= cid ) {
      max_seen_cluster=cid+1;
    }
  }
  return max_seen_cluster;
}


ClusterStructure cluster_structure;

ClusterStructure& cluster_analysis() {
  return cluster_structure;
}

void ClusterStructure::set_criterion(NeighborCriterion* c) {
  if (nc)
  {
    delete nc;
    nc=0;
  }
  nc=c;
}

   

 

