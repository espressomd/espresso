// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file statistics_cluster.c
 *
 *  This file contains the necklace cluster algorithm. It can be used
 *  to identify the substructures 'pearls' and 'strings' on a linear
 *  chain.
 *  See also \ref statistics_cluster.h
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 */


#include "statistics_cluster.h"

/** \name Data structures */
/************************************************************/
/*@{*/


struct StructClusterElement{  
  /** Element identity */
  int id;
  /** Pointer to next cluster element.*/
  struct StructClusterElement *next;
};

typedef struct StructClusterElement ClusterElement;

struct StructCluster {  
  /** Cluster identity */
  int id;
  /** Cluster size */
  int size;
  /** Pointer to first cluster element.*/
  ClusterElement *first;
  /** Pointer to next cluster */
  struct StructCluster *next;
  /** Pointer to previous cluster */
  struct StructCluster *prev;
};

typedef struct StructCluster Cluster;

/*@}*/

/** \name Variables */
/************************************************************/
/*@{*/

/** NULL terminated linked list of elements of a cluster (indices in particle list) */
ClusterElement *element;
/** Double linked list of clusters */
Cluster        *cluster;
/** first cluster in list of \ref cluster */
Cluster        *first_cluster;
/** last cluster in list of \ref cluster */
Cluster        *last_cluster;

/** parameter of necklace cluster algorithm */
int    backbone_distance;
/** parameter of necklace cluster algorithm */
double space_distance2;
/** parameter of necklace cluster algorithm */
int    pearl_treshold;

/*@}*/

/** \name Routines */
/************************************************************/
/*@{*/

/* perform step 1 of necklace cluster algorithm */
void cluster_init(Particle *part, int size) {
  int i;

  element       = (ClusterElement *)malloc(size*sizeof(ClusterElement));
  cluster       = (Cluster *)malloc(size*sizeof(Cluster));
  first_cluster = &cluster[0];
  last_cluster  = &cluster[size-1];
  /* we start with chain_length clusters of size 1 */
  for(i=0;i<size;i++) {
    cluster[i].id    = i;
    cluster[i].size  = 1;
    cluster[i].first = &element[i];
    cluster[i].next  = &cluster[(i+1)%size];
    cluster[i].prev  = &cluster[(size+i-1)%size];
    element[i].id    = i;
    element[i].next  = NULL;
  }
}

void cluster_free() {
  free(element);
  free(cluster);
}

/* join two clusters (Join cluster2 to cluster1 and eliminate cluster2) */
void clusters_join(Cluster *cluster1, Cluster *cluster2)
{
  ClusterElement *element;

  /* join cluster2 to cluster1 */
  element = cluster1->first;
  while(element->next != NULL) element = element->next;
  element->next  = cluster2->first;
  cluster1->size += cluster2->size;
  /* remove cluster2 */
  cluster2->prev->next = cluster2->next;
  cluster2->next->prev = cluster2->prev;
  if(cluster2 == last_cluster)  last_cluster = cluster2->prev;
  if(cluster2 == first_cluster) first_cluster = cluster2->next;
}

/* test wether two clusters are connected (criterion of step 2 of the necklace cluster algorithm) */
int clusters_connected(Particle *part, Cluster cluster1, Cluster cluster2)
{
  int id1,id2;
  ClusterElement *element1,*element2;

  element1 = cluster1.first;
  while ( element1 != NULL ) {
    id1      = element1->id;
    element2 = cluster2.first;
    while ( element2 != NULL ) {
      id2 = element2->id;
      if ( abs( id1 - id2 ) > backbone_distance && 
	   distance2(part[id1].r.p,part[id2].r.p) < space_distance2 ) return 1;
      element2 = element2->next;
    }
    element1 = element1->next;
  }
  return 0;
}

/* test wether two clusters interpenetrate (criterion 4 of the necklace cluster algorithm) */
int clusters_interpenetrate(Cluster cluster1, Cluster cluster2)
{
  int min1,min2,max1,max2;
  ClusterElement *element;

  element = cluster1.first;
  min1 = max1 = element->id;
  while(element->next != NULL) {
    element = element->next;
    if(element->id < min1) min1 = element->id;
    if(element->id > max1) max1 = element->id;
  }

  element = cluster2.first;
  min2 = max2 = element->id;
  while(element->next != NULL) {
    element = element->next;
    if(element->id < min2) min2 = element->id;
    if(element->id > max2) max2 = element->id;
  }

  /* test if one cluster is completelly on left side of the other ... */
  if( ((min1 < min2) && (max1 < min2)) ||
      ((min2 < min1) && (max2 < min1)) ) return 0;
  /* ... otherwise they interpenetrate. */
  return 1;
}

/* perform step 2 to 4 of the necklace cluster algorithm */
int cluster_joincicle(Particle *part) {
  int joins=1, cicles=0;
  Cluster *cluster1,*cluster2;

  while(joins > 0) {
    joins  =0;
    cicles ++;
    cluster1 = first_cluster;
    /* loop cluster 1 until last cluster but one */
    while(cluster1->next != first_cluster) {
      cluster2 = cluster1->next;
      /* loop cluster 2 until lust cluster */
      while(cluster2->prev != last_cluster) {
	if( clusters_connected(part,*cluster1,*cluster2) || 
	    clusters_interpenetrate(*cluster1,*cluster2) ) {
	  clusters_join(cluster1,cluster2);
	  joins++;
	}
	cluster2 = cluster2->next;
      }
      cluster1 = cluster1->next;
    }
  }

  return cicles;
}

/* perform step 5 to 7 of the necklace cluster algorithm 
   \return number of pearls
*/
int cluster_join_to_substructures() 
{
  int n_pearls = 0, previous_is_pearl=1;
  Cluster *cluster1,*cluster2;

  /* search first pearl */
  cluster1 = first_cluster;
  while(cluster1->next != last_cluster  && 
	cluster1->size < pearl_treshold ) cluster1 = cluster1->next;

  /* if no pearl found join all clusters and return 0 ... */
  if(cluster1->size < pearl_treshold) {
    cluster2 = cluster1->next;
    while(cluster2->next != last_cluster) {
      clusters_join(cluster1,cluster2);
      cluster2 = cluster2->next;
    }
    return 0;
  }

  /* ... join all previous clusters to first pearl */
  n_pearls++;
  cluster2 = cluster1->prev;
  while(cluster2->next != first_cluster) {
    clusters_join(cluster1,cluster2);
    cluster2 = cluster2->prev;
  }
  
  /* now go through the rest */
  cluster1 = cluster1->next;
  while(cluster1->prev != last_cluster) {
    /* join actual cluster to previous cluster, if they are of the same type */
    if(cluster1->size < pearl_treshold) {
      /* actual cluster is not a pearl */
      if(previous_is_pearl == 0) clusters_join(cluster1->prev,cluster1);
      else previous_is_pearl = 0;
    }
    else {
      /* actual cluster is a pearl */
      if(previous_is_pearl == 1) clusters_join(cluster1->prev,cluster1);
      else {
	previous_is_pearl = 1;
	n_pearls++;
      }
    }
    cluster1 = cluster1->next;
  }

  /* if last cluster is not a pearl join it to previous cluster */
  while(last_cluster->size < pearl_treshold) clusters_join(last_cluster->prev,last_cluster);
  return n_pearls;
}

/** Perform necklace cluster algorithm. 
    \param   part    pointer to first particle
    \param   np      Number of particles
    \return          Number of pearls in necklace structure
*/
int analyze_necklace(Particle *part, int np) {
  int n_pearls;

  /* initialize: step 1 in necklace cluster analyzation.*/
  cluster_init(part,np);
  /* perform step 2-4 in necklace cluster analyzation.*/
  cluster_joincicle(part);
  /* perform step 5-7 in necklace cluster analyzation.*/
  n_pearls = cluster_join_to_substructures();

  return n_pearls;
}

/* parser for necklace cluster analyzation:
   analyze necklace <pearl_treshold> <back_dist> <space_dist> <first> <length> 
 */
int parse_necklace_analyzation(Tcl_Interp *interp, int argc, char **argv)
{
  double space_dist;
  int first,length;
  Particle *part;
  Cluster *cluster;
  char buffer[TCL_INTEGER_SPACE];
  int n_pearls;

  /* check # of parameters */
  if (argc < 5) {
    Tcl_AppendResult(interp, "analyze necklace needs 5 parameters:\n", (char *)NULL);
    Tcl_AppendResult(interp, "<pearl_treshold> <back_dist> <space_dist> <first> <length>", (char *)NULL);
    return TCL_ERROR;
  }

  /* check parameter types */
  if( (! ARG_IS_I(0, pearl_treshold)) ||
      (! ARG_IS_I(1, backbone_distance))  ||
      (! ARG_IS_D(2, space_dist))     ||
      (! ARG_IS_I(3, first))          ||
      (! ARG_IS_I(4, length)) ) {
    Tcl_AppendResult(interp, "analyze necklace needs 5 parameters of type and meaning:\n", (char *)NULL);
    Tcl_AppendResult(interp, "INT INT DOUBLE INT INT\n", (char *)NULL);
    Tcl_AppendResult(interp, "<pearl_treshold> <back_dist> <space_dist> <first> <length>", (char *)NULL);
    return TCL_ERROR;
  }

  /* check parameter values */
  if( pearl_treshold < 10 ) {
    Tcl_AppendResult(interp, "analyze necklace: pearl_treshold should be >= 10", (char *)NULL);
    return TCL_ERROR;
  }
  if( backbone_distance < 2 ) {
    Tcl_AppendResult(interp, "analyze necklace: backbone_dist should be >= 2", (char *)NULL);
    return TCL_ERROR;
  }
  if( space_dist <= 0.0 ) {
    Tcl_AppendResult(interp, "analyze necklace: space_dist must be positive", (char *)NULL);
    return TCL_ERROR;
  }
  if( first < 0 ) {
    Tcl_AppendResult(interp, "analyze necklace: identity of first particle can not be negative", (char *)NULL);
    return TCL_ERROR;
  }
  if( first+length > n_total_particles+1) {
    Tcl_AppendResult(interp, "analyze necklace: identity of last particle out of partCfg array", (char *)NULL);
    return TCL_ERROR;
  }

  /* preparation */
  space_distance2 = SQR(space_dist);
  part = &partCfg[first];

  /* perform necklace cluster algorithm */
  n_pearls = analyze_necklace(part, length) ;

  /* Append result to tcl interpreter */
  sprintf(buffer,"%d",n_pearls);
  Tcl_AppendResult(interp, buffer, " pearls { ", (char *)NULL);
  cluster = first_cluster;
  while(cluster->prev != last_cluster) { 
    sprintf(buffer,"%d",cluster->size);
    Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
    cluster = cluster->next;
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

   /* free analyzation memory */
  cluster_free();
  
  return (TCL_OK);
}
/*@}*/
