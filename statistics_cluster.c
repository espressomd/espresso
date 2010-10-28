/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file statistics_cluster.c
 *
 *  This file contains the necklace cluster algorithm. It can be used
 *  to identify the substructures 'pearls' and 'strings' on a linear
 *  chain.
 *  See also \ref statistics_cluster.h
 */


#include "statistics_cluster.h"

/** \name Data structures */
/************************************************************/
/*@{*/

/** Structure for cluster element */
struct StructClusterElement{  
  /** Element identity */
  int id;
  /** Pointer to next cluster element.*/
  struct StructClusterElement *next;
};

/** Cluster Element */
typedef struct StructClusterElement ClusterElement;

/** Structure for Cluster */
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

/** Cluster */
typedef struct StructCluster Cluster;

/*@}*/

/** \name Variables */
/************************************************************/
/*@{*/

/** NULL terminated linked list of elements of a cluster (indices in particle list) */
ClusterElement *element;
/** Double linked list of \ref statistics_cluster::Cluster */
Cluster        *cluster;
/** first cluster in list of \ref statistics_cluster::Cluster */
Cluster        *first_cluster;
/** last cluster in list of \ref statistics_cluster::Cluster */
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
  // fprintf(stderr," analyze_necklace:\n");
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
  sortPartCfg();
  part = &partCfg[first];

  /* perform necklace cluster algorithm */
  n_pearls = analyze_necklace(part, length) ;

  /* Append result to tcl interpreter */
  sprintf(buffer,"%d",n_pearls);
  Tcl_AppendResult(interp, buffer, " pearls { ", (char *)NULL);
  if( n_pearls > 0 ) {
    cluster = first_cluster;
    sprintf(buffer,"%d",cluster->size);
    Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
    cluster = cluster->next;
    while(cluster->prev != last_cluster) { 
      sprintf(buffer,"%d",cluster->size);
      Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
      cluster = cluster->next;
    }
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

   /* free analyzation memory */
  cluster_free();
  
  return (TCL_OK);
}

/* HOLE CLUSTER ALGORITHM */

/** test if a mesh point belongs to free (return -1) or occupied (return -2) volume.
Needs feature LENNARD_JONES compiled in. */
int test_mesh_element(double pos[3], int probe_part_type) 
{
#ifdef LENNARD_JONES
  int i;
  double dist,vec[3];

  for (i=0; i<n_total_particles; i++) {
    IA_parameters *ia_params = get_ia_param(partCfg[i].p.type,probe_part_type);
    get_mi_vector(vec, pos, partCfg[i].r.p);
    dist = sqrt(sqrlen(vec));

    if ( dist < (ia_params->LJ_cut+ia_params->LJ_offset) ) return -2;
    
  }
#endif
  return -1;
} 


/** Test which mesh points belong to the free and occupied volume. 
    Free volume is marked by -1 and occupied volume by -2.
    Needs feature LENNARD_JONES compiled in. */
void create_free_volume_grid(IntList mesh, int dim[3], int probe_part_type)
{
  int i,ix=0,iy=0,iz=0;
  double pos[3];
  double mesh_c[3];

  for ( i=0; i<3; i++) mesh_c[i] = box_l[i] / (double)dim[i];

  for ( i=0; i<(dim[0]*dim[1]*dim[2]); i++) {
    
    pos[0] = (ix+0.5)*mesh_c[0];
    pos[1] = (iy+0.5)*mesh_c[1];
    pos[2] = (iz+0.5)*mesh_c[2];

    mesh.e[i] = test_mesh_element(pos, probe_part_type);

    ix++; 
    if ( ix >= dim[0]) { ix = ix - dim[0]; iy++; }
    if ( iy >= dim[1]) { iy = iy - dim[1]; iz++; }   
  }  


}

void cluster_neighbors(int point, int dim[3], int neighbors[6])
{
  int x,y,z,a;
  get_grid_pos(point, &x, &y, &z, dim);

  a = x-1; if ( a<0 ) a = dim[0]-1;
  neighbors[0] = get_linear_index(a, y, z, dim);
  a = x+1; if ( a==dim[0] ) a = 0;
  neighbors[1] = get_linear_index(a, y, z, dim);

  a = y-1; if ( a<0 ) a = dim[1]-1;
  neighbors[2] = get_linear_index(x, a, z, dim);
  a = y+1; if ( a==dim[1] ) a = 0;
  neighbors[3] = get_linear_index(x, a, z, dim);

  a = z-1; if ( a<0 ) a = dim[2]-1;
  neighbors[4] = get_linear_index(x, y, a, dim);
  a = z+1; if ( a==dim[2] ) a = 0;
  neighbors[5] = get_linear_index(x, y, a, dim);

}

/** hole cluster algorithm. 
    returns the number of holes and a list of mesh points belonging to each of them */
int cluster_free_volume_grid(IntList mesh, int dim[3], int ***holes)
{
  int i=0,j, k, n=-1, li;
  int neighbors[6];
  
  int *tmp = (int *) malloc( sizeof(int)* (dim[0]*dim[1]*dim[2]));
  int *sizes = (int *) malloc( sizeof(int)* (dim[0]*dim[1]*dim[2]));
 

  // step 1 go through all mesh points
  while ( i < (dim[0]*dim[1]*dim[2]) ) {
    // step 2 test if mesh point is occupied or allready assigned
    if ( mesh.e[i] == -2 || mesh.e[i] >= 0 ) { i++; }
    else {
      // step 3 mesh point is free, create a new cluster
      n++;
      mesh.e[i] = n;
      // start new cluster, put mesh point as first element, and set list pointer on first element
      sizes[n] = 1;
      tmp[0] = i;
      li = 0;
      // step 4 go through all elements of the cluster
      while ( li < sizes[n] ) {
	// step 5 go through all neighbors
	j =  tmp[li];
	cluster_neighbors(j, dim, neighbors);
	for ( k=0; k<6; k++ ) {
	  // step 5 b test if status is free and append it to the cluster
	  if ( mesh.e[ neighbors[k] ] == -1 ) {
	    mesh.e[ neighbors[k] ] = n;
	    // append mesh point as element in the list
	    tmp[ sizes[n] ] = neighbors[k];
	    sizes[n]++;
	  }
	  if (  mesh.e[ neighbors[k] ] > 0 &&  mesh.e[ neighbors[k] ]<n ) {
	    fprintf(stderr,"cfvg:error: i=%d, li=%d, n=%d, mesh.e[x]=%d, x=%d\n",i,li,n,mesh.e[ neighbors[k] ],neighbors[k]); fflush(stderr);
	  }
	}
	li++;
      }
    }
    
 }


  // allocate list space
  (*holes) = (int **) malloc ( sizeof(int *)*(n+1) );
  for ( i=0; i<=n; i++ ) { 
    (*holes)[i] = (int *) malloc ( sizeof(int)*(sizes[i]+1) );
    (*holes)[i][0] = 0;
  }

  for ( i=0; i<(dim[0]*dim[1]*dim[2]); i++ ) {
    j = mesh.e[i];
    if ( j >= 0 ) { 
      (*holes)[j][0] ++;
      (*holes)[j][ (*holes)[j][0] ] = i;
    }
  }

  free(tmp);
  free(sizes);

  return n;
}

/** Calculates the surface to volume ratios of the holes */
void cluster_free_volume_surface(IntList mesh, int dim[3], int nholes, int **holes, int *surface) 
{
  int i, j, n, neighbors[6], inner;

  for ( i=0; i<nholes; i++ ) surface[i] = 0;

  // go through all ellements
  for ( i=0; i<(dim[0]*dim[1]*dim[2]); i++ ) {
    n = mesh.e[i];
    // check all cluster elements
    if ( n >= 0 ) {
      inner = 1;
      cluster_neighbors(i, dim, neighbors);
      // check if all neighbors belong to the same cluster
      for ( j=0; j<6; j++ ) {
	if ( mesh.e[ neighbors[j] ] != n ) inner = 0;
      }
      if ( inner == 0 ) { surface[n]+=1.0; }
    }
  }

}

/* parser for hole cluster analyzation:
   analyze holes <prob_part_type_number> <mesh_size>.
   Needs feature LENNARD_JONES compiled in. */
int parse_hole_cluster_analyzation(Tcl_Interp *interp, int argc, char **argv)
{
  int i,j;
  int probe_part_type;
  int mesh_size=1, meshdim[3];
  double freevol=0.0;
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];

  IntList mesh;

  int n_holes;
  int **holes;
  int max_size=0;
  int *surface;

#ifndef LENNARD_JONES
   Tcl_AppendResult(interp, "analyze holes needs feature LENNARD_JONES compiled in.\n", (char *)NULL);
    return TCL_ERROR;
#endif

  /* check # of parameters */
  if (argc < 2) {
    Tcl_AppendResult(interp, "analyze holes needs 2 parameters:\n", (char *)NULL);
    Tcl_AppendResult(interp, "<prob_part_type_number> <mesh_size>", (char *)NULL);
    return TCL_ERROR;
  }

  /* check parameter types */
  if( (! ARG_IS_I(0, probe_part_type)) ||
      (! ARG_IS_I(1, mesh_size))  ) {
    Tcl_AppendResult(interp, "analyze holes needs 2 parameters of type and meaning:\n", (char *)NULL);
    Tcl_AppendResult(interp, "INT INT\n", (char *)NULL);
    Tcl_AppendResult(interp, "<prob_part_type_number> <mesh_size>", (char *)NULL);
    return TCL_ERROR;
  }

  /* check parameter values */
  if( probe_part_type > n_particle_types || probe_part_type < 0 ) {
    Tcl_AppendResult(interp, "analyze holes: probe particle type number does not exist", (char *)NULL);
    return TCL_ERROR;
  }
  if( mesh_size < 1  ) {
    Tcl_AppendResult(interp, "analyze holes: mesh size must be positive (min=1)", (char *)NULL);
    return TCL_ERROR;
  }

  /* preparation */
  updatePartCfg(WITHOUT_BONDS);
  meshdim[0]=mesh_size;
  meshdim[1]=mesh_size;
  meshdim[2]=mesh_size;
  alloc_intlist(&mesh, (meshdim[0]*meshdim[1]*meshdim[2]));

  /* perform free space identification*/
  create_free_volume_grid(mesh, meshdim, probe_part_type);
  /* perfrom hole cluster algorithm */
  n_holes = cluster_free_volume_grid(mesh, meshdim, &holes);
  /* surface to volume ratio */
  surface = (int *) malloc(sizeof(int)*(n_holes+1));
  cluster_free_volume_surface(mesh, meshdim, n_holes, holes, surface);
  /* calculate accessible volume / max size*/
  for ( i=0; i<=n_holes; i++ ) { 
    freevol += holes[i][0];
    if ( holes[i][0]> max_size ) max_size = holes[i][0];
  }

  /* Append result to tcl interpreter */
  Tcl_AppendResult(interp, "{ n_holes mean_hole_size max_hole_size free_volume_fraction { sizes } { surfaces }  { element_lists } } ", (char *)NULL);

  Tcl_AppendResult(interp, "{", (char *)NULL);

  /* number of holes */
  sprintf(buffer,"%d ",n_holes+1); Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  /* mean hole size */
  sprintf(buffer,"%f ",freevol/(n_holes+1.0)); Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  /* max hole size */
  sprintf(buffer,"%d ",max_size); Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  /* free volume fraction */
  sprintf(buffer,"%f ",freevol/(meshdim[0]*meshdim[1]*meshdim[2]));
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  /* hole sizes */
  Tcl_AppendResult(interp, "{ ", (char *)NULL);
  for ( i=0; i<=n_holes; i++ ) { 
    sprintf(buffer,"%d ",holes[i][0]); Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  /* hole surfaces */
  Tcl_AppendResult(interp, "{ ", (char *)NULL);
  for ( i=0; i<=n_holes; i++ ) { 
    sprintf(buffer,"%d ",surface[i]); Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  /* hole elements */ 
  Tcl_AppendResult(interp, "{ ", (char *)NULL);
  for ( i=0; i<=n_holes; i++ ) { 
    Tcl_AppendResult(interp, "{ ", (char *)NULL);
    for ( j=1; j <= holes[i][0]; j++ ) {
      sprintf(buffer,"%d",holes[i][j]);
      Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  
  Tcl_AppendResult(interp, "}", (char *)NULL);

  /* free allocated memory */
  realloc_intlist(&mesh, 0);
  free(surface);
  for ( i=0; i<=n_holes; i++ ) { free(holes[i]); }
  free(holes);

  return (TCL_OK);
}

/*@}*/
