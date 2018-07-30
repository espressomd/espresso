/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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


#include "statistics_cluster.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"

/** nullptr terminated linked list of elements of a cluster (indices in particle list) */
ClusterElement *element;
/** Double linked list of \ref statistics_cluster::Cluster */
Cluster        *cluster;
/** first cluster in list of \ref statistics_cluster::Cluster */
Cluster        *first_cluster;
/** last cluster in list of \ref statistics_cluster::Cluster */
Cluster        *last_cluster;


/** \name Routines */
/************************************************************/
/*@{*/


/* HOLE CLUSTER ALGORITHM */

/** test if a mesh point belongs to free (return -1) or occupied (return -2) volume.
Needs feature LENNARD_JONES compiled in. */
int test_mesh_elements(PartCfg & partCfg, double pos[3], int probe_part_type) 
{
#ifdef LENNARD_JONES
  double dist,vec[3];

  for (auto &p: partCfg) {
    IA_parameters *ia_params = get_ia_param(p.p.type,probe_part_type);
    get_mi_vector(vec, pos, p.r.p);
    dist = sqrt(sqrlen(vec));

    if ( dist < (ia_params->LJ_cut+ia_params->LJ_offset) ) return -2;
    
  }
#endif
  return -1;
} 


/** Test which mesh points belong to the free and occupied volume. 
    Free volume is marked by -1 and occupied volume by -2.
    Needs feature LENNARD_JONES compiled in. */
void create_free_volume_grid(PartCfg & partCfg, IntList mesh, int dim[3], int probe_part_type)
{
  int i,ix=0,iy=0,iz=0;
  double pos[3];
  double mesh_c[3];

  for ( i=0; i<3; i++) mesh_c[i] = box_l[i] / (double)dim[i];

  for ( i=0; i<(dim[0]*dim[1]*dim[2]); i++) {
    
    pos[0] = (ix+0.5)*mesh_c[0];
    pos[1] = (iy+0.5)*mesh_c[1];
    pos[2] = (iz+0.5)*mesh_c[2];

    mesh.e[i] = test_mesh_elements(partCfg, pos, probe_part_type);

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
  
  int *tmp = (int *) Utils::malloc( sizeof(int)* (dim[0]*dim[1]*dim[2]));
  int *sizes = (int *) Utils::malloc( sizeof(int)* (dim[0]*dim[1]*dim[2]));
 

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
  (*holes) = (int **) Utils::malloc ( sizeof(int *)*(n+1) );
  for ( i=0; i<=n; i++ ) { 
    (*holes)[i] = (int *) Utils::malloc ( sizeof(int)*(sizes[i]+1) );
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

/*@}*/
