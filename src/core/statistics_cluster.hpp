/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef STATISTICS_CLUSTER_H
#define STATISTICS_CLUSTER_H
/** \file statistics_cluster.hpp
 *
 *  1: This file contains the necklace cluster algorithm. It can be used
 *  to identify the substructures 'pearls' and 'strings' on a linear
 *  chain.
 *
 *  2: mesh based cluster algorithm to identify hole spaces 
 *  (see thesis chapter 3 of H. Schmitz for details) 
 */

#include "interaction_data.hpp"
#include "particle_data.hpp"

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
extern ClusterElement *element;
/** Double linked list of \ref statistics_cluster::Cluster */
extern Cluster        *cluster;
/** first cluster in list of \ref statistics_cluster::Cluster */
extern Cluster        *first_cluster;
/** last cluster in list of \ref statistics_cluster::Cluster */
extern Cluster        *last_cluster;

/** parameter of necklace cluster algorithm */
extern int    backbone_distance;
/** parameter of necklace cluster algorithm */
extern double space_distance2;
/** parameter of necklace cluster algorithm */
extern int    pearl_treshold;

/*@}*/

void cluster_free();
void create_free_volume_grid(IntList mesh, int dim[3], int probe_part_type);
int analyze_necklace(Particle *part, int np);
int cluster_free_volume_grid(IntList mesh, int dim[3], int ***holes);
void cluster_free_volume_surface(IntList mesh, int dim[3], int nholes, int **holes, int *surface);

#endif
