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
/** \file statistics_cluster.cpp
 *
 *  This file contains the necklace cluster algorithm. It can be used
 *  to identify the substructures 'pearls' and 'strings' on a linear
 *  chain.
 *  See also \ref statistics_cluster.hpp
 */


#include "statistics_cluster.hpp"
#include "statistics_cluster_tcl.hpp"
#include "parser.hpp"

/** \name Routines */
/************************************************************/
/*@{*/

/* parser for necklace cluster analyzation:
   analyze necklace <pearl_treshold> <back_dist> <space_dist> <first> <length> 
 */
int tclcommand_analyze_parse_necklace(Tcl_Interp *interp, int argc, char **argv)
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
  if( first+length > n_part+1) {
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

/* parser for hole cluster analyzation:
   analyze holes <prob_part_type_number> <mesh_size>.
   Needs feature LENNARD_JONES compiled in. */
int tclcommand_analyze_parse_holes(Tcl_Interp *interp, int argc, char **argv)
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
