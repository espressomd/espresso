/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

/** \file topology.cpp
 *
 *  This file contains functions for handling the system topology.
 *
 *  For more information see topology.hpp
 *   */

#include "utils.hpp"
#include "topology.hpp"
#include "statistics_chain.hpp"
#include "particle_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "molforces.hpp"

int     n_molecules = -1;
Molecule *topology = NULL;
int topo_part_info_synced = 0;

void realloc_topology(int size)
{
  int m;

  for(m = size ; m < n_molecules; m++) {
    realloc_intlist(&topology[m].part, 0);
  }
  
  topology = (Molecule*)realloc(topology, size*sizeof(Molecule));

  if (n_molecules < 0)
    n_molecules = 0;
  for(m = n_molecules; m < size; m++) {
    init_intlist(&topology[m].part);
#ifdef MOLFORCES
    topology[m].trap_flag = 32; 
    topology[m].noforce_flag = 32;
    topology[m].favcounter = -1;
    topology[m].fav[0] = 0;
    topology[m].fav[1] = 0;
    topology[m].fav[2] = 0;
    topology[m].trap_force[0] = 0;
    topology[m].trap_force[1] = 0;
    topology[m].trap_force[2] = 0;
#endif /*MOLFORCES*/
  }
  n_molecules = size;
  
  topo_part_info_synced = 0;

}

// Parallel function for synchronising topology and particle data
void sync_topo_part_info() {
  int i,j;
  Particle* p;
  for ( i = 0 ; i < n_molecules ; i ++ ) {
    for ( j = 0 ; j < topology[i].part.n ; j++ ) {
      p = local_particles[topology[i].part.e[j]];
      if(!p) { 
	/* Do nothing */ 
      } 
      else {
	p->p.mol_id = i;
      }
    }
  }

  topo_part_info_synced = 1;

}


int set_molecule_trap(int mol_num, int trap_flag,DoubleList *trap_center,double spring_constant, double drag_constant, int noforce_flag, int isrelative) {
#ifdef MOLFORCES
  int i;
  if ( mol_num < n_molecules ) {
#ifdef EXTERNAL_FORCES
    topology[mol_num].trap_flag &= ~COORDS_FIX_MASK;
    topology[mol_num].noforce_flag &= ~COORDS_FIX_MASK;
#endif
    /* set new values */
    topology[mol_num].trap_flag |= trap_flag;
    topology[mol_num].noforce_flag |= noforce_flag;

    for ( i = 0 ; i < trap_center->max ; i++){
      topology[mol_num].trap_center[i] = trap_center->e[i];
    }
    topology[mol_num].trap_spring_constant = spring_constant;
    topology[mol_num].drag_constant = drag_constant;
    topology[mol_num].isrelative = isrelative;
    /* check to see if any molecules are trapped */
    if ((topology[i].trap_flag != 32) && (topology[i].noforce_flag != 32)) {
      IsTrapped = 1;
    }
    return ES_OK;
  }
#endif
  return ES_ERROR;
}




