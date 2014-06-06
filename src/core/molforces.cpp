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

#include <mpi.h>
#include "utils.hpp"
#include "grid.hpp"
#include "molforces.hpp"
#include "topology.hpp"
#include "particle_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "forces.hpp"

/** \file molforces.cpp
 *  Routines for calculating and applying trap forces upon molecules.
 *  This trap force can be set to
 *  - a harmonic potential with a restlength of zero on the molecular centre of mass
 *  - a drag on the molecular velocity
 *  - a cancelation of the total force on the molecule (including thermostat forces)
 *  The centre of mass can be fixed to an absolute position or to a relative position in the
 *  simulation box.
 *  The molecular trap forces is distributed evenly upon all particles in a molecule.
 *  (see \ref topology.cpp and \ref molforces.cpp)  
 */

#ifdef MOLFORCES


/* global variable which indicates whether any molecules are trapped */
/* set in  mpi_sync_topo_part_info_slave and set_molecule_trap */
int IsTrapped = 0;

void apply_mol_constraints()
{
  Particle *p;
  int i, np, c, mi;
  Cell *cell;
  int j;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      mi = p[i].p.mol_id;
      for(j = 0; j < 3; j++) {
	/* Applies the trap force for this coordinate and this particle. */  
	p[i].f.f[j] += topology[mi].trap_force[j];
      }
    }
  }
}

/* calculates the force applies by traps on all molecules */
void calc_trap_force()
{
  Molecule *m;
  int mi,j;
#ifdef EXTERNAL_FORCES
  double trappos;
#endif
  

  if ( !topo_part_info_synced ) {
    char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{ 093 can't calculate moltrap: must execute analyse set topo_part_sync first }");
    return;
  } else {
    
    m = &topology[0];
    
    for (mi = 0; mi < n_molecules; mi++) {
      m = &topology[mi];
      for(j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
	if (m->trap_flag & COORD_FIXED(j)) {
	  /* Is the molecule trapped at a certain position in space? */
	  if (m->isrelative == 1) {
	    /* Is the trap set to absolute coordinates... */
	    trappos = m->trap_center[j]*box_l[j];
	  } else {
	    /* or to relative ones? */
	    trappos = m->trap_center[j];
	  }
	  m->trap_force[j] = 0;
	  /* the trap_force holding the molecule to the set position in calculated */
	  m->trap_force[j] += -((m->com[j]-trappos)*m->trap_spring_constant)/(double)(m->part.n);
	  /* the drag force applied to the molecule is calculated */
	  m->trap_force[j] += -(m->v[j]*m->drag_constant)/(double)(m->part.n);
	  /* the force applies by the traps is added to fav */
	  /* favcounter counts how many times we have added the force to fav since last time "analyze mol force" was called */
	  /* upon Espresso initialization it is set to -1 because in the first call of "integrate" there is an extra initial time step */
	  /* calling "analyze mol force" resets favcounter to 0 */
	  if (m->favcounter > -1) m->fav[j] -= m->v[j]*m->drag_constant + (m->com[j]-trappos)*m->trap_spring_constant;
	}
	if (m->noforce_flag & COORD_FIXED(j)) {
	  /* the trap force required to cancel out the total force acting on the molecule is calculated */
	  m->trap_force[j] -= m->f[j]/(double)m->part.n;
	  if (m->favcounter > -1) m->fav[j] -= m->f[j];
	}
#endif
      }
      m->favcounter++;
    }
  }
}

/* A list of trapped molecules present on this node is created (local_trapped_mols)*/
void get_local_trapped_mols (IntList *local_trapped_mols)
{
  int c, i, mol, j, fixed;

  for (c = 0; c < local_cells.n; c++) {
    for(i = 0; i < local_cells.cell[c]->n; i++) {
      mol = local_cells.cell[c]->part[i].p.mol_id;
      if ( mol >= n_molecules ) {
	char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{ 094 can't calculate molforces no such molecule as %d }",mol);
	return;
      }

      /* Check to see if this molecule is fixed */
      fixed =0;
      for(j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
	if (topology[mol].trap_flag & COORD_FIXED(j)) fixed = 1;
	if (topology[mol].noforce_flag & COORD_FIXED(j)) fixed = 1;
#endif
      }  
      if (fixed) {
	/* if this molecule isn't already in local_trapped_mols then add it in */
	if (!intlist_contains(local_trapped_mols,mol)) {
	  realloc_intlist(local_trapped_mols, local_trapped_mols->max + 1);
	  local_trapped_mols->e[local_trapped_mols->max-1] = mol;
	  local_trapped_mols->n = local_trapped_mols->max;
	}
      }
    }
  }
}

/* Calculate forces, mass,  and unnormalized center of mass and velocity*/
/* This is only done for the trapped molecules to save time */
void calc_local_mol_info (IntList *local_trapped_mols)
{
  int mi, i,j, mol;
  Particle *p;
  int np, c;
  Cell *cell;
  int lm;
  int fixed;

  /* First reset all molecule masses,forces,centers of mass*/
  for ( mi = 0 ; mi < n_molecules ; mi++ ) {
    topology[mi].mass = 0;
    for ( i = 0 ; i < 3 ; i++) {
      topology[mi].f[i] = 0.0;
      topology[mi].com[i] = 0.0;
      topology[mi].v[i] = 0.0;
    }
  }

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      mol = p[i].p.mol_id;
      if ( mol >= n_molecules ) {
	char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{ 094 can't calculate molforces no such molecule as %d }",mol);
	return;
      }

      /* Check to see if this molecule is fixed */
      fixed =0;
      for(j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
	if (topology[mol].trap_flag & COORD_FIXED(j)) fixed = 1;
	if (topology[mol].noforce_flag & COORD_FIXED(j)) fixed = 1;
#endif
      }  
      if (fixed) {
	topology[mol].mass += PMASS(p[i]);
	/* Unfold the particle */
	unfold_position(p[i].r.p,p[i].l.i);
	for ( j = 0 ; j < 3 ; j++ ) {
	  topology[mol].f[j] += p[i].f.f[j];
	  topology[mol].com[j] += p[i].r.p[j]*PMASS(p[i]); 
	  topology[mol].v[j] += p[i].m.v[j]*PMASS(p[i]); 
	}
	/* Fold the particle back */
	fold_position(p[i].r.p,p[i].l.i);

      }
    }
  }

  /* Final normalisation of centers of mass and velocity*/
  for ( lm = 0 ; lm < local_trapped_mols->n; lm++ ) {
    mi = local_trapped_mols->e[lm];
    for ( i = 0 ; i < 3 ; i++) {
      topology[mi].com[i] = topology[mi].com[i]/(double)(topology[mi].mass);
      topology[mi].v[i] = topology[mi].v[i]/(double)(topology[mi].mass);
    }
  }

}

/* Receives molecule information from the slave nodes. Combines this information,
   calculates trap forces, and returns information to slave nodes */

void mpi_comm_mol_info(IntList *local_trapped_mols) {
  int i, j, k, mol, count;
  double com[3] = {0,0,0};
  double v[3] = {0,0,0};
  double f[3] = {0,0,0};
  double mass = 0;
  /* number of trapped molecules on each node */
  int *n_local_mols;
  /* sum of all elements of n_local_mols */
  int sum_n_local_mols;
  /* lists of which molecules are on each node in order of ascending node number */
  int *local_mols;
  MPI_Status status;

  n_local_mols = (int *) malloc(n_nodes*sizeof(int));
  sum_n_local_mols = 0;

  /* Everyone tells me how many trapped molecules are on their node */
  for (i=1; i <n_nodes; i++) {
    MPI_Recv(&(n_local_mols[i]),1,MPI_INT,i,99,comm_cart,&status);
  }

  for (i=1; i <n_nodes; i++) {
    sum_n_local_mols += n_local_mols[i];
  }
  local_mols = (int *) malloc(sum_n_local_mols*sizeof(int));

  /* Everyone tells me which trapped molecules are on their node */
  count = 0;
  for (i=1; i <n_nodes; i++) {
    MPI_Recv(&(local_mols[count]),n_local_mols[i],MPI_INT,i,99,comm_cart,&status);
    count += n_local_mols[i];
  }

  /* Initialise the centre of masses, velocities and forces to 0
     except for molecules present on master node which are initialized to the local values on the master node
     The centre of masses and velocities are weighted by the total mass on the master node */
  for (i = 0; i < n_molecules; i++) {
    mol =i;
    if (intlist_contains(local_trapped_mols,i)) {
      for (j = 0; j < 3; j++) {
	topology[mol].com[j] = topology[mol].com[j] * topology[mol].mass;
	topology[mol].v[j] = topology[mol].v[j] * topology[mol].mass;
	topology[mol].f[j] = topology[mol].f[j];
      }
    } else {
      topology[mol].mass = 0;
      for (j = 0; j < 3; j++) {
	topology[mol].com[j] = 0;
	topology[mol].v[j] = 0;
	topology[mol].f[j] = 0;
      }
    }
  }
  
  /* The masses, coms, velocities and forces for trapped molecules are received from the slave nodes.
     They are added into the running sums in topology[mol] */
  count = 0;
  for (i = 1; i < n_nodes; i++) {
    for (j = 0; j < n_local_mols[i]; j++) {
      mol = local_mols[count];
      count += 1;
      MPI_Recv(&mass,1,MPI_DOUBLE,i,99,comm_cart,&status);
      MPI_Recv(com,3,MPI_DOUBLE,i,99,comm_cart,&status);
      MPI_Recv(v,3,MPI_DOUBLE,i,99,comm_cart,&status);
      MPI_Recv(f,3,MPI_DOUBLE,i,99,comm_cart,&status);
      topology[mol].mass = topology[mol].mass + mass;
      for (k = 0; k< 3; k++) {
	topology[mol].com[k] += com[k]*mass;
	topology[mol].v[k] += v[k]*mass;
	topology[mol].f[k] += f[k];
      }
    }    
  }

  /* The centre of masses and velocities are renormalized by the total molecular weights */
  for (mol = 0; mol < n_molecules; mol++) {
    for (k=0;k <3; k ++) {
      topology[mol].com[k] = topology[mol].com[k]/topology[mol].mass;
      topology[mol].v[k] = topology[mol].v[k]/topology[mol].mass;
    }
  }

  /* The force exerted by the traps on the molecules are calculated */
  calc_trap_force();

  /* The molecule information and trap forces are sent back to the slave nodes. */
  count = 0;
  for (i = 1; i < n_nodes ; i++) {
    for (j = 0; j < n_local_mols[i]; j++) {
      mol = local_mols[count];
      count += 1;
      MPI_Send(&(topology[mol].mass),1,MPI_DOUBLE,i,99,comm_cart);
      MPI_Send(topology[mol].com,3,MPI_DOUBLE,i,99,comm_cart);
      MPI_Send(topology[mol].v,3,MPI_DOUBLE,i,99,comm_cart);
      MPI_Send(topology[mol].f,3,MPI_DOUBLE,i,99,comm_cart);
      MPI_Send(topology[mol].trap_force,3,MPI_DOUBLE,i,99,comm_cart);
    }
  }

  free(local_mols);
  free(n_local_mols);

}

/* Send molecule information to the master node.
   Recieve the combined molecule information and the trap forces */

void mpi_comm_mol_info_slave(IntList *local_trapped_mols) {
  int i, mol;
  MPI_Status status;

  /* Tells master how many trapped molecules are on this node */
  MPI_Send(&(local_trapped_mols->n),1,MPI_INT,0,99,comm_cart);

  /* Tells master which trapped molecules are on this node */
  MPI_Send(local_trapped_mols->e,local_trapped_mols->n,MPI_INT,0,99,comm_cart);

  for (i = 0; i < local_trapped_mols->n ; i++) {
    mol = local_trapped_mols->e[i];
    /* Send all the masses and coms of the local molecules to the master node */
    MPI_Send(&(topology[mol].mass),1,MPI_DOUBLE,0,99,comm_cart);
    MPI_Send(topology[mol].com,3,MPI_DOUBLE,0,99,comm_cart);
    MPI_Send(topology[mol].v,3,MPI_DOUBLE,0,99,comm_cart);
    MPI_Send(topology[mol].f,3,MPI_DOUBLE,0,99,comm_cart);
  }

  for (i = 0; i < local_trapped_mols->n ; i++) {
    /* Receive all the masses and coms of the local molecules to the master node including info from other nodes*/
    mol = local_trapped_mols->e[i];
    MPI_Recv(&(topology[mol].mass),1,MPI_DOUBLE,0,99,comm_cart,&status);
    MPI_Recv(topology[mol].com,3,MPI_DOUBLE,0,99,comm_cart,&status);
    MPI_Recv(topology[mol].v,3,MPI_DOUBLE,0,99,comm_cart,&status);
    MPI_Recv(topology[mol].f,3,MPI_DOUBLE,0,99,comm_cart,&status);
    MPI_Recv(topology[mol].trap_force,3,MPI_DOUBLE,0,99,comm_cart,&status);
  }

}

/** 
    Calculate the center of mass, total mass, velocity, total force, and trap force on all trapped molecules 
*/
void calc_mol_info () {

  /* list of trapped molecules on this node */
  IntList local_trapped_mols;

  /* check to see if all the topology information has been synced to the various slave nodes */
  if ( !topo_part_info_synced ) {
    char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{ 093 can't calculate molforces: must execute analyse set topo_part_sync first }");
    return;
  }

  init_intlist(&local_trapped_mols);

  /* Find out which trapped molecules are on this node */
  get_local_trapped_mols(&local_trapped_mols);

  /* Calculate the center of mass, mass, velocity, force of whatever fraction of each trapped molecule is on this node*/
  calc_local_mol_info(&local_trapped_mols);

  /* Communicate all this molecular information between nodes.
     It is all sent to the master node which combines it, calculates the trap forces,
     and sends the information back */
  if (this_node == 0) { 
    mpi_comm_mol_info(&local_trapped_mols);
  } else {
    mpi_comm_mol_info_slave(&local_trapped_mols);
  }

  realloc_intlist(&local_trapped_mols,0);
}

void calc_and_apply_mol_constraints ()
{
  if (IsTrapped) {
    /* the molecular information and trap forces are calculated */
    calc_mol_info();
    /* the trap forces are applied to the particles */
    apply_mol_constraints();
  }
}

#endif
