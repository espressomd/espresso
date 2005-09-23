// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.

#include "utils.h"
#include "grid.h"
/** \file molforces.h
 *  Routines to for calculating molecule center of mass forces. 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:cooke@mpip-mainz.mpg.de">Ira</a>
 */
#ifdef MOLFORCES

/** 
    Calculates and applies the force on a molecule that would arise
    from an interaction between the molecule center of mass and a trap
    (eg optical tweezers).  This force is assumed to be a harmonic
    potential with a restlength of zero.  User can adjust the
    parameters for this interaction by setting them for each molecule
    individually. (see \file topology.c ).
*/
MDINLINE void calc_moltrap()
{
  Particle *p;
  Molecule *m;
  int i, np, c;
  Cell *cell;
  int j;
  double trapforce;
  if ( !topo_part_info_synced ) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{ 093 can't calculate moltrap: must execute analyse set topo_part_sync first }");
    return;
  }
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      m = &topology[p[i].p.mol_id];
      for(j = 0; j < 3; j++) {
	/* Calculate the trap force for this coordinate and this
	   particle.  Note that since we are modelling a laser trap
	   the force is evenly distributed among atoms within the
	   molecule */
	trapforce = ((m->com[j]-m->trap_center[j])*m->trap_spring_constant)/(double)(m->part.n);
	if (m->trap_flag & COORD_FIXED(j)) {
	  p[i].f.f[j] -= trapforce/PMASS(p[i]);
	}
      }
    }
  }

}

/** 
    Calculate the center of mass, total mass and total force on all molecules 

*/
MDINLINE void calc_mol_forces_coms ( ) {
  Particle *p;
  int i, np, c;
  Cell *cell;
  int j;
  int mi;
  int mol;
  if ( !topo_part_info_synced ) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{ 093 can't calculate molforces: must execute analyse set topo_part_sync first }");
    return;
  }
  if ( n_nodes > 1 ) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{ 092 molecular forces cannot be calculated for more than 1 node }");
    return;
  }

  /* First reset all molecule masses,forces,centers of mass*/
  for ( mi = 0 ; mi < n_molecules ; mi++ ) {
    topology[mi].mass = 0;
    for ( i = 0 ; i < 3 ; i++) {
      topology[mi].f[i] = 0.0;
      topology[mi].com[i] = 0.0;
    }
  }

  /* Calculate forces, mass and unnormalized center of mass */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      mol = p[i].p.mol_id;
      if ( mol >= n_molecules ) {
	char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{ 094 can't calculate molforces no such molecule as %d }",mol);
	return;
      }
      topology[mol].mass += PMASS(p[i]);
      /* Unfold the particle */
      unfold_position(p[i].r.p,p[i].l.i);
      for ( j = 0 ; j < 3 ; j++ ) {
	topology[mol].f[j] += p[i].f.f[j];
	topology[mol].com[j] += p[i].r.p[j]*PMASS(p[i]); 
      }
      /* Fold the particle back */
      fold_position(p[i].r.p,p[i].l.i);
    }
  }

  /* Final normalisation of centers of mass */
  for ( mi = 0 ; mi < n_molecules ; mi++ ) {
    for ( i = 0 ; i < 3 ; i++) {
      topology[mi].com[i] = topology[mi].com[i]/(double)(topology[mi].mass);
    }
  }
}

#endif
