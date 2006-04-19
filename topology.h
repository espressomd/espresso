// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

/** \file topology.h
 *
 *  This file contains functions for handling the system topology.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
 */

#include <tcl.h>
#include "utils.h"

/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/

/** Structure holding information about a molecule */
typedef struct {
  /** Type of the molecule */
  int type;
  /** List of particle identities contained in that molecule */
  IntList part;

#ifdef MOLFORCES
  /** Total force on the molecule */
  double f[3];
  /** Sum of forces on molecule over the last favcounter time steps*/
  double fav[3];
  /** counter for fav */
  int favcounter;
  /** Total mass of the molecule */
  double mass;
  /** Center of mass position */
  double com[3];
  /** velocity of particle*/
  double v[3];
  /** Whether to trap motion in a direction with a harmonic well*/
  int trap_flag;
  /** Location of a harmonic trap for this molecule */
  double trap_center[3];
  /** Trap stiffness */
  double trap_spring_constant;
  /** viscous drag applied to this molecule */
  double drag_constant;
  /** Whether to adjust forces on particles so that net force on molecule is 0 */
  int noforce_flag;
  /** whether trap_center is relative (i.e. fraction of box_length) (1)  or absolute (0)*/
  int isrelative;
  /** the force applied by the trap on the molecule */
  double trap_force[3];
#endif

} Molecule;

/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Number of molecules in the system */
extern int     n_molecules;
/** List of molecules. */
extern Molecule *topology;

extern int topo_part_info_synced;

/*@{*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** reallocate the topology information. All unnecessary entries
    will be freed correctly, all new entries have a zero size
    particle list. */    
void realloc_topology(int new_size);

int parse_analyze_set_topology(Tcl_Interp *interp, int argc, char **argv);

void sync_topo_part_info();
/*@}*/


#endif
