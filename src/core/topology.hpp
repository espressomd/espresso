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
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

/** \file topology.hpp
 *
 *  This file contains functions for handling the system topology.
 */

#include "utils.hpp"

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

void sync_topo_part_info();

int set_molecule_trap(int mol_num, int trap_flag,DoubleList *trap_center,double spring_constant, double drag_constant, int noforce_flag, int isrelative);
/*@}*/


#endif
