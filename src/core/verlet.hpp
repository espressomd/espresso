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
#ifndef VERLET_H
#define VERLET_H
/** \file verlet.hpp   
 *
 *  This file contains routines to setup and handle interaction pair
 *  lists (verlet pair lists) for the non bonded interactions. 
 *
 *  For the non-bonded interactions, the domain decomposition force
 *  calculation uses verlet pair lists which contain all particle
 *  pairs with a distance smaller than their maximal interaction
 *  range, based on their types, plus the \ref skin.  This allows one
 *  to use these verlet pair lists for several time steps, as long no
 *  particle has moved further than \ref skin / 2.0. You can tune the
 *  verlet pair algorithm with the variable \ref skin which you can
 *  set via the \ref tclcommand_setmd command. You can also acces the
 *  average number of integration steps the verlet lists have been
 *  reused with \ref tclcommand_setmd \ref verlet_reuse.
 *
 *  The verlet algorithm uses the data type \ref PairList to store
 *  interacting particle pairs.
 *
 *  To use verlet pair lists for the force calculation you can either
 *  use the functions \ref build_verlet_lists and \ref
 *  calculate_verlet_ia or a combination of those two \ref
 *  build_verlet_lists_and_calc_verlet_ia.
 *
 *  For energy and pressure calculations using verlet pair lists use
 *  \ref calculate_verlet_energies and \ref calculate_verlet_virials.
 *
 *  For more information see \ref verlet.cpp "verlet.c".
 */
#include "particle_data.hpp"

/************************************************
 * data types
 ************************************************/

/** Verlet pair list. The verlet pair list array is resized using a
    sophisticated (we hope) algorithm to avoid unnecessary resizes.
    Access using \ref resize_verlet_list.
*/
typedef struct {
  /** The pair payload (two pointers per pair) */
  Particle **pair;
  /** Number of pairs contained */
  int n;
  /** Number of pairs that fit in until a resize is needed */
  int max;
} PairList;

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize a Pair List.
 *  Use with care and ONLY for initialization! */
void init_pairList(PairList *list);

/** Free a Pair List . */
void free_pairList(PairList *list);

/** Fill verlet tables. */
void build_verlet_lists();

/** Nonbonded and bonded force calculation using the verlet list */
void calculate_verlet_ia();

/** Fill verlet tables and Calculate nonbonded and bonded forces. This
    is a combination of \ref build_verlet_lists and
    \ref calculate_verlet_ia.
*/
void build_verlet_lists_and_calc_verlet_ia();

/** Nonbonded and bonded energy calculation using the verlet list */
void calculate_verlet_energies();

/** Nonbonded and bonded pressure calculation using the verlet list
    @param v_comp flag which enables (1) compensation of the velocities required 
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst;
		  naturally it doesn't make sense to use it without NpT. */
void calculate_verlet_virials(int v_comp);

/*@}*/



#endif
