// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef VERLET_H
#define VERLET_H
/** \file verlet.h   
 *
 *  This file contains routines to setup and handle interaction pair
 *  lists (verlet pair lists) for the non bonded interactions. 
 *
 *  For the non-bonded interactions, the integrator uses verlet pair
 *  lists which contain all particle pairs with a distance smaller
 *  than \ref max_range_non_bonded = \ref max_cut_non_bonded + \ref
 *  skin. This allows one to use these verlet pair lists for several
 *  time steps, as long no particle has moved further than \ref skin /
 *  2.0. You can tune the verlet pair algorithm with the variable \ref
 *  skin which you can set via the \ref setmd command. You can also
 *  acces the average number of integration steps the verlet lists
 *  have been reused with \ref setmd \ref verlet_reuse.
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
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  For more information see \ref verlet.c "verlet.c".
 */
#include <tcl.h>
#include "particle_data.h"

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


/** \name Exported Variables */
/************************************************************/
/*@{*/

/** If non-zero, the verlet list has to be rebuilt. */
extern int rebuild_verletlist;

/*@}*/

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

/** spread the verlet criterion across the nodes. */
void announce_rebuild_vlist();

/** Callback for integrator flag tcl:verletflag c:rebuild_verletlist (= 0 or 1).
    <ul>
    <li> 1 means the integrator rebuilds the verlet list befor the
    first integration step.
    <li> 0 means the integrator reuses the verlet list that it remembers 
    from the last integration step.
    </ul>
    \return TCL status.
*/
int rebuild_vlist_callback(Tcl_Interp *interp, void *_data);
/*@}*/



#endif
