#ifndef VERLET_H
#define VERLET_H
/** \file verlet.h   Verlet list.
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

/** initialize a Pair List.
 *  Use with care and ONLY for initialization! */
void init_pairList(PairList *list);

/** fill verlet tables (old version of build_verlet_lists_and_force_calc()). */
void build_verlet_lists();

/** fill verlet tables and calculate forces. */
void build_verlet_lists_and_force_calc();

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
