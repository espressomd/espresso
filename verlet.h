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

/************************************************
 * data types
 ************************************************/

/** Verlet pair list. The verlet pair list array is resized using a
    sophisticated (we hope) algorithm to avoid unnecessary resizes.
    Access using \ref realloc_verletList.
*/
typedef struct {
  /** The pair payload (two integers per pair) */
  int *pair;
  /** Number of pairs contained */
  int n_pairs;
  /** Number of pairs that fit in until a resize is needed */
  int max_pairs;
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

/** initialize verlet list structure. */
void verlet_init();

/** fill the verlet table. */
void build_verlet_list();

/** exit verlet list structure. */
void verlet_exit();

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
