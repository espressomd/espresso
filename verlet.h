#ifndef VERLET_H
#define VERLET_H
/** \file verlet.h
    Verlet list. Header file for \ref verlet.c "verlet.c".
 */

/************************************************
 * exported variables
 ************************************************/

/** \name Exported Variables */
/*@{*/
/** Actual number of pairs in the verlet list. */
extern int   n_verletList;
/** Maximal number of pairs in the verlet list. */
extern int max_verletList;
/** Verlet list. */
extern int    *verletList;
/** If non-zero, the verlet list has to be rebuilt. */
extern int rebuild_verletlist;
/*@}*/

/************************************************
 * functions
 ************************************************/

/** \name Exported Functions */
/*@{*/
/** initialize verlet list structure. */
void verlet_init();

/** fill the verlet table. */
void build_verlet_list();

/** exit verlet list structure. */
void verlet_exit();
/*@}*/

#endif
