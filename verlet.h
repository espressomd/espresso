#ifndef VERLET_H
#define VERLET_H
/** \file verlet.h   Verlet list.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  For more information see \ref verlet.c "verlet.c".
 */

/** \name Exported Variables */
/************************************************************/
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

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** initialize verlet list structure. */
void verlet_init();

/** fill the verlet table. */
void build_verlet_list();

/** exit verlet list structure. */
void verlet_exit();
/*@}*/

#endif
