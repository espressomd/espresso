/** \file config.h 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/


/** if defined, the code will be slower, but with the \ref periodic
    array you can choose which coordinates are bound to p.b.c and
    which are not. If not defined, all coordinates are bound to
    p.b.c. 

    Has effect on: \ref per_callback, \ref find_node, \ref fields, 
    \ref cells_init and \ref sort_particles_into_cells.
*/
/* #define PARTIAL_PERIODIC */

/** if defined, you will get a warning when particles approach nearer than
    0.9 sigma, because then it's likely the integration will blow up.
*/
/* #define LJ_WARN_WHEN_CLOSE */

#define ELECTROSTATICS


/** Flag to enable external forces. E.g. apply a fixed external force
    to a particle or fix a particle in space. */
#define EXTERNAL_FORCES

/************************************************/
/** \name Default Parameter Settings */
/************************************************/
/*@{*/

/** CELLS: Default value for the maximal number of cells per node. */
#define CELLS_MAX_NUM_CELLS 512

/** P3M: Default for number of interpolation points of the charge
    assignment function. */
#define P3M_N_INTERPOL 32768

/** P3M: Default for baundary condition: Epsilon of the surrounding
    medium. */
#define P3M_EPSILON 0.0

/** P3M: Default for offset of first mesh point from the origin (left
    down corner of the simulation box. */
#define P3M_MESHOFF 0.5

/** P3M: Default for the number of Brillouin zones taken into account
    in the calculation of the optimal influence function (aliasing
    sums). */
#define P3M_BRILLOUIN 1

/** Precision for capture of round off errors. */
#define ROUND_ERROR_PREC 1.0e-14

/*@}*/
