/** \file myconfig-sample.h

    This is a sample for the file myconfig.h.

    Uncomment any of the following lines to myconfig.h to activate
    the corresponding feature of Espresso. It is recommended to turn
    only those features on that you actually need to optimize the
    performance of Espresso for your problem. For details on these
    features see the user's guide, Sec. 2.6 "Configuration options".

    To access the information on the compilation status of the code
    you are working with in your Espresso Tcl-script, use the
    corresponding \ref tcl_features "Tcl-commands".

    If you add a new feature to Espresso, you also have to add the
    corresponding lines in the function \ref compilation_callback and
    to add documentation in <tt>doc/text/features.doc</tt>.
 
    <b>Responsible:</b> 
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

*/

/**********************************************************************/
/*                       general core features                        */
/**********************************************************************/

/* #define PARTIAL_PERIODIC */
/* #define ELECTROSTATICS */
/* #define ROTATION */
/* #define DIPOLES */
/* #define EXTERNAL_FORCES */
/* #define CONSTRAINTS */
/* #define MASS */
/* #define EXCLUSIONS */
/* #define COMFORCE */
/* #define COMFIXED */
/* #define MOLFORCES */
/* #define BOND_CONSTRAINT */

/**********************************************************************/
/*                        integrator features                         */
/**********************************************************************/

/* #define NEMD */
/* #define NPT */ 
/* #define DPD */
/* #define LB */

/**********************************************************************/
/*                           interactions                             */
/**********************************************************************/

/* #define TABULATED */
/* #define LENNARD_JONES */
/* #define SMOOTH_STEP */
/* #define LJ_WARN_WHEN_CLOSE */
/* #define MORSE */
/* #define LJCOS */
/* #define LJCOS2 */
/* #define BUCKINGHAM */
/* #define SOFT_SPHERE */

/* Note: Activate ONLY ONE bonded angle potential out of the following! */
/* #define BOND_ANGLE_HARMONIC */
/* #define BOND_ANGLE_COSINE */
/* #define BOND_ANGLE_COSSQUARE */

/**********************************************************************/
/*                            debugging                               */
/**********************************************************************/

/* #define ADDITIONAL_CHECKS */

/* #define COMM_DEBUG */
/* #define EVENT_DEBUG */
/* #define INTEG_DEBUG */
/* #define CELL_DEBUG */
/* #define GHOST_DEBUG */
/* #define LATTICE_DEBUG */
/* #define HALO_DEBUG */
/* #define GRID_DEBUG */
/* #define VERLET_DEBUG */
/* #define PARTICLE_DEBUG */
/* #define P3M_DEBUG */
/* #define EWALD_DEBUG */
/* #define FFT_DEBUG */
/* #define RANDOM_DEBUG */
/* #define FORCE_DEBUG */
/* #define THERMO_DEBUG */ 
/* #define LJ_DEBUG */
/* #define MORSE_DEBUG */
/* #define ESR_DEBUG */
/* #define ESK_DEBUG */
/* #define FENE_DEBUG */
/* #define GHOST_FORCE_DEBUG */
/* #define ONEPART_DEBUG 13 */
/* #define STAT_DEBUG */ 
/* #define POLY_DEBUG */
/* #define MOLFORCES_DEBUG */
/* #define MEM_DEBUG */
/* #define MAGGS_DEBUG */
/* #define LB_DEBUG */

/* #define ASYNC_BARRIER */

/* #define MPI_CORE */
/* #define FORCE_CORE */
