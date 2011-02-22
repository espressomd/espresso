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
    corresponding lines in the function \ref tclcallback_compilation and
    to add documentation in <tt>doc/text/features.doc</tt>.
 */

/**********************************************************************/
/*                       general core features                        */
/**********************************************************************/

/* #define PARTIAL_PERIODIC */
/* #define ELECTROSTATICS */
/* #define ROTATION */
/* #define ROTATIONAL_INTERIA */
/* #define DIPOLES */
/* #define MAGNETOSTATICS */
/* #define EXTERNAL_FORCES */
/* #define CONSTRAINTS */
/* #define MASS */
/* #define EXCLUSIONS */
/* #define COMFORCE */
/* #define COMFIXED */
/* #define MOLFORCES */
/* #define BOND_CONSTRAINT */
/* #define MODES */
/* #define BOND_VIRTUAL */

// To use address, also activate VIRTUAL_SITES and VIRTUAL_SITES_COM
/* #define ADRESS*/

/* #define DAWAANR */
/* #define MAGNETIC_DIPOLAR_DIRECT_SUM */
/* #define MDLC */
/* #define METADYNAMICS */
/* #define OVERLAPPED */

/* Note: Activate only one virtual sites implementation! */
/* #define VIRTUAL_SITES_COM */
/* #define VIRTUAL_SITES_RELATIVE */

/* #define VIRTUAL_SITES_NO_VELOCITY */
/* #define VIRTUAL_SITES_THERMOSTAT */
/* #define THERMOSTAT_IGNORE_NON_VIRTUAL */


/**********************************************************************/
/*                        integrator features                         */
/**********************************************************************/

/* #define NEMD */
/* #define NPT */ 
/* #define DPD */
/* #define TRANS_DPD */
/* #define INTER_DPD */

/*Note: Activate ONLY ONE dpd mass  out of the following! */
/* #define DPD_MASS_RED */
/* #define DPD_MASS_LIN */

/* #define LB */
/* #define LB_ELECTROHYDRODYNAMICS */

/**********************************************************************/
/*                           interactions                             */
/**********************************************************************/

/* #define TABULATED */
/* #define LENNARD_JONES */
/* #define LJ_WARN_WHEN_CLOSE */
/* #define LENNARD_JONES_GENERIC */
/* #define LJCOS */
/* #define LJCOS2 */
/* #define LJ_ANGLE */
/* #define GAY_BERNE */
/* #define SMOOTH_STEP */
/* #define HERTZIAN */
/* #define BMHTF_NACL */
/* #define MORSE */
/* #define BUCKINGHAM */
/* #define SOFT_SPHERE */
/* #define INTER_DPD */
/* #define INTER_RF */
/* #define MOL_CUT */
/* #define TUNABLE_SLIP */
/* #define NO_INTRA_NB */

/* Note: Activate ONLY ONE bonded angle potential out of the following! */
/* #define BOND_ANGLE_HARMONIC */
/* #define BOND_ANGLE_COSINE */
/* #define BOND_ANGLE_COSSQUARE */

/* #define BOND_ANGLEDIST */
/* #define BOND_ENDANGLEDIST */

/* activate the old dihedral form. Only uncomment this if you need
   to run old code using the previous phase definition. */
/* #define OLD_DIHEDRAL */

/**********************************************************************/
/* rarely used constants. Normally, you don't need to touch these.    */
/* Change at own risk!                                                */
/**********************************************************************/

/** P3M: Number of Brillouin zones taken into account
    in the calculation of the optimal influence function (aliasing
    sums). */
//#define P3M_BRILLOUIN 1
/** P3M: Maximal mesh size that will be checked. The current setting
         limits the memory consumption to below 1GB, which is probably
	 reasonable for a while. */
//#define P3M_MAX_MESH 128

/** Whether to use the approximation of Abramowitz/Stegun
    AS_erfc_part() for \f$\exp(d^2) erfc(d)\f$, or the C function erfc
    in P3M and Ewald summation. */
//#define USE_ERFC_APPROXIMATION 1

/** Precision for capture of round off errors. */
//#define ROUND_ERROR_PREC 1.0e-14

/** Tiny angle cutoff for sinus calculations */
//#define TINY_SIN_VALUE 1e-10
/** Tiny angle cutoff for cosine calculations */
//#define TINY_COS_VALUE 0.9999999999
/** Tiny length cutoff */
//#define TINY_LENGTH_VALUE 0.0001

/** maximal number of iterations in the RATTLE algorithm before it bails out. */
//#define SHAKE_MAX_ITERATIONS 1000

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
/* #define PTENSOR_DEBUG*/
/* #define MEM_DEBUG */
/* #define MAGGS_DEBUG */
/* #define LB_DEBUG */
/* #define VIRTUAL_SITES_DEBUG */

/* #define ASYNC_BARRIER */

/* #define MPI_CORE */
/* #define FORCE_CORE */

/* #define OLD_RW_VERSION*/
