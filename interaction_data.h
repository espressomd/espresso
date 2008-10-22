// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef IA_DATA_H
#define IA_DATA_H
/** \file interaction_data.h
    Various procedures concerning interactions between particles.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    interaction_data.h contains the parser \ref #inter for the
    \ref tcl_inter Tcl command. Therefore the parsing of bonded and nonbonded
    interaction definition commands both is done here. It also contains
    procedures for low-level interactions setup and some helper functions.
    Moreover it contains code for the treatments of constraints.
*/

#include <tcl.h>
#include "utils.h"
#include "particle_data.h" /* needed for constraints */

/** \name Type codes of bonded interactions
    Enumeration of implemented bonded interactions.
*/
/************************************************************/
/*@{*/

/** This bonded interaction was not set. */
#define BONDED_IA_NONE     -1
/** Type of bonded interaction is a FENE potential 
    (to be combined with Lennard Jones). */
#define BONDED_IA_FENE      0
/** Type of bonded interaction is a HARMONIC potential. */
#define BONDED_IA_HARMONIC  1
/** Type of bonded interaction is a bond angle potential. */
#define BONDED_IA_ANGLE     2
/** Type of bonded interaction is a dihedral potential. */
#define BONDED_IA_DIHEDRAL  3
/** Type of tabulated bonded interaction potential, 
    may be of bond length, of bond angle or of dihedral type. */
#define BONDED_IA_TABULATED 4
/** Type of bonded interaction is a (-LJ) potential. */
#define BONDED_IA_SUBT_LJ   5
/** Type of a Rigid/Constrained bond*/
#define BONDED_IA_RIGID_BOND  6
/** Type of a virtual bond*/
#define BONDED_IA_VIRTUAL_BOND  7

/* Specify tabulated bonded interactions  */
#define TAB_UNKNOWN          0
#define TAB_BOND_LENGTH      1
#define TAB_BOND_ANGLE       2
#define TAB_BOND_DIHEDRAL    3

/*@}*/

/** \name Type codes for the type of Coulomb interaction
    Enumeration of implemented methods for the electrostatic
    interaction.
*/
/************************************************************/
/*@{*/

/** Coulomb interation switched off (NONE). */
#define COULOMB_NONE    0
/** Coulomb method is Debye-Hueckel. */
#define COULOMB_DH      1
/** Coulomb method is Debye-Hueckel with parallel separate calculation. */
#define COULOMB_DH_PW   2
/** Coulomb method is P3M. */
#define COULOMB_P3M     3
/** Coulomb method is one-dimensional MMM */
#define COULOMB_MMM1D   4
/** Coulomb method is two-dimensional MMM */
#define COULOMB_MMM2D   5
/** Coulomb method is "Maggs" */
#define COULOMB_MAGGS   6
/** Coulomb method is standard Ewald */
#define COULOMB_EWALD   7
/** Coulomb method is P3M plus ELC. */
#define COULOMB_ELC_P3M 8
/** Coulomb method is Reaction-Field. */
#define COULOMB_RF 9
/** Coulomb method is Reaction-Field BUT as interactions */
#define COULOMB_INTER_RF 10

/*@}*/


/** \name Type codes for constraints
    Enumeration of implemented constraint types.
*/
/************************************************************/
/*@{*/

/** No constraint applied */
#define CONSTRAINT_NONE 0
/** wall constraint applied */
#define CONSTRAINT_WAL 1
/** spherical constraint applied */
#define CONSTRAINT_SPH 2
/** (finite) cylinder shaped constraint applied */
#define CONSTRAINT_CYL 3
/** Rod-like charge. */
#define CONSTRAINT_ROD 4
/** Plate-like charge. */
#define CONSTRAINT_PLATE 5
/** maze-like constraint applied */
#define CONSTRAINT_MAZE 6
/** pore constraint applied */
#define CONSTRAINT_PORE 7
//ER
/** External magnetic field constraint applied */
#define CONSTRAINT_EXT_MAGN_FIELD 8
//end ER
/** Constraint for tunable-lsip boundary conditions */
#define CONSTRAINT_PLANE 9
/*@}*/

/* Data Types */
/************************************************************/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
#ifdef LENNARD_JONES
  /** \name Lennard-Jones with shift */
  /*@{*/
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;
  double LJ_capradius;
  /*@}*/
#endif

#ifdef LENNARD_JONES_GENERIC
  /** \name Generic Lennard-Jones with shift */
  /*@{*/
  double LJGEN_eps;
  double LJGEN_sig;
  double LJGEN_cut;
  double LJGEN_shift;
  double LJGEN_offset;
  double LJGEN_capradius;
  int LJGEN_a1;
  int LJGEN_a2;
  double LJGEN_b1;
  double LJGEN_b2;
  /*@}*/
#endif

#ifdef LJ_ANGLE
  /** \name Directional Lennard-Jones */
  /*@{*/
  double LJANGLE_eps;
  double LJANGLE_sig;
  double LJANGLE_cut;
  /* Locate bonded partners */
  int LJANGLE_bonded1type; 
  int LJANGLE_bonded1pos;
  int LJANGLE_bonded1neg;
  int LJANGLE_bonded2pos;
  int LJANGLE_bonded2neg;
  /* Cap */
  double LJANGLE_capradius;
  /*@}*/
#endif

#ifdef SMOOTH_STEP
  /** \name smooth step potential */
  /*@{*/
  double SmSt_eps;
  double SmSt_sig;
  double SmSt_cut;
  double SmSt_d;
  int    SmSt_n;
  double SmSt_k0;
  /*@}*/
#endif

#ifdef BMHTF_NACL
  /** \name BMHTF NaCl potential */
  /*@{*/
  double BMHTF_A;
  double BMHTF_B;
  double BMHTF_C;
  double BMHTF_D;
  double BMHTF_sig;
  double BMHTF_cut;
  double BMHTF_computed_shift;
  /*@}*/
#endif

#ifdef MORSE 
  /** \name Morse potential */
  /*@{*/
  double MORSE_eps;
  double MORSE_alpha;
  double MORSE_rmin;
  double MORSE_cut;
  double MORSE_rest;
  double MORSE_capradius;
  /*@}*/
#endif

#ifdef BUCKINGHAM
  /** \name Buckingham potential */
  /*@{*/
  double BUCK_A;
  double BUCK_B;
  double BUCK_C;
  double BUCK_D;
  double BUCK_cut;
  double BUCK_discont;
  double BUCK_shift;
  double BUCK_capradius;
  double BUCK_F1;
  double BUCK_F2;
  /*@}*/
#endif

#ifdef SOFT_SPHERE
  /** \name soft-sphere potential */
  /*@{*/
  double soft_a;
  double soft_n;
  double soft_cut;
  double soft_offset;
  /*@}*/
#endif

#ifdef LJCOS
  /** \name Lennard-Jones+Cos potential */
  /*@{*/
  double LJCOS_eps;
  double LJCOS_sig;
  double LJCOS_cut;
  double LJCOS_offset;
  double LJCOS_alfa;
  double LJCOS_beta;
  double LJCOS_rmin;
  /*@}*/
#endif

#ifdef LJCOS2
  /** \name Lennard-Jones with a different Cos potential */
  /*@{*/
  double LJCOS2_eps;
  double LJCOS2_sig;
  double LJCOS2_cut;
  double LJCOS2_offset;
  double LJCOS2_w;
  double LJCOS2_rchange;
  double LJCOS2_capradius;
  /*@}*/
#endif
  
#ifdef ROTATION
  /** \name Gay-Berne potential */
  /*@{*/
  double GB_eps;
  double GB_sig;
  double GB_cut;
  double GB_k1;
  double GB_k2;
  double GB_mu;
  double GB_nu;
  double GB_chi1;
  double GB_chi2;
  /*@}*/  
#endif

#ifdef TABULATED
  /** \name Tabulated potential */
  /*@{*/
  int TAB_npoints;
  int TAB_startindex;
  double TAB_minval;
  double TAB_minval2;
  double TAB_maxval;
  double TAB_maxval2;
  double TAB_stepsize;
  /** The maximum allowable filename length for a tabulated potential file*/
#define MAXLENGTH_TABFILE_NAME 256
  char TAB_filename[MAXLENGTH_TABFILE_NAME];
  /*@}*/  
#endif

#ifdef COMFORCE
  /** \name center of mass directed force */
  /*@{*/
  int COMFORCE_flag;
  int COMFORCE_dir;
  double COMFORCE_force;
  double COMFORCE_fratio;
  /*@}*/
#endif
  
#ifdef COMFIXED
  /** \name center of mass directed force */
  /*@{*/
  int COMFIXED_flag;
  /*@}*/
#endif

#ifdef INTER_DPD
  /** \name DPD as interaction */
  /*@{*/
  double dpd_gamma;
  double dpd_r_cut;
  int dpd_wf;
  double dpd_pref1;
  double dpd_pref2;
  double dpd_tgamma;
  double dpd_tr_cut;
  int dpd_twf;
  double dpd_pref3;
  double dpd_pref4;
  /*@}*/  
#endif

#ifdef INTER_RF
  int rf_on;
#endif

#ifdef MOL_CUT
  int mol_cut_type;
  double mol_cut_cutoff;
#endif
  
#ifdef ADRESS
  /** \name Tabulated potential */
  /*@{*/
  int ADRESS_TAB_npoints;
  int ADRESS_TAB_startindex;
  int ADRESS_IC_npoints;
  double ADRESS_TAB_minval;
  double ADRESS_TAB_minval2;
  double ADRESS_TAB_maxval;
  double ADRESS_TAB_maxval2;
  double ADRESS_TAB_stepsize;
  /** The maximum allowable filename length for a tabulated potential file*/
#define MAXLENGTH_ADRESSTABFILE_NAME 256
  char ADRESS_TAB_filename[MAXLENGTH_ADRESSTABFILE_NAME];
  /*@}*/    
#endif

#ifdef TUNABLE_SLIP
  double TUNABLE_SLIP_temp;
  double TUNABLE_SLIP_gamma;
  double TUNABLE_SLIP_r_cut;
  double TUNABLE_SLIP_time;
  double TUNABLE_SLIP_vx;
  double TUNABLE_SLIP_vy;
  double TUNABLE_SLIP_vz;
#endif

} IA_parameters;

/** \name Compounds for Coulomb interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
typedef struct {
  /** Bjerrum length. */
  double bjerrum;
  /** bjerrum length times temperature. */
  double prefactor;
  
  /** Method to treat coulomb interaction. See \ref COULOMB_NONE "Type codes for Coulomb" */
  int method;

} Coulomb_parameters;

/*@}*/

/** Defines parameters for a bonded interaction. */
typedef struct {
  /** bonded interaction type. See \ref BONDED_IA_FENE "Type code for bonded" */
  int type;
  /** (Number of particles - 1) interacting for that type */ 
  int num;
  /** union to store the different bonded interaction parameters. */
  union {
    /** Parameters for FENE bond Potential.
	k - spring constant.
	drmax - maximal bond streching.
	r0 - equilibrium bond length.
	drmax2 - square of drmax (internal parameter). 
    */
    struct {
      double k;
      double drmax;
      double r0;
      double drmax2;
      double drmax2i;
    } fene;
    /** Parameters for harmonic bond Potential */
    struct {
      double k;
      double r;
      double r_cut;
    } harmonic;
    /** Parameters for three body angular potential (bond-angle potentials). 
	ATTENTION: Note that there are different implementations of the bond angle
	potential which you may chose with a compiler flag in the file \ref config.h !
	bend - bending constant.
	phi0 - equilibrium angle (default is 180 degrees / Pi) */
    struct {
      double bend;
      double phi0;
#ifdef BOND_ANGLE_COSINE
      double cos_phi0;
      double sin_phi0;
#endif
#ifdef BOND_ANGLE_COSSQUARE
      double cos_phi0;
#endif
    } angle;
    /** Parameters for four body angular potential (dihedral-angle potentials). */
    struct {
      double mult;
      double bend;
      double phase;
    } dihedral;
#ifdef TABULATED
    /** Parameters for n-body tabulated potential (n=2,3,4). */
    struct {
      char   *filename;
      int    type;
      int    npoints;
      double minval;
      double maxval;
      double invstepsize;
      double *f;
      double *e;
    } tab;
#endif
    /** Dummy parameters for -LJ Potential */
    struct {
      double k;
      double r;
      double r2;
    } subt_lj;  
    /**Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM*/
    struct {
      /**Length of rigid bond/Constrained Bond*/
      //double d;
      /**Square of the length of Constrained Bond*/
      double d2;
      /**Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE iterations during position corrections*/
      double p_tol;
      /**Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations during velocity corrections */
      double v_tol;
    } rigid_bond;
  } p;
} Bonded_ia_parameters;

#ifdef CONSTRAINTS
/** \name Compounds for constraints */
/*@{*/

/** Parameters for a WALL constraint (or a plane if you like that more). */
typedef struct {
  /** normal vector on the plane. */
  double n[3];
  /** distance of the wall from the origin. */
  double d;
} Constraint_wall;

/** Parameters for a SPHERE constraint. */
typedef struct {
  /** sphere center. */
  double pos[3];
  /** sphere radius. */
  double rad;  
  /** sphere direction. (+1 outside -1 inside interaction direction)*/
  double direction;
} Constraint_sphere;

/** Parameters for a CYLINDER constraint. */
typedef struct {
  /** center of the cylinder. */
  double pos[3];
  /** Axis of the cylinder .*/
  double axis[3];
  /** cylinder radius. */
  double rad;
  /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
  double length;
  /** cylinder direction. (+1 outside -1 inside interaction direction)*/
  double direction;
} Constraint_cylinder;

/** Parameters for a PORE constraint. */
typedef struct {
  /** center of the cylinder. */
  double pos[3];
  /** Axis of the cylinder .*/
  double axis[3];
  /** cylinder radius. */
  double rad;
  /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
  double length;
} Constraint_pore;

/** Parameters for a ROD constraint. */
typedef struct {
  /** center of the cylinder in the x-y plane. */
  double pos[2];
  /** line charge density. Only makes sense if the axis along the rod is
      periodically replicated and only with MMM1D. */
  double lambda;
} Constraint_rod;

/** Parameters for a PLATE constraint. */
typedef struct {
  /** height of plane in z-axis. */
  double pos;
  /** charge density. Only makes sense if the axis along the rod is
      periodically replicated and only with MMM2D. */
  double sigma;
} Constraint_plate;

/** Parameters for a MAZE constraint. */
typedef struct {
  /** number of spheres. */
  double nsphere;
  /** dimension of the maze. */
  double dim;
  /** sphere radius. */
  double sphrad;
  /** cylinder (connecting the spheres) radius*/
  double cylrad;
} Constraint_maze;

//ER
/** Parameters for a EXTERNAL MAGNETIC FIELD constraint */
typedef struct{
  /** vector (direction and magnitude) of the external magnetic field */
  double ext_magn_field[3];
} Constraint_ext_magn_field;
//end ER

/** Parameters for a plane constraint which is needed for tunable-slip boundary conditions. */
typedef struct {
  /** Position of the plain. Negative values mean non-existing in that direction. */
  double pos[3];
} Constraint_plane;

/** Structure to specify a constraint. */
typedef struct {
  /** type of the constraint. */
  int type;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_rod rod;
    Constraint_plate plate;
    Constraint_maze maze;
    Constraint_pore pore;
    //ER
    Constraint_ext_magn_field emfield;
    //end ER
    Constraint_plane plane;
  } c;

  /** particle representation of this constraint. Actually needed are only the identity,
      the type and the force. */
  Particle part_rep;
} Constraint;
/*@}*/
#endif


/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int n_particle_types;
/* Number of nonbonded (short range) interactions. Not used so far.*/
extern int n_interaction_types;

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

/** number of bonded interactions. Not used so far. */
extern int n_bonded_ia;
/** Field containing the paramters of the bonded ia types */
extern Bonded_ia_parameters *bonded_ia_params;

/** Array containing all tabulated forces*/
extern DoubleList tabulated_forces;
/** Array containing all tabulated energies*/
extern DoubleList tabulated_energies;

#ifdef ADRESS
/** Array containing all adress tabulated forces*/
extern DoubleList adress_tab_forces;
/** Array containing all adress tabulated energies*/
extern DoubleList adress_tab_energies;
#endif

/** Maximal interaction cutoff (real space/short range interactions). */
extern double max_cut;
/** Maximal interaction cutoff (real space/short range non-bonded interactions). */
extern double max_cut_non_bonded;

/** For the warmup you can cap the singularity of the Lennard-Jones
    potential at r=0. look into the warmup documentation for more
    details (who wants to wite that?).*/
extern double lj_force_cap;

/** For the warmup you can cap the singularity of the directionnal LJ
    potential at r=0. look into the warmup documentation for more
    details (who wants to write that?).*/
extern double ljangle_force_cap;

/** For the warmup you can cap the singularity of the Morse
    potential at r=0. look into the warmup documentation for more
    details (who wants to wite that?).*/
extern double morse_force_cap;

/** For warm up integration, the maximum force between any two particles
    interacting via Buckingham potential can be set and this magnitude of max
    force is stored in buck_force_cap*/
extern double buck_force_cap;

/** For the warmup you can cap any tabulated potential at the value
    tab_force_cap.  This works for most common potentials where a
    singularity in the force occurs at small separations.  If you have
    more specific requirements calculate a separate lookup table for
    each stage of the warm up.

    \note If the maximum value of the tabulated force at small
    separations is less than the force cap then a warning will be
    issued since the user should provide tabulated values in the range
    where particle interactions are expected.  Even so the program
    will still run and a linear extrapolation will be used at small
    separations until the force reaches the capped value or until zero
    separation */
extern double tab_force_cap;

#ifdef CONSTRAINTS
/** numnber of constraints. */
extern int n_constraints;
/** field containing constraints. */
extern Constraint *constraints;
#endif
/**Switch for nonbonded interaction exclusion*/
extern int ia_excl;

/************************************************
 * exported functions
 ************************************************/
/** Function for initializing force and energy tables */
void force_and_energy_tables_init();

#ifdef ADRESS
/** Function for initializing adress force and energy tables */
void adress_force_and_energy_tables_init();
#endif

/** Implementation of the tcl command \ref tcl_inter. This function
    allows the interaction parameters to be modified.
 */
int inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implementation of the Tcl function constraint. This function
    allows to set and delete constraints.
 */
int constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);

/** Callback for setmd niatypes. */
int niatypes_callback(Tcl_Interp *interp, void *data);

/** get interaction parameters between particle sorts i and j */
MDINLINE IA_parameters *get_ia_param(int i, int j) {
  extern IA_parameters *ia_params;
  extern int n_particle_types;
  return &ia_params[i*n_particle_types + j];
}

/** Makes sure that ia_params is large enough to cover interactions
    for this particle type. The interactions are initialized with values
    such that no physical interaction occurs. */
void make_particle_type_exist(int type);

/** Makes sure that \ref bonded_ia_params is large enough to cover the
    parameters for the bonded interaction type. Attention: 1: There is
    no initialization done here. 2: Use only in connection with
    creating new or overwriting old bond types*/
void make_bond_type_exist(int type);

/** This function increases the LOCAL ia_params field
    to the given size. Better use
    \ref make_particle_type_exist since it takes care of
    the other nodes.  */
void realloc_ia_params(int nsize);

/** calculates the maximal cutoff of all real space
    interactions. these are: bonded, non bonded + real space
    electrostatics. The result is stored in the global variable
    max_cut. The maximal cutoff of the non-bonded + real space
    electrostatic interactions is stored in max_cut_non_bonded. This
    value is used in the verlet pair list algorithm (see \ref
    verlet.h). */
void calc_maximal_cutoff();

/** check whether all force calculation routines are properly initialized. */
int check_obs_calc_initialized();

/**  check if a non bonded interaction is defined */
int checkIfInteraction(IA_parameters *data);

/** check if the types of particles i and j have any non bonded
    interaction defined. */
MDINLINE int checkIfParticlesInteract(int i, int j) {
  return checkIfInteraction(get_ia_param(i, j));
}

char *get_name_of_bonded_ia(int i);
#endif
