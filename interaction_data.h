// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
#include "particle_data.h"

/** \name Type codes of bonded interactions
    Enumeration of implemented bonded interactions.
*/
/************************************************************/
/*@{*/

/** Type of bonded interaction is a FENE potential. */
#define BONDED_IA_FENE     0
/** Type of bonded interaction is a angle potential. */
#define BONDED_IA_ANGLE    1
/** Type of bonded interaction is a dihedral potential. */
#define BONDED_IA_DIHEDRAL 2
/** Type of bonded interaction is a HARMONIC potential. */
#define BONDED_IA_HARMONIC 3
/** Type of bonded interaction is a (HARMONIC-LJ) potential. */
#define BONDED_IA_SUBT_LJ_HARM 99
/** Type of bonded interaction is a (FENE-LJ) potential. */
#define BONDED_IA_SUBT_LJ_FENE 98
/** This bonded interaction was not set. */
#define BONDED_IA_NONE    -1

/*@}*/


/** \name Type codes for the type of Coulomb interaction
    Enumeration of implemented methods for the electrostatic
    interaction.
*/
/************************************************************/
/*@{*/

/** Coulomb interation switched off (NONE). */
#define COULOMB_NONE  0
/** Coulomb method is Debye-Hueckel. */
#define COULOMB_DH    1
/** Coulomb method is Debye-Hueckel with parallel separate calculation. */
#define COULOMB_DH_PW 2
/** Coulomb method is P3M. */
#define COULOMB_P3M   3
/** Coulomb method is one-dimensional MMM */
#define COULOMB_MMM1D 4

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
/** Rod-like constraint applied. In addition to the general cylinder the rod also allows for a line charge. */
#define CONSTRAINT_ROD 4
/** maze-like constraint applied */
#define CONSTRAINT_MAZE 5

/*@}*/

/* Data Types */
/************************************************************/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
  /** \name Lennard-Jones with shift */
  /*@{*/
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;
  double LJ_capradius;
  /*@}*/

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

  /** \name Tabulated potential */
  /*@{*/
  int TAB_npoints;
  int TAB_startindex;
  double TAB_minval;
  double TAB_minval2;
  double TAB_maxval;
  double TAB_maxval2;
  double TAB_stepsize;
  /*@}*/  
  
} IA_parameters;

/** \name Compounds for Coulomb interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
typedef struct {
  /** Bjerrum length. */
  double bjerrum;
  /** Method to treat coulomb interaction. See \ref COULOMB_NONE "Type codes for Coulomb" */
  int method;
} Coulomb_parameters;

/** Structure to hold Debye-Hueckel Parameters. */
typedef struct {
  /** bjerrum length. */
  double bjerrum;
  /** Cutoff for Debey-Hueckel interaction. */
  double r_cut;
  /** Debye kappa (inverse Debye length) . */
  double kappa;
  /** bjerrum length times temperature. */
  double prefac;
} Debye_hueckel_params;

/*@}*/

/** Defines parameters for a bonded interaction. */
typedef struct {
  /** bonded interaction type. See \ref BONDED_IA_FENE "Type code for bonded" */
  int type;
  /** (Number of particles - 1) interacting for that type */ 
  int num;
  /** union to store the different bonded interaction parameters. */
  union {
    /** Parameters for FENE Potential */
    struct {
      double k;
      double r;
      double r2;
    } fene;
    /** Parameters for Cosine bend potential */
    struct {
      double bend;
    } angle;
    /** Parameters for dihedral potential */
    struct {
      int dummy;
    } dihedral;
    /** Parameters for HARMONIC Potential */
    struct {
      double k;
      double r;
      double r2;
    } harmonic;
   /** Parameters for HARMONIC-LJ Potential */
    struct {
      double k;
      double r;
      double r2;
    } subt_lj_harm; 
   /** Parameters for FENE-LJ Potential */
    struct {
      double k;
      double r;
      double r2;
    } subt_lj_fene; 
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

/** Parameters for a ROD constraint. */
typedef struct {
  /** center of the cylinder in the x-y plane. */
  double pos[2];
  /** cylinder radius. */
  double rad;
  /** line charge density. Only makes sense if the axis along the rod is
      periodically replicated */
  double lambda;
} Constraint_rod;

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

/** Structure to specify a constraint. */
typedef struct {
  /** type of the constraint. */
  int type;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_rod rod;
    Constraint_maze maze;
  } c;

  /** particle representation of this constraint. Actually needed are only the identity,
      the type and the force. */
  Particle part_rep;
} Constraint;
#endif
/*@}*/


/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int n_particle_types;
/* Number of nonbonded (short range) interactions. Not used so far.*/
extern int n_interaction_types;
/** Array of the interaction parameters. Should be accessed only via
    \ref get_ia_param  */
extern IA_parameters *ia_params;

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

/** Structure containing the Debye-Hueckel parameters. */
extern Debye_hueckel_params dh_params;

/** number of bonded interactions. Not used so far. */
extern int n_bonded_ia;
/** Field containing the paramters of the bonded ia types */
extern Bonded_ia_parameters *bonded_ia_params;

/** Array containing all tabulated forces*/
extern DoubleList tabulated_forces;
/** Array containing all tabulated energies*/
extern DoubleList tabulated_energies;

/** Maximal interaction cutoff (real space/short range interactions). */
extern double max_cut;

/** For the warmup you can cap the singularity of the Lennard-Jones
    potential at r=0. look into the warmup documentation for more
    details (who wants to wite that?).*/
extern double lj_force_cap;

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

/************************************************
 * exported functions
 ************************************************/
/** Function for initializing force and energy tables */
void force_and_energy_tables_init();

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

/** get interaction particles between particle sorts i and j */
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
    parameters for the bonded interaction type. Attention: There is no
    initialization done here. */
void make_bond_type_exist(int type);

/** This function increases the LOCAL ia_params field
    to the given size. Better use
    \ref make_particle_type_exist since it takes care of
    the other nodes.  */
void realloc_ia_params(int nsize);

/** calculates the maximal cutoff of all real space interactions. 
    these are: bonded, non bonded + real space electrostatics. */
void calc_maximal_cutoff();

int checkIfParticlesInteract(int i, int j);
 
char *get_name_of_bonded_ia(int i);
#endif
