#ifndef IA_DATA_H
#define IA_DATA_H
/** \file interaction_data.h

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    undocumented.
 */

#include <tcl.h>
#include "particle_data.h"

/** \name Defines */
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
/** This bonded interaction was not set. */
#define BONDED_IA_NONE    -1

/** Coulomb interation switched off (NONE). */
#define COULOMB_NONE  0
/** Coulomb method is Debye-Hueckel. */
#define COULOMB_DH    1
/** Coulomb method is P3M. */
#define COULOMB_P3M   2
/** Coulomb method is one-dimensional MMM */
#define COULOMB_MMM1D 3

#define CONSTRAINT_NONE 0
#define CONSTRAINT_WAL 1
#define CONSTRAINT_SPH 2
#define CONSTRAINT_CYL 3
#define CONSTRAINT_ROD 4

/*@}*/

/** \name Data Types */
/************************************************************/
/*@{*/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
  /** Lennard-Jones with shift */
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;
  double LJ_capradius;

  /* relaxation potential */
  double ramp_cut;
  double ramp_force;

  /** Lennard-Jones+Cos potential */
  double LJCOS_eps;
  double LJCOS_sig;
  double LJCOS_cut;
  double LJCOS_offset;
  double LJCOS_alfa;
  double LJCOS_beta;
  double LJCOS_rmin;
  
  /** Gay-Berne potential */
  double GB_eps;
  double GB_sig;
  double GB_cut;
  double GB_k1;
  double GB_k2;
  double GB_mu;
  double GB_nu;
  double GB_chi1;
  double GB_chi2;
  
  
} IA_parameters;

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
typedef struct {
  /** Bjerrum length. */
  double bjerrum;
  /** Method to treat coulomb interaction. 
      So far implemented: <ul>
      <li> COULOMB_P3M see \ref p3m.h
      <li> COULOMB_DH  see \ref debye_hueckel.h
      </ul>
   */
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

/** Defines parameters for a bonded interaction. */
typedef struct {
  /** bonded interaction type:  0 = FENE, 1 = ANGLE, 2 = DIHEDRAL, 3 = HARMONIC */
  int type;
  /** (Number of particles - 1) interacting for that type */ 
  int num;
  /** union to store the different bonded interaction parameters. */
  union {
    /* Parameters for FENE Potential */
    struct {
      double k;
      double r;
      double r2;
    } fene;
    /* Parameters for Cosine bend potential */
    struct {
      double bend;
    } angle;
    /* Parameters for dihedral potential */
    struct {
      int dummy;
    } dihedral;
    /* Parameters for HARMONIC Potential */
    struct {
      double k;
      double r;
      double r2;
    } harmonic;
	 

  } p;
} Bonded_ia_parameters;

#ifdef CONSTRAINTS
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

/** Structure to specify a constraint. */
typedef struct {
  /** type of the constraint. */
  int type;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_rod rod;
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

/** Maximal interaction cutoff (real space/short range interactions). */
extern double max_cut;

/** For the warmup you can cap the singularity of the Lennard-Jones
    potential at r=0. look into the warmup documentation for more
    details (who wants to wite that?).*/
extern double lj_force_cap;

#ifdef CONSTRAINTS
/** numnber of constraints. */
extern int n_constraints;
/** field containing constraints. */
extern Constraint *constraints;
#endif

/************************************************
 * exportet functions
 ************************************************/

/** Implementation of the tcl command \ref tcl_inter. This function
    allows to modify the interaction parameters.
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
 

#endif
