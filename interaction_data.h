#ifndef IA_DATA_H
#define IA_DATA_H
/** \file interaction_data.h

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    undocumented.
 */

#include <tcl.h>
#include "particle_data.h"

/************************************************
 * data types
 ************************************************/

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
} IA_parameters;

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
typedef struct {
  /** Bjerrum length. */
  double bjerrum;
  /** Method to treat coulomb interaction. */
  char method[8];
} Coulomb_parameters;

/** Type of bonded interaction is a FENE potential. */
#define BONDED_IA_FENE     0
/** Type of bonded interaction is a angle potential. */
#define BONDED_IA_ANGLE    1
/** Type of bonded interaction is a dihedral potential. */
#define BONDED_IA_DIHEDRAL 2
/** This bonded interaction was not set. */
#define BONDED_IA_NONE    -1

/** Defines parameters for a bonded interaction. */
typedef struct {
  /** bonded interaction type:  0 = FENE, 1 = ANGLE, 2 = DIHEDRAL */
  int type;
  /** (Number of particles - 1) interacting for that type */ 
  int num;
  /** union to store the different bonded interaction parameters. */
  union {
    /* Parameters for FENE Potential */
    struct {
      double k_fene;
      double r_fene;
    } fene;
    /* Parameters for Cosine bend potential */
    struct {
      double bend;
    } angle;
    /* Parameters for dihedral potential */
    struct {
      int dummy;
    } dihedral;

  } p;
} Bonded_ia_parameters;

/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int n_particle_types;
/** Array of the interaction parameters. Should be accessed only via
    \ref get_ia_param  */
extern IA_parameters *ia_params;
/* Number of nonbonded (short range) interactions. Not used so far.*/
extern int n_interaction_types;

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

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

/************************************************
 * exportet functions
 ************************************************/

/** Implementation of the Tcl function inter. This function
    allows to modify the interaction parameters for nonbonded
    interactions like Lennard-Jones.
 */
int inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Callback for setmd niatypes. */
int niatypes_callback(Tcl_Interp *interp, void *data);

/** Callback for setmd lj_force_cap. */
int lj_force_cap_callback(Tcl_Interp *interp, void *data);

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

#endif
