#ifndef IA_DATA_H
#define IA_DATA_H
#include <tcl.h>

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

  /* relaxation potential */
  double ramp_cut;
  double ramp_force;
} IA_parameters;

/** Type of bonded interaction is a dihedral potential. */
#define BONDED_IA_DIHEDRAL 0
/** Type of bonded interaction is a angle potential. */
#define BONDED_IA_ANGLE    1

/** Defines a parameters for a bonded interaction. Not used so far. */
typedef struct {
  int bonded_ia_type;
  union {
    struct {
      int dummy;
    } dihedral;
    struct {
      int dummy;
    } angle;
  } body;
} Bonded_ia_parameters;

/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int n_particle_types;
/** Array of the interaction parameters. Should be accessed only via
    \ref get_ia_param or \ref safe_get_ia_param. */
extern IA_parameters *ia_params;
/* Number of nonbonded (short range) interactions. Not used so far.*/
extern int n_interaction_types;

/** number of bonded interactions. Not used so far. */
extern int n_bonded_ia;
/** Field containing the paramters of the bonded ia types */
extern Bonded_ia_parameters *bonded_ia_params;

/************************************************
 * functions
 ************************************************/

/** tcl procedure for interaction access */
int inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
/* callback for setmd niatypes */
int niatypes_callback(Tcl_Interp *interp, void *data);

/** get interaction particles between particle sorts i and j */
MDINLINE IA_parameters *get_ia_param(int i, int j) {
  extern IA_parameters *ia_params;
  extern int n_particle_types;
  return &ia_params[i*n_particle_types + j];
}

/** get interaction particles between particle sorts i and j.
    returns NULL if i or j < 0, allocates if necessary */
IA_parameters *safe_get_ia_param(int i, int j);

/** realloc n_particle_types */
void realloc_ia_params(int nsize);

/** initialize interaction parameters */
MDINLINE void initialize_ia_params(IA_parameters *params) {
  params->LJ_eps =
    params->LJ_sig =
    params->LJ_cut =
    params->LJ_shift =
    params->LJ_offset = 0;
  
  params->ramp_cut =
    params->ramp_force = 0;
}

/** copy interaction parameters */
MDINLINE void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
  dst->LJ_eps = src->LJ_eps;
  dst->LJ_sig = src->LJ_sig;
  dst->LJ_cut = src->LJ_cut;
  dst->LJ_shift = src->LJ_shift;
  dst->LJ_offset = src->LJ_offset;

  dst->ramp_cut = src->ramp_cut;
  dst->ramp_force = src->ramp_force;  
}
#endif
