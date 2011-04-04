
#ifndef STATISTICS_OBSERVABLE_H
#define STATISTICS_OBSERVABLE_H

//#include "statistics_correlation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tcl.h"
#include "parser.h"


typedef struct {
  void* args;
  int n;
  int (*fun)  ( void* A_args, double* A, unsigned int dim_A);
  double* last_value;
  int last_update;
} observable;

extern observable** observables;
extern int n_observables; 

int* observable_calc(observable* self, double* A);
int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv);
int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ); 


//int parse_observable(Tcl_Interp* interp, int argc, char** argv, int* change, int (**A_fun)  ( void* A_args, double* A, unsigned int dim_A), int* dim_A, void** A_args);

/* Here we have the particular observables listed */


int observable_particle_velocities(void* idlist, double* A, unsigned int n_A);
int observable_com_velocity(void* idlist, double* A, unsigned int n_A); 
/** Obtain the particle positions.
 * TODO: Folded or unfolded?
 */ 
int observable_particle_positions(void* typelist, double* A, unsigned int n_A);

#ifdef ELECTROSTATICS
int observable_particle_currents(void* typelist, double* A, unsigned int n_A);
int observable_currents(void* typelist, double* A, unsigned int n_A);
#endif


/** Calculate structure factor from positions and scattering length */
int observable_structure_factor(void* params, double* A, unsigned int n_A);
typedef struct {
  int order;
  int dim_sf; // number of q vectors
  int *q_vals; // values of q vectors
  double *q_density; // number of q vectors per bin
  // entries for spherical averaging
} observable_sf_params;

/** See if particles from idList1 interact with any of the particles in idList2 
input parameters are passed via struct iw_params
*/
int observable_interacts_with(void* params, double* A, unsigned int n_A);
typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
} iw_params;


/** Do nothing */
int observable_obs_nothing (void* params, double* A, unsigned int n_A);

int observable_flux_profile(void* params, double* A, unsigned int n_A);
typedef struct { 
  IntList* id_list;
  double minz;
  double maxz;
  int nbins;
} profile_data;

int observable_density_profile(void* params, double* A, unsigned int n_A);

int observable_lb_velocity_profile(void* params, double* A, unsigned int n_A);

int observable_radial_density_profile(void* params, double* A, unsigned int n_A);
typedef struct {
  IntList* id_list;
  double maxr;
  double minz;
  double maxz;
  double center[3];
  int nbins;
} radial_profile_data;



#endif
