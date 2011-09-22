
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
int observable_com_position(void* idlist, double* A, unsigned int n_A);

#ifdef ELECTROSTATICS
int observable_particle_currents(void* typelist, double* A, unsigned int n_A);
int observable_currents(void* typelist, double* A, unsigned int n_A);
int observable_dipole_moment(void* typelist, double* A, unsigned int n_A);
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

/** For each particle from ids1 get the nearest interaction partner out of 
  those among ids2; 
  If no partner is found within cutoff, set it to -1. 
  Increment the condition when changing from some interaction (before) 
  to some interaction (now).
  Input parameters are passed via struct nn_cond_params
*/
int observable_nearest_neighbour_conditional(void* params, double* A, unsigned int n_A);
typedef struct {
  double cutoff;
  // maximum difference between ids whic makes physical sense
  int chain_length; 
  IntList *ids1;
  IntList *ids2;
  IntList *prev_partners;
  IntList *conditions;
} nn_cond_params;


/** Do nothing */
int observable_obs_nothing (void* params, double* A, unsigned int n_A);

int observable_flux_density_profile(void* params, double* A, unsigned int n_A);
typedef struct { 
  IntList* id_list;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  int xbins;
  int ybins;
  int zbins;
} profile_data;

int observable_density_profile(void* params, double* A, unsigned int n_A);

int observable_lb_velocity_profile(void* params, double* A, unsigned int n_A);

int observable_radial_density_profile(void* params, double* A, unsigned int n_A);
int observable_radial_flux_density_profile(void* params, double* A, unsigned int n_A);
int observable_lb_radial_velocity_profile(void* params, double* A, unsigned int n_A);
typedef struct {
  IntList* id_list;
  double minr;
  double maxr;
  double minphi;
  double maxphi;
  double minz;
  double maxz;
  double center[3];
  double axis[3];
  int phibins;
  int rbins;
  int zbins;
} radial_profile_data;



#endif
