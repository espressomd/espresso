#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/* fetch the callbacks */
#include "grid.h"
#include "communication.h"

/**********************************************
 * local settings
 **********************************************/

/* increment size of particle buffer */
#define PART_INCREMENT 256

/* decrement size of bonded ia
 * if difference between old length
 * (from recycling) and new length is
 * larger than this, reduce to new size */
#define BONDED_REDUCE 64

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** read only callback. If you choose this, the
    variable cannot be changed from Tcl */
int ro_callback(Tcl_Interp *interp, void *data);

/* callback for npart */
int npart_callback(Tcl_Interp *interp, void *data);
/* callback for nptypes */
int nptypes_callback(Tcl_Interp *interp, void *data);
/* callback for niatypes */
int niatypes_callback(Tcl_Interp *interp, void *data);
/* callback for box_l */
int boxl_callback(Tcl_Interp *interp, void *_data);

const Datafield fields[] = {
  {&nprocs,    TYPE_INT,    1, "nprocs",    ro_callback }, /* communication.c */
  {processor_grid, TYPE_INT, 3, "procgrid", pgrid_callback }, /* grid.c */
  {local_box_l, TYPE_DOUBLE, 3, "local_box_l", ro_callback }, /* global.c */
  {box_l, TYPE_DOUBLE, 3, "box_l", boxl_callback },
  {&n_total_particles, TYPE_INT, 1, "nparticles", npart_callback },
  {&n_particle_types, TYPE_INT, 1, "nptypes", nptypes_callback },
  {&n_interaction_types, TYPE_INT, 1, "niatypes", niatypes_callback },
  {&time_step, TYPE_DOUBLE, 1, "time_step", ro_callback }, /* integrator.c */
  {&max_cut, TYPE_DOUBLE,   1, "max_cut", ro_callback },
  {&skin, TYPE_DOUBLE,   1, "skin", ro_callback },
  {&max_range, TYPE_DOUBLE,   1, "max_range", ro_callback },
  { NULL, 0, 0, NULL, NULL }
};

/**********************************************
 * variables
 **********************************************/

/* simulation box and domain decompostion */ 
double box_l[3]       = {1, 1, 1};
double local_box_l[3] = {1, 1, 1};
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

/* particles */
int n_total_particles = 0;

int     n_particles = 0;
int   max_particles = 0;
Particle *particles = NULL;
int min_free_particle = -1;

int     n_ghosts = 0;
int   max_ghosts = 0;
Particle *ghosts = NULL;
int min_free_ghost = -1;

int *local_index;

/* nonbonded (short range) interactions */
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

/**********************************************
 * procedures
 **********************************************/

void init_data()
{
  max_particles = PART_INCREMENT;
  particles = (Particle *)malloc(sizeof(Particle)*PART_INCREMENT);
}

int got_particle(int part)
{
  int i;
  for (i = 0; i < n_particles; i++)
    if (particles[i].identity == part)
      break;
  if (i == n_particles)
    i = -1;
  return i;
}

int add_particle(int part)
{
  int index = got_particle(part);

  /* particle exists, return */
  if (index >= 0)
    return index;

  if (min_free_particle == -1) {
    /* cannot recycle an old particle */
    index = n_particles++;
    if (max_particles > n_particles) {
      max_particles += PART_INCREMENT;
      particles = (Particle *)
	realloc(particles, sizeof(Particle)*max_particles);
    }

    particles[index].max_bonds = 0;
    particles[index].bonds  = NULL;
  }
  else {
    /* recycling */
    index = min_free_particle;
    for (min_free_particle++; min_free_particle < n_total_particles;
	 min_free_particle++) 
      if (particles[min_free_particle].identity == -1)
	break;
    if (min_free_particle == n_total_particles)
      min_free_particle = -1;
  }

  particles[index].identity = part;

  return index;
}

void free_particle(int index)
{
  particles[index].identity = -1;

  if ((min_free_particle == -1) ||
      (min_free_particle > index))
    min_free_particle = index;
}

void realloc_bonds(int index, int size)
{
  /* reallocate if either too small or too large */
  if ((size  > particles[index].max_bonds) ||
      (size <= particles[index].max_bonds - BONDED_REDUCE)) {
    // is a free
    particles[index].bonds = (int *)
      realloc(particles[index].bonds, sizeof(int)*size);
  }
  particles[index].n_bonds = size;
}

int ro_callback(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

int npart_callback(Tcl_Interp *interp, void *data)
{
  int new_part = *(int *)data;
  if (n_total_particles > new_part) {
    int i;
    /* clean up excess particles. Done on each node */
    for (i = 0; i <  n_particles; i++)
      if (particles[i].identity >= new_part)
	free_particle(i);
  }
  n_total_particles = new_part;
  return (TCL_OK);
}

int nptypes_callback(Tcl_Interp *interp, void *data)
{
  n_particle_types = *(int *)data;
  return (TCL_OK);
}

int niatypes_callback(Tcl_Interp *interp, void *data)
{
  n_interaction_types = *(int *)data;
  return (TCL_OK);
}

int boxl_callback(Tcl_Interp *interp, void *_data)
{
  double *data = _data;

  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  box_l[0] = data[0];
  box_l[1] = data[1];
  box_l[2] = data[2];

  changed_topology();

  return (TCL_OK);
}

void changed_topology()
{
  int i;
  for(i = 0; i < 3; i++) {
    local_box_l[i] = box_l[i]/(double)processor_grid[i]; 
  }
  rebuild_verletlist = 1;
}

IA_parameters *get_ia_param(int i, int j)
{
  return &ia_params[i*n_particle_types + j];
}
