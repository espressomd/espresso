#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "debug.h"

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

/* callback for nptypes */
int nptypes_callback(Tcl_Interp *interp, void *data);
/* callback for niatypes */
int niatypes_callback(Tcl_Interp *interp, void *data);
/* callback for box_l */
int boxl_callback(Tcl_Interp *interp, void *_data);

/* do not change order !!! */
const Datafield fields[] = {
  {&nprocs,    TYPE_INT,    1, "nprocs",    ro_callback }, /* communication.c */
#define FIELD_NPROCS 0
  {processor_grid, TYPE_INT, 3, "procgrid", pgrid_callback }, /* grid.c */
#define FIELD_PGRID  1
  {local_box_l, TYPE_DOUBLE, 3, "local_box_l", ro_callback }, /* global.c */
#define FIELD_LBOXL  2
  {box_l, TYPE_DOUBLE, 3, "box_l", boxl_callback },
#define FIELD_BOXL   3
  {&n_total_particles, TYPE_INT, 1, "nparticles", ro_callback },
#define FIELD_NTOTAL 4
  {&n_particle_types, TYPE_INT, 1, "nptypes", ro_callback },
#define FIELD_NPTYPE 5
  {&n_interaction_types, TYPE_INT, 1, "niatypes", niatypes_callback },
#define FIELD_NITYPE 6
  {&time_step, TYPE_DOUBLE, 1, "time_step", ro_callback }, /* integrator.c */
#define FIELD_TSTEP  7
  {&max_cut, TYPE_DOUBLE,   1, "max_cut", ro_callback },
#define FIELD_MCUT   8
  {&skin, TYPE_DOUBLE,   1, "skin", ro_callback },
#define FIELD_SKIN   9
  {&max_range, TYPE_DOUBLE,   1, "max_range", ro_callback },
#define FIELD_RANGE 10
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
int *particle_node = NULL;

int     n_particles = 0;
int   max_particles = 0;
int        n_ghosts = 0;
Particle *particles = NULL;

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

void map_particle_node(int part, int node)
{
  int old_max = n_total_particles, i;
  if (part >= n_total_particles) {
    n_total_particles = PART_INCREMENT*((part + PART_INCREMENT)/PART_INCREMENT);
    particle_node = realloc(particle_node, sizeof(int)*n_total_particles);
    for (i = old_max; i < n_total_particles; i++)
      particle_node[i] = -1;
  }
  PART_TRACE(fprintf(stderr, "mapping %d -> %d (%d)\n", part, node, n_total_particles));
  particle_node[part] = node;
}

void build_particle_node()
{
  if (n_total_particles == 0)
    return;
  if (particle_node)
    free(particle_node);
  particle_node = malloc(n_total_particles*sizeof(int));
  mpi_who_has();
}

void realloc_particles(int size)
{
  int old_max = max_particles, i;
  if (size < max_particles) {
    /* shrink not as fast, just lose half, rounded up */
    max_particles = PART_INCREMENT*(((n_particles + size)/2 +
				     PART_INCREMENT - 1)/PART_INCREMENT);
  }
  else
    /* round up */
    max_particles = PART_INCREMENT*((size + PART_INCREMENT - 1)/PART_INCREMENT);
  if (max_particles != old_max)
    particles = (Particle *) realloc(particles, sizeof(Particle)*max_particles);
  for (i = old_max; i < max_particles; i++)
    particles[i].identity = -1;
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
  int index;
  if ((index = got_particle(part)) != -1)
    return index;
  index = alloc_particle();
  particles[index].identity = part;
  return index;
}

int alloc_particle()
{
  int i,index;

  /* add at end */
  index = n_particles++;

  realloc_particles(n_particles);
    
  particles[index].n_bonds = 0;
  particles[index].max_bonds = 0;
  particles[index].bonds  = NULL;
  for(i = 0; i < 3; i++)
    particles[index].i[i] = 0;
  return index;
}

void fold_particle(double pos[3],int image_box[3])
{
  int i;
  for(i=0;i<3;i++) {
    image_box[i] += floor(pos[i]/box_l[i]);
    pos[i]       = pos[i] - image_box[i]*box_l[i];    
  }
}

void unfold_particle(double pos[3],int image_box[3])
{
  int i;
  for(i=0;i<3;i++) {
    pos[i]       = pos[i] + image_box[i]*box_l[i];    
    image_box[i] = 0;
  }
}

void particle_finalize_data()
{
  /* if this is zero, n_total_particles didn't change */
  if (!particle_node)
    return;

  /* calculate n_total_particles */
  while (particle_node[n_total_particles - 1] == -1)
    n_total_particles--;

  mpi_bcast_parameter(FIELD_NTOTAL);

  /* invalidate particle->node data */
  if (particle_node) {
    free(particle_node);
    particle_node = NULL;
  }
}

void realloc_bonds(int index, int size)
{
  /* reallocate if either too small or too large */
  if ((size  > particles[index].max_bonds) ||
      (size <= particles[index].max_bonds - BONDED_REDUCE)) {
    particles[index].max_bonds = size;
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

  mpi_bcast_parameter(FIELD_LBOXL);

  rebuild_verletlist = 1;
}

IA_parameters *get_ia_param(int i, int j)
{
  return &ia_params[i*n_particle_types + j];
}

IA_parameters *safe_get_ia_param(int i, int j)
{
  if ((i < 0) || (j < 0))
    return NULL;

  /* expand array if necessary */
  realloc_ia_params(((i > j) ? i : j) + 1);

  return &ia_params[i*n_particle_types + j];
}

void realloc_ia_params(int nsize)
{
  int i, j;
  IA_parameters *new_params;

  if (nsize <= n_particle_types)
    return;

  if (ia_params)
    free(ia_params);
  new_params = (IA_parameters *) malloc(nsize*nsize*sizeof(IA_parameters));
  for (i = 0; i < nsize; i++)
    for (j = 0; j < nsize; j++) {
      if ((i < n_particle_types) && (i < n_particle_types))
	copy_ia_params(&new_params[i*n_particle_types + j],
			   &ia_params[i*n_particle_types + j]);
      else
	initialize_ia_params(&new_params[i*n_particle_types + j]);
    }

  n_particle_types = nsize;
  ia_params = new_params;
}

void copy_ia_params(IA_parameters *dst, IA_parameters *src)
{
  /* change if any allocating ia types enter !!!! */
  memcpy(dst, src, sizeof(IA_parameters));
}

void initialize_ia_params(IA_parameters *params)
{
  /* change if an interaction needs non-0 param to do nothing !!!!! */
  memset(params, sizeof(IA_parameters), 0);
}
