#include "particle_data.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include "debug.h"

/************************************************
 * defines
 ************************************************/

/* increment size of particle buffer */
#define PART_INCREMENT 256

/* decrement size of bonded ia
 * if difference between old length
 * (from recycling) and new length is
 * larger than this, reduce to new size */
#define BONDED_REDUCE 64

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
  particles[index].type = 0;
  particles[index].q    = 0;
  particles[index].p[0] = 0;
  particles[index].p[1] = 0;
  particles[index].p[2] = 0;
  particles[index].v[0] = 0;
  particles[index].v[1] = 0;
  particles[index].v[2] = 0;
  particles[index].f[0] = 0;
  particles[index].f[1] = 0;
  particles[index].f[2] = 0;
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
  int tmp;
  for(i=0;i<3;i++) {
    image_box[i] += (tmp = (int)floor(pos[i]/box_l[i]));
    /*pos[i]       = pos[i] - image_box[i]*box_l[i];    */
    pos[i]       = pos[i] - tmp*box_l[i];    
    if(pos[i] < 0. || pos[i] > box_l[i])
      {
	fprintf(stderr,"Warning fold_particle: Particle out of range image_box[%d] = %d\n",i,image_box[i]);
	fprintf(stderr,"exiting!\n");
	exit(1);
      }
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

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;
  int node, j;
  char buffer[50 + TCL_DOUBLE_SPACE];

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <part num> ?what? ?value?\"", (char *) NULL);
    return (TCL_ERROR);
  }
  if (!processor_grid_is_set())
    setup_processor_grid();

  part_num = atol(argv[1]);
  if (part_num < 0) {
    Tcl_AppendResult(interp, "illegal particle", (char *) NULL);
    return (TCL_ERROR);
  }

  if (!particle_node)
    build_particle_node();

  node = (part_num < n_total_particles) ? particle_node[part_num] : -1;
 
  /* print out particle information */
  if (argc == 2) {
    Particle part;
    /* retrieve particle data */
    if (node == -1) {
      sprintf(buffer, "particle %d does not exist", part_num);
      Tcl_AppendResult(interp, buffer, (char *) NULL);
      return (TCL_ERROR);
    }
    mpi_recv_part(node, part_num, &part);
    Tcl_PrintDouble(interp, part.p[0], buffer);
    Tcl_AppendResult(interp, "p ", buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.p[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.p[2], buffer);
    Tcl_AppendResult(interp, buffer, " type ", (char *)NULL);
    sprintf(buffer, "%d", part.type);
    Tcl_AppendResult(interp, buffer, " q ", (char *)NULL);
    Tcl_PrintDouble(interp, part.q, buffer);
    Tcl_AppendResult(interp, buffer, " v ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[0], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[2], buffer);
    Tcl_AppendResult(interp, buffer, " f ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[0], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    /* FIXME: print bonding structure here */
    free(part.bonds);
    return (TCL_OK);
  }
  
  /* set particle data */
  argc -= 2;
  argv += 2;
  while (argc > 0) {
    if (!strncmp(argv[0], "pos", strlen(argv[0]))) {
      double pos[3];
      if (argc < 4) {
	Tcl_AppendResult(interp, "pos requires 3 arguments", (char *) NULL);
	return (TCL_ERROR);
      }
      /* set position */
      for (j = 0; j < 3; j++) {
	if (Tcl_GetDouble(interp, argv[1 + j], &pos[j]) == TCL_ERROR)
	  return (TCL_ERROR);
      }
    
      if (node == -1) {
	/* spatial position */
	node = find_node(pos);
	mpi_attach_particle(part_num, node);
	map_particle_node(part_num, node);
      }

      mpi_send_pos(node, part_num, pos);

      argc -= 4;
      argv += 4;
    }
    else {
      if (node == -1) {
	Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
	return (TCL_ERROR);
      }
      
      if (!strncmp(argv[0], "q", strlen(argv[0]))) {
	double q;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "q requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set charge */
	if (Tcl_GetDouble(interp, argv[1], &q) == TCL_ERROR)
	  return (TCL_ERROR);
	
	mpi_send_q(node, part_num, q);

	argc -= 2;
	argv += 2;
      }
      else if (!strncmp(argv[0], "v", strlen(argv[0]))) {
	double v[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "v requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &v[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &v[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &v[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_v(node, part_num, v);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "f", strlen(argv[0]))) {
	double f[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "f requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &f[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &f[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &f[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_f(node, part_num, f);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "type", strlen(argv[0]))) {
	int type;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "type requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set type */
	if (Tcl_GetInt(interp, argv[1], &type) == TCL_ERROR)
	  return (TCL_ERROR);

	if (type < 0) {
	  Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
	  return (TCL_ERROR);	  
	} 

	/* make sure type exists */
	realloc_ia_params(type);

	mpi_send_type(node, part_num, type);

	argc -= 2;
	argv += 2;
      }
      else {
	Tcl_AppendResult(interp, "unknown particle parameter \"", argv[0],"\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }

  return (TCL_OK);
}
