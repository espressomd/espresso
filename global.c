#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/* increment size of particle buffer */
#define PART_INCREMENT 100

/******************* variables ****************/
double box_l[3]       = {1, 1, 1};
double local_box_l[3] = {1, 1, 1};
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

int n_total_particles = 0;

int     n_particles = 0;
int   max_particles;
Particle *particles;

int *local_index;

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
  if (index < 0) {
    index = n_particles++;
    if (max_particles > n_particles) {
      max_particles += PART_INCREMENT;
      particles = (Particle *)
	realloc(particles, sizeof(Particle)*max_particles);
    }
  }

  particles[index].identity = part;
  particles[index].type = 0;

  particles[index].p[0] =
    particles[index].p[1] =
    particles[index].p[2] = 0;
  particles[index].q = 0;

  particles[index].v[0] =
    particles[index].v[1] =
    particles[index].v[2] = 0;
  particles[index].f[0] =
    particles[index].f[1] =
    particles[index].f[2] = 0;

  particles[index].n_bond = 0;
  particles[index].bonds  = NULL;  
  
  return index;
}
