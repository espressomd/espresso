// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file particle_data.c
    This file contains everything related to particle storage. If you want to add a new
    property to the particles, it is probably a good idea to modify \ref Particle to give
    scripts access to that property. You always have to modify two positions: first the
    print section, where you should add your new data at the end, and second the read
    section where you have to find a nice and short name for your property to appear in
    the Tcl code. Then you just parse your part out of argc and argv.

    The corresponding header file is particle_data.h.
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "utils.h"
#include "particle_data.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include "interaction_data.h"
#include "integrate.h"
#include "cells.h"
#include "parser.h"
#include "rotation.h"


/************************************************
 * defines
 ************************************************/

/** granularity of the particle buffers in particles */
#define PART_INCREMENT 8

/** my magic MPI code for send/recv_particles */
#define REQ_SNDRCV_PART 0xaa

/************************************************
 * variables
 ************************************************/

int max_seen_particle = -1;
int n_total_particles = 0;
int max_particle_node = 0;
int *particle_node = NULL;
int max_local_particles = 0;
Particle **local_particles = NULL;
Particle *partCfg = NULL;
int partCfgSorted = 0;

/** bondlist for partCfg, if bonds are needed */
IntList partCfg_bl = { NULL, 0, 0 };

/************************************************
 * local functions
 ************************************************/

/** Remove bond from particle if possible */
int try_delete_bond(Particle *part, int *bond);

/** Remove exclusion from particle if possible */
void try_delete_exclusion(Particle *part, int part2);

/** Insert an exclusion if not already set */
void try_add_exclusion(Particle *part, int part2);

/** Automatically add the next \<distance\> neighbors in each molecule to the exclusion list.
    This uses the bond topology obtained directly from the particles, since only this contains
    the full topology, in contrast to \ref topology::topology. To easily setup the bonds, all data
    should be on a single node, therefore the \ref partCfg array is used. With large amounts
    of particles, you should avoid this function and setup exclusions manually. */
void auto_exclusion(int distance);

/************************************************
 * particle initialization functions
 ************************************************/

void init_particle(Particle *part) 
{
  /* ParticleProperties */
  part->p.identity = -1;
  part->p.type     = 0;
  part->p.mol_id   = -1;

#ifdef MASS
  part->p.mass     = 1.0;
#endif

#ifdef ELECTROSTATICS
  part->p.q        = 0.0;
#endif

  /* ParticlePosition */
  part->r.p[0]     = 0.0;
  part->r.p[1]     = 0.0;
  part->r.p[2]     = 0.0;

#ifdef BOND_CONSTRAINT
  part->r.p_old[0] = 0.0;
  part->r.p_old[1] = 0.0;
  part->r.p_old[2] = 0.0;
#endif

#ifdef ROTATION
  part->r.quat[0]  = 1.0;
  part->r.quat[1]  = 0.0;
  part->r.quat[2]  = 0.0;
  part->r.quat[3]  = 0.0;
#endif

#ifdef DIPOLES
  part->r.dip[0]    = 0.0;
  part->r.dip[1]    = 0.0;
  part->r.dip[2]    = 0.0;
  part->p.dipm      = 0.0;
#endif

  /* ParticleMomentum */
  part->m.v[0]     = 0.0;
  part->m.v[1]     = 0.0;
  part->m.v[2]     = 0.0;
#ifdef ROTATION
  part->m.omega[0] = 0.0;
  part->m.omega[1] = 0.0;
  part->m.omega[2] = 0.0;
#endif

  /* ParticleForce */
  part->f.f[0]     = 0.0;
  part->f.f[1]     = 0.0;
  part->f.f[2]     = 0.0;
#ifdef ROTATION
  part->f.torque[0] = 0.0;
  part->f.torque[1] = 0.0;
  part->f.torque[2] = 0.0;
#endif

  /* ParticleLocal */
  part->l.p_old[0]   = 0.0;
  part->l.p_old[1]   = 0.0;
  part->l.p_old[2]   = 0.0;
  part->l.i[0]       = 0;
  part->l.i[1]       = 0;
  part->l.i[2]       = 0;

#ifdef EXTERNAL_FORCES
  part->l.ext_flag   = 0;
  part->l.ext_force[0] = 0.0;
  part->l.ext_force[1] = 0.0;
  part->l.ext_force[2] = 0.0;
#endif

  init_intlist(&(part->bl));
#ifdef EXCLUSIONS
  init_intlist(&(part->el));
#endif

#ifdef VIRTUAL_SITES
  part->p.isVirtual      = 0;
#endif
}

void free_particle(Particle *part) {
  realloc_intlist(&(part->bl), 0);
#ifdef EXCLUSIONS
  realloc_intlist(&(part->el), 0);
#endif
}


/************************************************
 * organizational functions
 ************************************************/

void updatePartCfg(int bonds_flag)
{
  int j;

  if(partCfg)
    return;

  partCfg = malloc(n_total_particles*sizeof(Particle));
  if (bonds_flag != WITH_BONDS)
    mpi_get_particles(partCfg, NULL);
  else
    mpi_get_particles(partCfg,&partCfg_bl);
 
  for(j=0; j<n_total_particles; j++)
    unfold_position(partCfg[j].r.p,partCfg[j].l.i);

  partCfgSorted = 0;
}

int sortPartCfg()
{
  int i;
  Particle *sorted;

  if (!partCfg)
    updatePartCfg(WITHOUT_BONDS);

  if (partCfgSorted)
    return 1;

  if (n_total_particles != max_seen_particle + 1)
    return 0;

  sorted = malloc(n_total_particles*sizeof(Particle));
  for(i = 0; i < n_total_particles; i++)
    memcpy(&sorted[partCfg[i].p.identity], &partCfg[i], sizeof(Particle));
  free(partCfg);
  partCfg = sorted;

  partCfgSorted = 1;

  return 1;
}

void freePartCfg()
{
  free(partCfg);
  partCfg = NULL;
  realloc_intlist(&partCfg_bl, 0);
}

/** resize \ref local_particles.
    \param part the highest existing particle
*/
void realloc_local_particles(int part)
{
  if (part >= max_local_particles) {
    /* round up part + 1 in granularity PART_INCREMENT */
    max_local_particles = PART_INCREMENT*((part + PART_INCREMENT)/PART_INCREMENT);
    local_particles = (Particle **)realloc(local_particles, sizeof(Particle *)*max_local_particles);
  }
}

/** resize \ref particle_node.
    This procedure is only used on the master node in Tcl mode.
    \param part the highest existing particle
*/
static void realloc_particle_node(int part)
{
  if (part >= max_particle_node) {
    /* round up part + 1 in granularity PART_INCREMENT */
    max_particle_node = PART_INCREMENT*((part + PART_INCREMENT)/PART_INCREMENT);
    particle_node = (int *)realloc(particle_node, sizeof(int)*max_particle_node);
  }
}

void particle_invalidate_part_node()
{
  /* invalidate particle->node data */
  if (particle_node) {
    free(particle_node);
    particle_node = NULL;
    max_particle_node = 0;
  }
}

void build_particle_node()
{
  realloc_particle_node(max_seen_particle);
  mpi_who_has();
}

void init_particlelist(ParticleList *pList)
{
  pList->n    = 0;
  pList->max  = 0;
  pList->part = NULL;
}

int realloc_particlelist(ParticleList *l, int size)
{
  int old_max = l->max;
  Particle *old_start = l->part;

  PART_TRACE(fprintf(stderr, "%d: realloc_particlelist %p: %d/%d->%d\n", this_node,
		     l, l->n, l->max, size));

  if (size < l->max) {
    if (size == 0)
      /* to be able to free an array again */
      l->max = 0;
    else
      /* shrink not as fast, just lose half, rounded up */
      l->max = PART_INCREMENT*(((l->max + size + 1)/2 +
				PART_INCREMENT - 1)/PART_INCREMENT);
  }
  else
    /* round up */
    l->max = PART_INCREMENT*((size + PART_INCREMENT - 1)/PART_INCREMENT);
  if (l->max != old_max)
    l->part = (Particle *) realloc(l->part, sizeof(Particle)*l->max);
  return l->part != old_start;
}

void update_local_particles(ParticleList *pl)
{
  Particle *p = pl->part;
  int n = pl->n, i;
  for (i = 0; i < n; i++)
    local_particles[p[i].p.identity] = &p[i];
}

Particle *got_particle(ParticleList *l, int id)
{
  int i;

  for (i = 0; i < l->n; i++)
    if (l->part[i].p.identity == id)
      break;
  if (i == l->n)
    return NULL;
  return &(l->part[i]);
}

Particle *append_unindexed_particle(ParticleList *l, Particle *part)
{
  Particle *p;

  realloc_particlelist(l, ++l->n);
  p = &l->part[l->n - 1];

  memcpy(p, part, sizeof(Particle));
  return p;
}

Particle *append_indexed_particle(ParticleList *l, Particle *part)
{
  int re;
  Particle *p;

  re = realloc_particlelist(l, ++l->n);
  p  = &l->part[l->n - 1];

  memcpy(p, part, sizeof(Particle));

  if (re)
    update_local_particles(l);
  else
    local_particles[p->p.identity] = p;
  return p;
}

Particle *move_unindexed_particle(ParticleList *dl, ParticleList *sl, int i)
{
  Particle *dst, *src, *end;

  realloc_particlelist(dl, ++dl->n);
  dst = &dl->part[dl->n - 1];
  src = &sl->part[i];
  end = &sl->part[sl->n - 1];
  memcpy(dst, src, sizeof(Particle));
  if ( src != end )
    memcpy(src, end, sizeof(Particle));
  sl->n -= 1;
  realloc_particlelist(sl, sl->n);
  return dst;
}

Particle *move_indexed_particle(ParticleList *dl, ParticleList *sl, int i)
{
  int re = realloc_particlelist(dl, ++dl->n);
  Particle *dst = &dl->part[dl->n - 1];
  Particle *src = &sl->part[i];
  Particle *end = &sl->part[sl->n - 1];

  memcpy(dst, src, sizeof(Particle));
  if (re) {
    //fprintf(stderr, "%d: m_i_p: update destination list after realloc\n",this_node);
    update_local_particles(dl); }
  else {
    //fprintf(stderr, "%d: m_i_p: update loc_part entry for moved particle (id %d)\n",this_node,dst->p.identity);
    local_particles[dst->p.identity] = dst;
  }
  if ( src != end ) {
    //fprintf(stderr, "%d: m_i_p: copy end particle in source list (id %d)\n",this_node,end->p.identity);
    memcpy(src, end, sizeof(Particle));

  }
  if (realloc_particlelist(sl, --sl->n)) {
    //fprintf(stderr, "%d: m_i_p: update source list after realloc\n",this_node);
    update_local_particles(sl); }
  else if ( src != end ) {
    //fprintf(stderr, "%d: m_i_p: update loc_part entry for end particle (id %d)\n",this_node,src->p.identity);
    local_particles[src->p.identity] = src; }
  return dst;
}

#ifdef ROTATION
void part_print_omega(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->m.omega[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.omega[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.omega[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_torque(Particle *part, char *buffer, Tcl_Interp *interp)
{
  double torque[3];
//in Espresso torques are in body-fixed frames. We should convert they to the space-fixed coordinates.
convert_torques_body_to_space(part, torque);

  Tcl_PrintDouble(interp, torque[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, torque[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, torque[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_quat(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->r.quat[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.quat[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.quat[2], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.quat[3], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

#ifdef DIPOLES
void part_print_dip(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->r.dip[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.dip[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.dip[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

#ifdef VIRTUAL_SITES
void part_print_isVirtual(Particle *part, char *buffer, Tcl_Interp *interp)
{
  sprintf(buffer,"%i", part->p.isVirtual);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
}
#endif

void part_print_v(Particle *part, char *buffer, Tcl_Interp *interp)
{
  /* unscale velocities ! */
  Tcl_PrintDouble(interp, part->m.v[0]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.v[1]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.v[2]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_f(Particle *part, char *buffer, Tcl_Interp *interp)
{
  /* unscale forces ! */
  Tcl_PrintDouble(interp, part->f.f[0]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.f[1]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.f[2]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_position(Particle *part, char *buffer, Tcl_Interp *interp)
{
  double ppos[3];
  int img[3];
  memcpy(ppos, part->r.p, 3*sizeof(double));
  memcpy(img, part->l.i, 3*sizeof(int));
  unfold_position(ppos, img);
  Tcl_PrintDouble(interp, ppos[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_folded_position(Particle *part, char *buffer, Tcl_Interp *interp)
{
  double ppos[3];
  int img[3];
  memcpy(ppos, part->r.p, 3*sizeof(double));
  memcpy(img, part->l.i, 3*sizeof(int));
  fold_position(ppos, img);

  Tcl_PrintDouble(interp, ppos[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void part_print_bonding_structure(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int i,j,size;
  Tcl_AppendResult(interp, "{ ", (char *)NULL);
  i = 0;
  while(i<part->bl.n) {
    size = bonded_ia_params[part->bl.e[i]].num;
    sprintf(buffer, "{%d ", part->bl.e[i]); i++;
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    for(j=0;j<size-1;j++) {
      sprintf(buffer, "%d ", part->bl.e[i]); i++;
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    sprintf(buffer, "%d} ", part->bl.e[i]); i++;
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
}

/* keep a unique list for particle i. Particle j is only added if it is not i
   and not already in the list. */
MDINLINE void add_partner(IntList *il, int i, int j, int distance)
{
  int k;
  if (j == i) return;
  for (k = 0; k < il->n; k += 2)
    if (il->e[k] == j)
      return;
  realloc_intlist(il, il->n + 2);
  il->e[il->n++] = j;
  il->e[il->n++] = distance;
}

/* Add a link of size size consisting of (link[l], p) */
MDINLINE void add_link(IntList *il, IntList *link, int l, int p, int size)
{
  int i;
  realloc_intlist(il, il->n + size);
  for ( i = 0; i < size-1; i++ ) {
    il->e[il->n++] = link->e[l+i];
  } 
  il->e[il->n++] = p;
}
 
/** Return all bond partners of a particle including bonds that are not stored at the particle itself up to a certain distance in numbers of bonds. Return a tcl list to the interpreter c*/
void part_print_bond_partners(Particle *part, char *buffer, Tcl_Interp *interp, int distance)
{
  int c, p, i, p1, p2;
  Bonded_ia_parameters *ia_params;
  Particle *part1;
  /* partners is a list containing the currently found bond partners and their distance in number of 
     bonds for each particle */
  IntList *partners;
  /* link list for all linear connections up to a certain number of bonds for the particle requiested */
  IntList *links;

  /* setup bond partners and distance list. Since we need to identify particles via their identity,
     we use a full sized array */
  partners    = malloc((max_seen_particle + 1)*sizeof(IntList));
  for (p = 0; p <= max_seen_particle; p++) init_intlist(&partners[p]);
  updatePartCfg(WITH_BONDS);

  /* determine initial connectivity */
  for (p = 0; p < n_total_particles; p++) {
    part1 = &partCfg[p];
    p1    = part1->p.identity;
    for (i = 0; i < part1->bl.n;) {
      ia_params = &bonded_ia_params[part1->bl.e[i++]];
      if (ia_params->num == 1) {
	p2 = part1->bl.e[i++];
	add_partner(&partners[p1], p1, p2, 1);
	add_partner(&partners[p2], p2, p1, 1);
      }
      else
	i += ia_params->num;
    }
  }

  /* Create links to particle */
  distance++;
  links    = malloc((distance+1)*sizeof(IntList));
  for( c = 0; c <= distance; c++)  init_intlist(&links[c]);

  p1 = part->p.identity;
  add_link(&links[0], &links[0],0, p1, 1);

  for (p = 0; p < partners[p1].n; p+=2) {
     add_link(&links[1], &links[0], 0, partners[p1].e[p], 2);
  }
 
  for (c = 2; c <= distance; c++) {
    for (i = 0; i < links[c-1].n; i+=c) {
      p2 = links[c-1].e[c+i-1];
      for (p = 0; p < partners[p2].n; p+=2) {
	if( partners[p2].e[p] !=  links[c-1].e[c+i-2]) {
	  add_link(&links[c], &links[c-1], i, partners[p2].e[p],(c+1));
	}
      }
    }
  }

  p1 = part->p.identity;

  for (c = 2; c <= distance; c++) {
    Tcl_AppendResult(interp, " {", (char *)NULL);
    for (i = 0; i < links[c-1].n; i+=c ) {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      for (p = 1; p < c; p++) {
	sprintf(buffer, "%d ",links[c-1].e[i+p]);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
    Tcl_AppendResult(interp, " }", (char *)NULL);
  }

  /* free memory */
  for (p = 0; p <= max_seen_particle; p++) realloc_intlist(&partners[p], 0);
  free(partners);
  for(i=0;i<distance+1; i++) realloc_intlist(&links[i], 0);
  free(links);
}

#ifdef EXCLUSIONS
void part_print_exclusions(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int i;
  for (i = 0; i < part->el.n; i++) {
    sprintf(buffer, "%d ", part->el.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
}
#endif

#ifdef EXTERNAL_FORCES
void part_print_fix(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int i;
  for (i = 0; i < 3; i++) {
    if (part->l.ext_flag & COORD_FIXED(i))
      Tcl_AppendResult(interp, "1 ", (char *)NULL);
    else
	    Tcl_AppendResult(interp, "0 ", (char *)NULL);
  }
}

void part_print_ext_force(Particle *part, char *buffer, Tcl_Interp *interp)
{
  if(part->l.ext_flag & PARTICLE_EXT_FORCE) {
    Tcl_PrintDouble(interp, part->l.ext_force[0], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	  Tcl_PrintDouble(interp, part->l.ext_force[1], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	  Tcl_PrintDouble(interp, part->l.ext_force[2], buffer);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  else {
    Tcl_AppendResult(interp, "0.0 0.0 0.0 ", (char *)NULL);
  }
}
#endif

/** append particle data in ASCII form to the Tcl result.
    @param part_num the particle which data is appended
    @param interp   the Tcl interpreter to which result to add to */
int printParticleToResult(Tcl_Interp *interp, int part_num)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle part;

  if (get_particle_data(part_num, &part) == TCL_ERROR)
    return (TCL_ERROR);

  sprintf(buffer, "%d", part.p.identity);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_AppendResult(interp, " pos ", (char *)NULL);
  part_print_position(&part, buffer, interp);
  Tcl_AppendResult(interp, " type ", (char *)NULL);
  sprintf(buffer, "%d", part.p.type);
  if (part.p.mol_id > -1) {
	Tcl_AppendResult(interp, buffer, " molecule ", (char *)NULL);
  	sprintf(buffer, "%d", part.p.mol_id);
  }
#ifdef MASS
  Tcl_AppendResult(interp, buffer, " mass ", (char *)NULL);
  Tcl_PrintDouble(interp, part.p.mass, buffer);
#endif
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, buffer, " q ", (char *)NULL);
  Tcl_PrintDouble(interp, part.p.q, buffer);
#endif
  Tcl_AppendResult(interp, buffer, " v ", (char *)NULL);
  part_print_v(&part, buffer, interp);
  Tcl_AppendResult(interp, " f ", (char *)NULL);
  part_print_f(&part, buffer, interp);

#ifdef ROTATION
  /* print information about rotation */
  Tcl_AppendResult(interp, " quat ", (char *)NULL);
  part_print_quat(&part, buffer, interp);

  Tcl_AppendResult(interp, " omega ", (char *)NULL);
  part_print_omega(&part, buffer, interp);

  Tcl_AppendResult(interp, " torque ", (char *)NULL);
  part_print_torque(&part, buffer, interp);
#endif

#ifdef DIPOLES
  /* print information about dipoles */
  Tcl_AppendResult(interp, " dip ", (char *)NULL);
  part_print_dip(&part, buffer, interp);
  Tcl_AppendResult(interp, " dipm ", (char *)NULL);
  Tcl_PrintDouble(interp, part.p.dipm, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
#endif

#ifdef VIRTUAL_SITES
  /* print information about isVirtual */
  Tcl_AppendResult(interp, " virtual ", (char *)NULL);
  part_print_isVirtual(&part, buffer, interp);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
#endif

#ifdef EXCLUSIONS
  if (part.el.n > 0) {
    Tcl_AppendResult(interp, " exclude ", (char *)NULL);
    part_print_exclusions(&part, buffer, interp);
  }
#endif

#ifdef EXTERNAL_FORCES
  /* print external force information. */
  if (part.l.ext_flag & PARTICLE_EXT_FORCE) {
    Tcl_AppendResult(interp, " ext_force ", (char *)NULL);
    part_print_ext_force(&part, buffer, interp);
  }

  /* print fix information. */
  if (part.l.ext_flag & COORDS_FIX_MASK) {
    Tcl_AppendResult(interp, " fix ", (char *)NULL);
    part_print_fix(&part, buffer, interp);
  }
#endif

  /* print bonding structure */
  if (part.bl.n > 0) {
    Tcl_AppendResult(interp, " bond ", (char *)NULL);
    part_print_bonding_structure(&part, buffer, interp);
  }

  free_particle(&part);
  return (TCL_OK);
}

int part_print_all(Tcl_Interp *interp)
{
  int i = 0, start = 1;

  if (!particle_node)
    build_particle_node();

  PART_TRACE(fprintf(stderr, "max_seen %d\n", max_seen_particle));

  for (i = 0; i <= max_seen_particle ; i++) {

    PART_TRACE(fprintf(stderr, "particle %d\n", i));

    if (particle_node[i] != -1) {
      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);

      printParticleToResult(interp, i);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
  }

  return TCL_OK;
}

int part_parse_print(Tcl_Interp *interp, int argc, char **argv,
		     int part_num)
{

  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle part;

  if (part_num > max_seen_particle) {
    Tcl_AppendResult(interp, "na", (char *)NULL);

    return TCL_OK;
  }

  if (get_particle_data(part_num, &part) == TCL_ERROR) {
    Tcl_AppendResult(interp, "na", (char *)NULL);
    return TCL_OK;
  }

  while (argc > 0) {
    if (ARG0_IS_S("identity")) {
      sprintf(buffer, "%d", part.p.identity);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else if (ARG0_IS_S("position"))
      part_print_position(&part, buffer, interp);
    else if (ARG0_IS_S("force"))
      part_print_f(&part, buffer, interp);
    else if (ARG0_IS_S("folded_position"))
      part_print_folded_position(&part, buffer, interp);
    else if (ARG0_IS_S("type")) {
      sprintf(buffer, "%d", part.p.type);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else if (ARG0_IS_S("molecule_id")) {
      sprintf(buffer, "%d", part.p.mol_id);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#ifdef MASS
    else if (ARG0_IS_S("mass")) {
      Tcl_PrintDouble(interp, part.p.mass, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#endif

#ifdef ELECTROSTATICS
    else if (ARG0_IS_S("q")) {
      Tcl_PrintDouble(interp, part.p.q, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#endif
    else if (ARG0_IS_S("v"))
      part_print_v(&part, buffer, interp);

#ifdef ROTATION
    else if (ARG0_IS_S("quat"))
      part_print_quat(&part, buffer, interp);
    else if (ARG0_IS_S("omega"))
      part_print_omega(&part, buffer, interp);
    else if (ARG0_IS_S("torque"))
      part_print_torque(&part, buffer, interp);
#endif

#ifdef DIPOLES
    else if (ARG0_IS_S("dip"))
      part_print_dip(&part, buffer, interp);
    else if (ARG0_IS_S("dipm"))
      Tcl_PrintDouble(interp, part.p.dipm, buffer);
#endif

#ifdef VIRTUAL_SITES
    else if (ARG0_IS_S("virtual"))
      part_print_isVirtual(&part, buffer, interp);
#endif

#ifdef EXTERNAL_FORCES
    else if (ARG0_IS_S("ext_force"))
      part_print_ext_force(&part,buffer,interp);
    else if (ARG0_IS_S("fix"))
      part_print_fix(&part, buffer, interp);
#endif

#ifdef EXCLUSIONS
    else if (ARG0_IS_S("exclusions"))
      part_print_exclusions(&part, buffer, interp);
#endif

    else if (ARG0_IS_S("bonds"))
      part_print_bonding_structure(&part, buffer, interp);
    else if (ARG0_IS_S("connections")) {
      int distance = 1;
      if (argc ==2) {
	ARG1_IS_I(distance);
	argc--; argv++;
      }
      part_print_bond_partners(&part, buffer, interp, distance);
    }
    else {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "unknown particle data \"", argv[0], "\" requested", (char *)NULL);
      free_particle(&part);
      return TCL_ERROR;
    }
    if (argc > 1)
      Tcl_AppendResult(interp, " ", (char *)NULL);

    argc--;
    argv++;

  }

  free_particle(&part);
  return TCL_OK;
}

int part_cmd_delete(Tcl_Interp *interp, int argc, char **argv,
		   int part_num, int * change)
{
  *change = 0;

  if (remove_particle(part_num) == TCL_ERROR) {
    char buffer[32 + TCL_INTEGER_SPACE];
    sprintf(buffer, "particle %d", part_num);
    Tcl_AppendResult(interp, buffer, " does not exist and cannot be removed",
		     (char *)NULL);

    return TCL_ERROR;
  }
  if (argc > 0) {
    Tcl_AppendResult(interp, "part <id> delete does not expect parameters",
		     (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_pos(Tcl_Interp *interp, int argc, char **argv,
		   int part_num, int * change)
{
  double pos[3];
  int j;

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "pos requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }

  /* set position */
  for (j = 0; j < 3; j++)
    if (! ARG_IS_D(j, pos[j]))
      return TCL_ERROR;

  if (place_particle(part_num, pos) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle could not be set", (char *) NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

#ifdef MASS
int part_parse_mass(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    double mass;

    *change = 1;

    if (argc < 1) {
      Tcl_AppendResult(interp, "mass requires 1 argument", (char *) NULL);
      return TCL_ERROR;
    }

    /* set mass */
    if (! ARG0_IS_D(mass))
      return TCL_ERROR;

    if (set_particle_mass(part_num, mass) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}
#endif

#ifdef DIPOLES
int part_parse_dipm(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    double dipm;

    *change = 1;

    if (argc < 1) {
      Tcl_AppendResult(interp, "dipm requires 1 argument", (char *) NULL);
      return TCL_ERROR;
    }

    /* set dipole moment */
    if (! ARG0_IS_D(dipm))
      return TCL_ERROR;

    if (set_particle_dipm(part_num, dipm) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}

int part_parse_dip(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double dip[3];
  double quat[4];
  double dipm;
  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "dip requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set dipole orientation */
  if (! ARG_IS_D(0, dip[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, dip[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, dip[2]))
    return TCL_ERROR;

  dipm = sqrt(dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2]);

  if (dipm == 0) {
     Tcl_AppendResult(interp, "cannot set dipole with zero length", (char *)NULL);

    return TCL_ERROR;
  }

  convert_dip_to_quat_one(dip, quat);

  if (set_particle_dip(part_num, dip) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

   if (set_particle_dipm(part_num, dipm) == TCL_ERROR) {
     Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
   }

  if (set_particle_quat(part_num, quat) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

#endif

#ifdef VIRTUAL_SITES
int part_parse_isVirtual(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    int isVirtual;

    *change = 1;

    if (argc < 1) {
      Tcl_AppendResult(interp, "virtual requires 1 argument", (char *) NULL);
      return TCL_ERROR;
    }

    /* set isVirtual */
    if (! ARG0_IS_I(isVirtual))
      return TCL_ERROR;

    if (set_particle_isVirtual(part_num, isVirtual) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}
#endif

int part_parse_q(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    double q;

    *change = 1;

    if (argc < 1) {
      Tcl_AppendResult(interp, "q requires 1 argument", (char *) NULL);
      return TCL_ERROR;
    }

    /* set charge */
    if (! ARG0_IS_D(q))
      return TCL_ERROR;

    if (set_particle_q(part_num, q) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}

int part_parse_v(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
  double v[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "v requires 3 arguments", (char *) NULL);

    return TCL_ERROR;
  }

  /* set v */
  if (! ARG_IS_D(0, v[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, v[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, v[2]))
    return TCL_ERROR;

  /* scale velocity with time step */
  v[0] *= time_step;
  v[1] *= time_step;
  v[2] *= time_step;

  if (set_particle_v(part_num, v) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_f(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
  double f[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "f requires 3 arguments", (char *) NULL);

    return TCL_ERROR;
  }

  /* set f */
  if (! ARG_IS_D(0, f[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, f[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, f[2]))
    return TCL_ERROR;

  /* rescale forces */
  f[0] *= (0.5*time_step*time_step);
  f[1] *= (0.5*time_step*time_step);
  f[2] *= (0.5*time_step*time_step);

  if (set_particle_f(part_num, f) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_type(Tcl_Interp *interp, int argc, char **argv,
		    int part_num, int * change)
{
  int type;

  *change = 1;

  if (argc < 1 ) {
    Tcl_AppendResult(interp, "type requires 1 argument", (char *) NULL);
    return TCL_ERROR;
  }

  /* set type */
  if (! ARG0_IS_I(type))
    return TCL_ERROR;

  if (type < 0) {
    Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
    return TCL_ERROR;
  }

  if (set_particle_type(part_num, type) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_mol_id(Tcl_Interp *interp, int argc, char **argv,
		      int part_num, int * change)
{
  int mid;

  *change = 1;

  if (argc < 1 ) {
    Tcl_AppendResult(interp, "molecule_id requires 1 argument", (char *) NULL);
    return TCL_ERROR;
  }

  /* set mid */
  if (ARG0_IS_S("off"))
    mid = -1;
  else {
    /* number >= -1 (off)*/
    if (!ARG0_IS_I(mid) || mid < -1) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "invalid molecule id, must be integer or \"off\"", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if (set_particle_mol_id(part_num, mid) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

#ifdef ROTATION

int part_parse_quat(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double quat[4];

  *change = 4;

#ifdef DIPOLES
 //Here we should check if the dipole moment of the particle is already given
#endif

  if (argc < 4) {
    Tcl_AppendResult(interp, "quaternion requires 4 arguments", (char *) NULL);
    return TCL_ERROR;
  }

  /* set orientation */
  if (! ARG_IS_D(0, quat[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, quat[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, quat[2]))
    return TCL_ERROR;

  if (! ARG_IS_D(3, quat[3]))
    return TCL_ERROR;

  if (set_particle_quat(part_num, quat) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_omega(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double omega[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "omega requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set angular velocity */
  if (! ARG_IS_D(0, omega[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, omega[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, omega[2]))
    return TCL_ERROR;

   if (set_particle_omega(part_num, omega) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_torque(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double torque[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "torque requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set torque */
  if (! ARG_IS_D(0, torque[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, torque[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, torque[2]))
    return TCL_ERROR;

  if (set_particle_torque(part_num, torque) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }


  return TCL_OK;
}
#endif

#ifdef EXTERNAL_FORCES

int part_parse_ext_force(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double ext_f[3];
  int ext_flag;

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "ext_force requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set external force */
  if (! ARG_IS_D(0, ext_f[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, ext_f[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, ext_f[2]))
    return TCL_ERROR;

  if (ext_f[0] == 0 && ext_f[1] == 0 && ext_f[2] == 0)
    ext_flag = 0;
  else
    ext_flag = PARTICLE_EXT_FORCE;

  if (set_particle_ext(part_num, ext_flag, ext_f) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_fix(Tcl_Interp *interp, int argc, char **argv,
		   int part_num, int * change)
{
  int fixed_coord_flag[3] = {0, 0, 0};
  int ext_flag = 0;
  int i;

  if (argc == 0 || !ARG_IS_I(0,fixed_coord_flag[0])) {
    Tcl_ResetResult(interp);

    fixed_coord_flag[0] = 1;
    fixed_coord_flag[1] = 1;
    fixed_coord_flag[2] = 1;
    *change = 0;
  }
  else {
    if (argc < 3){
      Tcl_AppendResult(interp, "fix has either 0 or 3 arguments: { <fixed_coord> <fixed_coord> <fixed_coord> }", (char *)NULL);
      return TCL_ERROR;
    }
    /* first one already read in the if statement above */
    if (! ARG_IS_I(1,fixed_coord_flag[1] ))
      return TCL_ERROR;
    if (! ARG_IS_I(2,fixed_coord_flag[2] ))
      return TCL_ERROR;
    *change = 3;
  }

  for (i = 0; i < 3l; i++)
    if (fixed_coord_flag[i])
      ext_flag |= COORD_FIXED(i);

  if (set_particle_fix(part_num, ext_flag) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_unfix(Tcl_Interp *interp, int argc, char **argv,
		     int part_num, int * change)
{
  *change = 0;

  if (set_particle_fix(part_num, 0) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

#endif

int part_parse_bond(Tcl_Interp *interp, int argc, char **argv,
		    int part_num, int * change)
{
  int delete = 0;
  int type_num;
  int n_partners;
  /* Bond type number and the bond partner atoms are stored in this field. */
  int *bond;
  int j;

  *change = 0;

  /* check number of arguments */
  if (argc < 1) {
    Tcl_AppendResult(interp, "bond requires at least 1 arguments: "
		     "[delete] { <type_num> <partner> }", (char *) NULL);
    return TCL_ERROR;
  }

  /* parse away delete eventually */
  if (ARG0_IS_S("delete")) {
    delete = 1;
    argc--;
    argv++;
    *change += 1;
  }

  if (argc == 1) {
    /* check for the new, nested Tcl-list format */
    int param1, tmp_argc;
    char  **tmp_argv;
    bond = NULL;

    if (!particle_node)
      build_particle_node();
    
    Tcl_SplitList(interp, argv[0], &tmp_argc, &tmp_argv);

    for(param1 = 0 ; param1 < tmp_argc; param1++) {
      int tmp2_argc;
      char  **tmp2_argv;
      Tcl_SplitList(interp, tmp_argv[param1], &tmp2_argc, &tmp2_argv);
      if (tmp2_argc < 1) {
	Tcl_AppendResult(interp, "usage: part <p> bond <type> <partner>+\n", (char *) NULL);
	return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, tmp2_argv[0], &type_num) == TCL_ERROR) return TCL_ERROR;
      if(type_num < 0 || type_num >= n_bonded_ia) {
	Tcl_AppendResult(interp, "invalid bonded interaction type_num"
			 " (set bonded interaction parameters first)", (char *) NULL);
	return TCL_ERROR;
      }

      /* check partners */
      n_partners = bonded_ia_params[type_num].num;
      if(tmp2_argc < 1+n_partners) {
	char buffer[256 + 2*TCL_INTEGER_SPACE];
	sprintf(buffer, "bond type %d requires %d arguments.",
		type_num, n_partners+1);
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	return TCL_ERROR;
      }

      bond = (int *)malloc( (n_partners+1)*sizeof(int) );
      bond[0] = type_num;
      j=1;
      while(j <= n_partners) {
	if (Tcl_GetInt(interp, tmp2_argv[j], &bond[j]) == TCL_ERROR) {
	  free(bond);
	  return TCL_ERROR;
	}
	if(bond[j] < 0 || bond[j] > max_seen_particle || particle_node[bond[j]] == -1) {
	  char buffer[256 + 2*TCL_INTEGER_SPACE];
	  sprintf(buffer, "partner atom %d (identity %d) not known, set all partner atoms first",
		  j+1,bond[j]);
	  Tcl_AppendResult(interp, buffer , (char *) NULL);
	  free(bond);
	  return TCL_ERROR;
	}
	j++;
      }
      /* set/delete bond */
      if (change_particle_bond(part_num, bond, delete) != TCL_OK) {
	Tcl_AppendResult(interp, "bond to delete did not exist", (char *)NULL);
	free(bond);
	return TCL_ERROR;
      }
      Tcl_Free((char *)tmp2_argv);
    }
    free(bond);
    Tcl_Free((char *)tmp_argv);

    *change += 1;
    return TCL_OK;
  }
  /* check for the old, multiple numbers format */
  if (argc == 0) {
    if (delete == 0) {
      Tcl_AppendResult(interp, "usage: part <p> bond <type> <partner>+\n", (char *) NULL);
      return TCL_ERROR;
    }
    // delete all bonds
    bond = NULL;
    // since there is not even a type...
    n_partners = -1;
  }
  else if (! ARG0_IS_I(type_num))
    return TCL_ERROR;
  else {
    if(type_num < 0 || type_num >= n_bonded_ia) {
      Tcl_AppendResult(interp, "invalid bonded interaction type_num"
		       " (set bonded interaction parameters first)", (char *) NULL);
      return TCL_ERROR;
    }
    /* check partners */
    n_partners = bonded_ia_params[type_num].num;
    if(argc < 1+n_partners) {
      char buffer[256 + 2*TCL_INTEGER_SPACE];
      sprintf(buffer, "bond type %d requires %d arguments.",
	      type_num, n_partners+1);
      Tcl_AppendResult(interp, buffer, (char *) NULL);
      return TCL_ERROR;
    }

    if (!particle_node)
      build_particle_node();

    bond = (int *)malloc( (n_partners+1)*sizeof(int) );
    bond[0] = type_num;
    j=1;
    while(j <= n_partners) {
      if (! ARG_IS_I(j, bond[j])) {
	free(bond);
	return TCL_ERROR;
      }
      if(bond[j] < 0 || bond[j] > max_seen_particle || particle_node[bond[j]] == -1) {
	char buffer[256 + 2*TCL_INTEGER_SPACE];
	sprintf(buffer, "partner atom %d (identity %d) not known, set all partner atoms first",
		j+1,bond[j]);
	Tcl_AppendResult(interp, buffer , (char *) NULL);
	free(bond);
	return TCL_ERROR;
      }
      j++;
    }
  }
  /* set/delete bond */
  if (change_particle_bond(part_num, bond, delete) != TCL_OK) {
    Tcl_AppendResult(interp, "bond to delete did not exist", (char *)NULL);
    free(bond);
    return TCL_ERROR;
  }
  free(bond);

  *change += (1 + n_partners);

  return TCL_OK;
}

#ifdef EXCLUSIONS
int part_parse_exclusion(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  int delete = 0, partner;

  *change = 0;

  /* check number of arguments */
  if (argc < 1) {
    Tcl_AppendResult(interp, "exclusion requires at least 1 arguments: "
		     "[delete] { <partner 1> <partner 2> ... }", (char *) NULL);
    return TCL_ERROR;
  }

  /* parse away delete eventually */
  if (ARG0_IS_S("delete")) {
    delete = 1;
    argc--;
    argv++;
    (*change)++;

    if (argc < 1) {
      Tcl_AppendResult(interp, "exclusion requires at least 1 arguments: "
		       "[delete] { <partner 1> <partner 2> ... }", (char *) NULL);
      return TCL_ERROR;
    }
  }

  /* parse partners */
  if (!particle_node)
    build_particle_node();

  while (argc > 0) {
    if (!ARG0_IS_I(partner)) {
      /* seems to be the next property */
      Tcl_ResetResult(interp);
      break;
    }

    /* set/delete exclusion */
    if (change_exclusion(part_num, partner, delete) != TCL_OK) {
      if (delete)
	Tcl_AppendResult(interp, "exclusion to delete did not exist", (char *)NULL);
      else
	Tcl_AppendResult(interp, "particle to exclude from interaction does not exist or is the same", (char *)NULL);	
      return TCL_ERROR;
    }

    argc--; argv++;
    (*change)++;
  }
  return TCL_OK;
}
#endif

int part_parse_cmd(Tcl_Interp *interp, int argc, char **argv,
		   int part_num)
{
  int change = 0, err = TCL_OK;

#ifdef ADDITIONAL_CHECKS
  if (!particle_node)
    build_particle_node();
  mpi_bcast_event(CHECK_PARTICLES);
#endif

  if (ARG0_IS_S("print"))
    err = part_parse_print(interp, argc-1, argv+1, part_num);
  else while (argc > 0) {
    if (ARG0_IS_S("delete"))
      err = part_cmd_delete(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("pos"))
      err = part_parse_pos(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("type"))
      err = part_parse_type(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("molecule_id"))
      err = part_parse_mol_id(interp, argc-1, argv+1, part_num, &change);
#ifdef MASS
    else if (ARG0_IS_S("mass"))
      err = part_parse_mass(interp, argc-1, argv+1, part_num, &change);
#endif
    else if (ARG0_IS_S("q"))
      err = part_parse_q(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("v"))
      err = part_parse_v(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("f"))
      err = part_parse_f(interp, argc-1, argv+1, part_num, &change);

#ifdef ROTATION

    else if (ARG0_IS_S("quat"))
      err = part_parse_quat(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("omega"))
      err = part_parse_omega(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("torque"))
      err = part_parse_torque(interp, argc-1, argv+1, part_num, &change);

#endif

#ifdef DIPOLES

    else if (ARG0_IS_S("dip"))
      err = part_parse_dip(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("dipm"))
      err = part_parse_dipm(interp, argc-1, argv+1, part_num, &change);

#endif

#ifdef VIRTUAL_SITES

    else if (ARG0_IS_S("virtual"))
      err = part_parse_isVirtual(interp, argc-1, argv+1, part_num, &change);

#endif

#ifdef EXTERNAL_FORCES

    else if (ARG0_IS_S("unfix"))
      err = part_parse_unfix(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("fix"))
      err = part_parse_fix(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("ext_force"))
      err = part_parse_ext_force(interp, argc-1, argv+1, part_num, &change);

#endif

    else if (ARG0_IS_S("bond"))
      err = part_parse_bond(interp, argc-1, argv+1, part_num, &change);

#ifdef EXCLUSIONS
    else if (ARG0_IS_S("exclude"))
      err = part_parse_exclusion(interp, argc-1, argv+1, part_num, &change);
#endif

    else {
      Tcl_AppendResult(interp, "unknown particle parameter \"",
		       argv[0],"\"", (char *)NULL);
      err = TCL_ERROR;
    }
    
    if (err == TCL_ERROR)
      break;

    argc -= (change + 1);
    argv += (change + 1);
  }

#ifdef ADDITIONAL_CHECKS
  mpi_bcast_event(CHECK_PARTICLES);
#endif

  return mpi_gather_runtime_errors(interp, err);
}

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;

  /* if no further arguments are given, print out all stored particles */
  if (argc == 1)
    return part_print_all(interp);
  
  /* parse away all non particle number arguments */
  if (ARG1_IS_S("deleteall")) {
    if (argc != 2) {
      Tcl_AppendResult(interp, "part deleteall takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    remove_all_particles();
    return TCL_OK;
  }
#ifdef EXCLUSIONS
  else if (ARG1_IS_S("delete_exclusions")) {
    if (argc != 2) {
      Tcl_AppendResult(interp, "part delete_exclusions takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    remove_all_exclusions();
    return TCL_OK;
  }
  else if (ARG1_IS_S("auto_exclusions")) {
    int bond_neighbors = 2;
    if (argc > 3) {
      Tcl_AppendResult(interp, "usage: part auto_exclusions {distance (bond neighbors to exclude, defaults to 2)}", (char *)NULL);
      return TCL_ERROR;
    }
    if (argc == 3 && !ARG_IS_I(2, bond_neighbors))
      return TCL_ERROR;
    auto_exclusion(bond_neighbors);
    return TCL_OK;
  }
#endif

  /* only one argument is given */
  if (argc == 2) {

    if (ARG1_IS_I(part_num)) {
      if (printParticleToResult(interp, part_num) == TCL_ERROR)
	Tcl_AppendResult(interp, "na", (char *)NULL);
      return TCL_OK;
    }

    return TCL_ERROR;
  }

  /* from here we have at least two arguments */
  /* The first one has to be an integer */
  if (! ARG1_IS_I(part_num))
    return TCL_ERROR;

  if (part_num < 0) {
    Tcl_AppendResult(interp, "particle identities must be positive",
		     (char *)NULL);

    return TCL_ERROR;
  }

  return part_parse_cmd(interp, argc-2, argv+2, part_num);
}

int get_particle_data(int part, Particle *data)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;

  pnode = particle_node[part];
  if (pnode == -1)
    return TCL_ERROR;
  mpi_recv_part(pnode, part, data);
  return TCL_OK;
}

int place_particle(int part, double p[3])
{
  int new, i;
  int pnode, retcode = TCL_OK;

  if (part < 0)
    return TCL_ERROR;

  if (!particle_node)
    build_particle_node();

  pnode = (part <= max_seen_particle) ? particle_node[part] : -1;
  new = (pnode == -1);
  if (new) {
    /* new particle, node by spatial position */
    pnode = cell_structure.position_to_node(p);

    /* master node specific stuff */
    realloc_particle_node(part);
    particle_node[part] = pnode;

    /* fill up possible gap */
    for (i = max_seen_particle + 1; i < part; i++)
      particle_node[i] = -1;

    retcode = TCL_CONTINUE;
  }

  mpi_place_particle(pnode, part, new, p);

  return retcode;
}

int set_particle_v(int part, double v[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_v(pnode, part, v);
  return TCL_OK;
}

int set_particle_f(int part, double F[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_f(pnode, part, F);
  return TCL_OK;
}

#ifdef MASS
int set_particle_mass(int part, double mass)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_mass(pnode, part, mass);
  return TCL_OK;
}
#endif

#ifdef DIPOLES
int set_particle_dipm(int part, double dipm)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_dipm(pnode, part, dipm);
  return TCL_OK;
}

int set_particle_dip(int part, double dip[3])
{
  int pnode;
  double quat[4];

  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_dip(pnode, part, dip);

  convert_dip_to_quat_one(dip, quat);
  mpi_send_quat(pnode, part, quat);
  return TCL_OK;
}

#endif

#ifdef VIRTUAL_SITES
int set_particle_isVirtual(int part, int isVirtual)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_isVirtual(pnode, part, isVirtual);
  return TCL_OK;
}
#endif

int set_particle_q(int part, double q)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_q(pnode, part, q);
  return TCL_OK;
}

int set_particle_type(int part, int type)
{
  int pnode;

  make_particle_type_exist(type);

  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_type(pnode, part, type);
  return TCL_OK;
}

int set_particle_mol_id(int part, int mid)
{
  int pnode;

  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_mol_id(pnode, part, mid);
  return TCL_OK;
}

#ifdef ROTATION
int set_particle_quat(int part, double quat[4])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_quat(pnode, part, quat);
  return TCL_OK;
}

int set_particle_omega(int part, double omega[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_omega(pnode, part, omega);
  return TCL_OK;
}

int set_particle_torque(int part, double torque[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_torque(pnode, part, torque);
  return TCL_OK;
}

#endif

#ifdef EXTERNAL_FORCES
int set_particle_ext(int part, int flag, double force[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_ext(pnode, part, flag, PARTICLE_EXT_FORCE, force);
  return TCL_OK;
}

int set_particle_fix(int part,  int flag)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  mpi_send_ext(pnode, part, flag, COORDS_FIX_MASK, NULL);
  return TCL_OK;
}

#endif

int change_particle_bond(int part, int *bond, int delete)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return TCL_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return TCL_ERROR;
  if(delete != 0 || bond == NULL)
    delete = 1;

  if (bond != NULL) {
    if (bond[0] < 0 || bond[0] >= n_bonded_ia) {
      char *errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{048 invalid/unknown bonded interaction type %d}", bond[0]);
      return TCL_ERROR;
    }
  }
  return mpi_send_bond(pnode, part, bond, delete);
}

void remove_all_particles()
{
  mpi_remove_particle(-1, -1);
  realloc_particle_node(0);
}

int remove_particle(int part)
{
  int pnode;

  if (!particle_node)
    build_particle_node();

  if (part > max_seen_particle)
    return TCL_ERROR;

  pnode = particle_node[part];
  if (pnode == -1)
    return TCL_ERROR;

  particle_node[part] = -1;

  mpi_remove_particle(pnode, part);

  if (part == max_seen_particle) {
    while (max_seen_particle >= 0 && particle_node[max_seen_particle] == -1)
      max_seen_particle--;
    mpi_bcast_parameter(FIELD_MAXPART);
  }

  return TCL_OK;
}

void local_remove_particle(int part)
{
  int ind, c;
  Particle *p = local_particles[part];
  ParticleList *pl = NULL, *tmp;

  /* the tricky - say ugly - part: determine
     the cell the particle is located in by checking
     wether the particle address is inside the array */
  for (c = 0; c < local_cells.n; c++) {
    tmp = local_cells.cell[c];
    ind = p - tmp->part;
    if (ind >= 0 && ind < tmp->n) {
      pl = tmp;
      break;
    }
  }
  if (!pl) {
    fprintf(stderr, "%d: INTERNAL ERROR: could not find cell of particle %d, exiting\n",
	    this_node, part);
    errexit();
  }

  free_particle(p);

  /* remove local_particles entry */
  local_particles[p->p.identity] = NULL;

  if (&pl->part[pl->n - 1] != p) {
    /* move last particle to free position */
    memcpy(p, &pl->part[pl->n - 1], sizeof(Particle));
    /* update the local_particles array for the moved particle */
    local_particles[p->p.identity] = p;
  }

  pl->n--;
}

void local_place_particle(int part, double p[3], int new)
{
  Cell *cell;
  double pp[3];
  int i[3], rl;
  Particle *pt;

  i[0] = 0;
  i[1] = 0;
  i[2] = 0;
  pp[0] = p[0];
  pp[1] = p[1];
  pp[2] = p[2];
  fold_position(pp, i);
  
  if (new) {
    /* allocate particle anew */
    cell = cell_structure.position_to_cell(pp);
    if (!cell) {
      fprintf(stderr, "%d: INTERNAL ERROR: particle %d at %f(%f) %f(%f) %f(%f) does not belong on this node\n",
	      this_node, part, p[0], pp[0], p[1], pp[1], p[2], pp[2]);
      errexit();
    }
    rl = realloc_particlelist(cell, ++cell->n);
    pt = &cell->part[cell->n - 1];
    init_particle(pt);

    pt->p.identity = part;
    if (rl)
      update_local_particles(cell);
    else
      local_particles[pt->p.identity] = pt;
  }
  else
    pt = local_particles[part];

  PART_TRACE(fprintf(stderr, "%d: local_place_particle: got particle id=%d @ %f %f %f\n",
		     this_node, part, p[0], p[1], p[2]));

  memcpy(pt->r.p, pp, 3*sizeof(double));
  memcpy(pt->l.i, i, 3*sizeof(int));
#ifdef BOND_CONSTRAINT
  memcpy(pt->r.p_old, pp, 3*sizeof(double));
#endif
}

void local_remove_all_particles()
{
  Cell *cell;
  int c;
  n_total_particles = 0;
  max_seen_particle = -1;
  for (c = 0; c < local_cells.n; c++) {
    Particle *p;
    int i,   np;
    cell = local_cells.cell[c];
    p = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      realloc_intlist(&p[i].bl, 0);
    cell->n = 0;
  }
}

void local_rescale_particles(int dir, double scale) {
  Particle *p,*p1;
  int j, c;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    for(j = 0; j < cell->n; j++) {
      p1 = &p[j];
      if(dir < 3)
	p1->r.p[dir] *= scale;
      else {
	p1->r.p[0] *= scale;
	p1->r.p[1] *= scale;
	p1->r.p[2] *= scale;
      }
    }
  }
}


void added_particle(int part)
{
  int i;

  n_total_particles++;

  if (part > max_seen_particle) {
    realloc_local_particles(part);
    /* fill up possible gap. Part itself is ESSENTIAL!!!  */
    for (i = max_seen_particle + 1; i <= part; i++)
      local_particles[i] = NULL;
    max_seen_particle = part;
  }
}

int local_change_bond(int part, int *bond, int delete)
{
  IntList *bl;
  Particle *p;
  int bond_size;
  int i;

  p = local_particles[part];
  if (delete)
    return try_delete_bond(p, bond);

  bond_size = bonded_ia_params[bond[0]].num + 1;
  bl = &(p->bl);
  realloc_intlist(bl, bl->n + bond_size);
  for(i = 0; i < bond_size; i++)
    bl->e[bl->n++] = bond[i];
  return TCL_OK;
}

int try_delete_bond(Particle *part, int *bond)
{
  IntList *bl = &part->bl;
  int i, j, type, partners;

  if (!bond) {
    realloc_intlist(bl, bl->n = 0);
    return TCL_OK;
  }

  for (i = 0; i < bl->n;) {
    type = bl->e[i];
    partners = bonded_ia_params[type].num;
    if (type != bond[0])
      i += 1 + partners;
    else {
      for(j = 1; j <= partners; j++)
      { 
	if (bond[j] != bl->e[i + j])
	  break;
      }
      if (j > partners) {
	bl->n -= 1 + partners;
	memcpy(bl->e + i, bl->e + i + 1 + partners, sizeof(int)*(bl->n - i));
	realloc_intlist(bl, bl->n);
	return TCL_OK;
      }
      i += 1 + partners;
    }
  }
  return TCL_ERROR;
}

void remove_all_bonds_to(int identity)
{
  Cell *cell;
  int p, np, c;
  Particle *part;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    np = cell->n;
    part = cell->part;
    for (p = 0; p < np; p++) {
      IntList *bl = &part[p].bl;
      int i, j, partners;

      for (i = 0; i < bl->n;) {
	partners = bonded_ia_params[bl->e[i]].num;
	for(j = 1; j <= partners; j++)
	  if (bl->e[i + j] == identity)
	    break;
	if (j <= partners) {
	  bl->n -= 1 + partners;
	  memcpy(bl->e + i, bl->e + i + 1 + partners,
		 sizeof(int)*(bl->n - i));
	  realloc_intlist(bl, bl->n);
	}
	else
          i += 1 + partners;
      }
      if (i != bl->n) {
	fprintf(stderr, "%d: INTERNAL ERROR: bond information corrupt for particle %d, exiting...\n",
		this_node, part[p].p.identity);
	errexit();
      }
    }
  }
}

#ifdef EXCLUSIONS
void local_change_exclusion(int part1, int part2, int delete)
{
  Cell *cell;
  int p, np, c;
  Particle *part;

  if (part1 == -1 && part2 == -1) {
    /* delete all exclusions */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      np = cell->n;
      part = cell->part;
      for (p = 0; p < np; p++)
	realloc_intlist(&part[p].el, part[p].el.n = 0);
    }
    return;
  }

  /* part1, if here */
  part = local_particles[part1];
  if (part) {
    if (delete)
      try_delete_exclusion(part, part2);
    else
      try_add_exclusion(part, part2);
  }

  /* part2, if here */
  part = local_particles[part2];
  if (part) {
    if (delete)
      try_delete_exclusion(part, part1);
    else
      try_add_exclusion(part, part1);
  }
}

void try_add_exclusion(Particle *part, int part2)
{
  int i;
  for (i = 0; i < part->el.n; i++)
    if (part->el.e[i] == part2)
      return;
  
  realloc_intlist(&part->el, part->el.n + 1);
  part->el.e[part->el.n++] = part2;
}

void try_delete_exclusion(Particle *part, int part2)
{
  IntList *el = &part->el;
  int i;

  for (i = 0; i < el->n;) {
    if (el->e[i] == part2) {
      el->n--;
      memcpy(el->e + i, el->e + i + 1, sizeof(int)*(el->n - i));
      realloc_intlist(el, el->n);
      break;
    }
  }
}
#endif

void send_particles(ParticleList *particles, int node)
{
  int pc;
  /* Dynamic data, bonds and exclusions */
  IntList local_dyn;

  PART_TRACE(fprintf(stderr, "%d: send_particles %d to %d\n", this_node, particles->n, node));

  MPI_Send(&particles->n, 1, MPI_INT, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);
  MPI_Send(particles->part, particles->n*sizeof(Particle),
	   MPI_BYTE, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);

  init_intlist(&local_dyn);
  for (pc = 0; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    int size =  local_dyn.n + p->bl.n;
#ifdef EXCLUSIONS
    size += p->el.n;
#endif
    realloc_intlist(&local_dyn, size);
    memcpy(local_dyn.e + local_dyn.n, p->bl.e, p->bl.n*sizeof(int));
    local_dyn.n += p->bl.n;
#ifdef EXCLUSIONS
    memcpy(local_dyn.e + local_dyn.n, p->el.e, p->el.n*sizeof(int));
    local_dyn.n += p->el.n;
#endif
  }

  PART_TRACE(fprintf(stderr, "%d: send_particles sending %d bond ints\n", this_node, local_dyn.n));
  if (local_dyn.n > 0) {
    MPI_Send(local_dyn.e, local_dyn.n*sizeof(int),
	     MPI_BYTE, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);
    realloc_intlist(&local_dyn, 0);
  }

  /* remove particles from this nodes local list and free data */
  for (pc = 0; pc < particles->n; pc++) {
    local_particles[particles->part[pc].p.identity] = NULL;
    free_particle(&particles->part[pc]);
  }

  realloc_particlelist(particles, particles->n = 0);
}

void recv_particles(ParticleList *particles, int node)
{
  int transfer=0, read, pc;
  MPI_Status status;
  IntList local_dyn;

  PART_TRACE(fprintf(stderr, "%d: recv_particles from %d\n", this_node, node));

  MPI_Recv(&transfer, 1, MPI_INT, node,
	   REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);

  PART_TRACE(fprintf(stderr, "%d: recv_particles get %d\n", this_node, transfer));

  realloc_particlelist(particles, particles->n + transfer);
  MPI_Recv(&particles->part[particles->n], transfer*sizeof(Particle), MPI_BYTE, node,
	   REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);
  particles->n += transfer;

  init_intlist(&local_dyn);
  for (pc = particles->n - transfer; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    local_dyn.n += p->bl.n;
#ifdef EXCLUSIONS
    local_dyn.n += p->el.n;
#endif

    PART_TRACE(fprintf(stderr, "%d: recv_particles got particle %d\n", this_node, p->p.identity));
#ifdef ADDITIONAL_CHECKS
    if (local_particles[p->p.identity] != NULL) {
      fprintf(stderr, "%d: transmitted particle %d is already here...\n", this_node, p->p.identity);
      errexit();
    }
#endif
  }

  update_local_particles(particles);

  PART_TRACE(fprintf(stderr, "%d: recv_particles expecting %d bond ints\n", this_node, local_dyn.n));
  if (local_dyn.n > 0) {
    alloc_intlist(&local_dyn, local_dyn.n);
    MPI_Recv(local_dyn.e, local_dyn.n*sizeof(int), MPI_BYTE, node,
	     REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);
  }
  read = 0;
  for (pc = particles->n - transfer; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    if (p->bl.n > 0) {
      alloc_intlist(&p->bl, p->bl.n);
      memcpy(p->bl.e, &local_dyn.e[read], p->bl.n*sizeof(int));
      read += p->bl.n;
    }
    else
      p->bl.e = NULL;
#ifdef EXCLUSIONS
    if (p->el.n > 0) {
      alloc_intlist(&p->el, p->el.n);
      memcpy(p->el.e, &local_dyn.e[read], p->el.n*sizeof(int));
      read += p->el.n;
    }
    else
      p->el.e = NULL;
#endif
  }
  if (local_dyn.n > 0)
    realloc_intlist(&local_dyn, 0);
}

#ifdef EXCLUSIONS

int change_exclusion(int part1, int part2, int delete)
{
  if (!particle_node)
    build_particle_node();

  if (part1 < 0 || part1 > max_seen_particle ||
      part2 < 0 || part2 > max_seen_particle ||
      part1 == part2 ||
      particle_node[part1] == -1 ||
      particle_node[part2] == -1)
    return TCL_ERROR;

  mpi_send_exclusion(part1, part2, delete);
  return TCL_OK;
}

void remove_all_exclusions()
{
  mpi_send_exclusion(-1, -1, 1);
}

void auto_exclusion(int distance)
{
  int count, p, i, j, p1, p2, p3, dist1, dist2;
  Bonded_ia_parameters *ia_params;
  Particle *part1;
  /* partners is a list containing the currently found excluded particles for each particle,
     and their distance, as a interleaved list */
  IntList *partners;

  updatePartCfg(WITH_BONDS);

  /* setup bond partners and distance list. Since we need to identify particles via their identity,
     we use a full sized array */
  partners    = malloc((max_seen_particle + 1)*sizeof(IntList));
  for (p = 0; p <= max_seen_particle; p++)
    init_intlist(&partners[p]);

  /* determine initial connectivity */
  for (p = 0; p < n_total_particles; p++) {
    part1 = &partCfg[p];
    p1    = part1->p.identity;
    for (i = 0; i < part1->bl.n;) {
      ia_params = &bonded_ia_params[part1->bl.e[i++]];
      if (ia_params->num == 1) {
	p2 = part1->bl.e[i++];
	/* you never know what the user does, may bond a particle to itself...? */
	if (p2 != p1) {
	  add_partner(&partners[p1], p1, p2, 1);
	  add_partner(&partners[p2], p2, p1, 1);
	}
      }
      else
	i += ia_params->num;
    }
  }

  /* calculate transient connectivity. For each of the current neighbors,
     also exclude their close enough neighbors.
  */
  for (count = 1; count < distance; count++) {
    for (p1 = 0; p1 <= max_seen_particle; p1++) {
      for (i = 0; i < partners[p1].n; i += 2) {
	p2 = partners[p1].e[i];
	dist1 = partners[p1].e[i + 1];
	if (dist1 > distance) continue;
	/* loop over all partners of the partner */
	for (j = 0; j < partners[p2].n; j += 2) {
	  p3 = partners[p2].e[j];
	  dist2 = dist1 + partners[p2].e[j + 1];
	  if (dist2 > distance) continue;
	  add_partner(&partners[p1], p1, p3, dist2);
	  add_partner(&partners[p3], p3, p1, dist2);
	}
      }
    }
  }

  /* setup the exclusions and clear the arrays. We do not setup the exclusions up there,
     since on_part_change clears the partCfg, so that we would have to restore it
     continously. Of course this could be optimized by bundling the exclusions, but this
     is only done once and the overhead is as much as for setting the bonds, which
     the user apparently accepted.
  */
  for (p = 0; p <= max_seen_particle; p++) {
    for (j = 0; j < partners[p].n; j++)
      if (p < partners[p].e[j]) change_exclusion(p, partners[p].e[j], 0);
    realloc_intlist(&partners[p], 0);
  }
  free(partners);
}

int do_nonbonded(Particle *p1, Particle *p2)
{
  int i, i2;
  /* check for particle 2 in particle 1's exclusion list. The exclusion list is
     symmetric, so this is sufficient. */
  i2  = p2->p.identity;
  for (i = 0; i < p1->el.n; i++)
    if (i2 == p1->el.e[i]) return 0;
  return 1;
}
#endif
