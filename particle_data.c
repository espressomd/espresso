// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file particle_data.c
    This file contains everything related to particle storage. If you want to add a new
    property to the particles, it is probably a good idea to modify \ref Particle to give
    scripts access to that property. You always have to modify two positions: first the
    print section, where you should add your new data at the end, and second the read
    section where you have to find a nice and short name for your property to appear in
    the Tcl code. Then you just parse your part out of argc and argv.

    The corresponding header file is \ref particle_data.h "particle_data.h".
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "particle_data.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include "interaction_data.h" 
#include "integrate.h"
#include "debug.h"
#include "utils.h"
#include "cells.h"
#include "parser.h"
#include "debug.h"

/* cwz-build-command: make all
 */

/************************************************
 * defines
 ************************************************/

/** granularity of the particle buffers in particles */
#define PART_INCREMENT 32

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

/************************************************
 * particle initialization functions
 ************************************************/

void init_particle(Particle *part) 
{
  /* ParticleProperties */
  part->p.identity = -1;
  part->p.type     = 0;
  part->p.mol_id   = -1;
#ifdef ELECTROSTATICS
  part->p.q        = 0.0;
#endif

  /* ParticlePosition */
  part->r.p[0]     = 0.0;
  part->r.p[1]     = 0.0;
  part->r.p[2]     = 0.0;
#ifdef ROTATION
  part->r.quat[0]  = 1.0;
  part->r.quat[1]  = 0.0;
  part->r.quat[2]  = 0.0;
  part->r.quat[3]  = 0.0;
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
}

void free_particle(Particle *part) {
  realloc_intlist(&(part->bl), 0);
}


/************************************************
 * organizational functions
 ************************************************/

void updatePartCfg(int bonds_flag )
{
  int j;
  IntList bl;

  if(partCfg)
    return;


  partCfg = malloc(n_total_particles*sizeof(Particle));
  if ( bonds_flag != WITH_BONDS ){
    mpi_get_particles(partCfg, NULL);
  }
  else {
    mpi_get_particles(partCfg,&bl);
  }
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

int try_delete_bond(Particle *part, int *bond)
{
  IntList *bl = &part->bl;
  int i, j, type, partners;
  if (!bond) {
    realloc_intlist(bl, bl->n = 0);
    return TCL_OK;
  }

  for (i = 0; i < bl->n;) {
    type = bond[i];
    partners = bonded_ia_params[type].num;
    if (type != bond[0])
      i += 1 + partners;
    else {
      for(j = 1; j <= partners; j++)
	if (bond[j] != bl->e[i + j])
	  break;
      if (j > partners) {
	bl->n -= 1 + partners;
	memcpy(bl->e + i, bl->e + i + 1 + partners,
	       bl->n - i);
	realloc_intlist(bl, bl->n);
	return TCL_OK;
      }
      i += 1 + partners;
    }
  }
  return TCL_ERROR;
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
void part_print_omega(Particle part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part.m.omega[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.m.omega[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.m.omega[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_torque(Particle part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part.f.torque[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f.torque[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f.torque[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_quat(Particle part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part.r.quat[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.quat[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.quat[2], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.quat[3], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}
#endif

void part_print_v(Particle part, char *buffer, Tcl_Interp *interp)
{
  /* unscale velocities ! */
  Tcl_PrintDouble(interp, part.m.v[0]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.m.v[1]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.m.v[2]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_f(Particle part, char *buffer, Tcl_Interp *interp)
{
  /* unscale forces ! */
  Tcl_PrintDouble(interp, part.f.f[0]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f.f[1]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f.f[2]/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_position(Particle part, char *buffer, Tcl_Interp *interp)
{
  double ppos[3];
  int img[3];
  memcpy(ppos, part.r.p, 3*sizeof(double));
  memcpy(img, part.l.i, 3*sizeof(int));
  unfold_position(ppos, img);
  Tcl_PrintDouble(interp, ppos[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ppos[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_folded_position(Particle part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part.r.p[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.p[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.p[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return;
}

void part_print_bonding_structure(Particle part, char *buffer, Tcl_Interp *interp, IntList *bl)
{
  int i=0,j,size;
  Tcl_AppendResult(interp, " { ", (char *)NULL);
  while(i<bl->n) {
    size = bonded_ia_params[bl->e[i]].num;
    sprintf(buffer, "{%d ", bl->e[i]); i++;
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    for(j=0;j<size-1;j++) {
      sprintf(buffer, "%d ", bl->e[i]); i++;
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    sprintf(buffer, "%d} ", bl->e[i]); i++;
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  return;
}

#ifdef EXTERNAL_FORCES
void part_print_fix(Particle part, char *buffer, Tcl_Interp *interp)
{
  int i;
  for (i = 0; i < 3; i++) {
    if (part.l.ext_flag & COORD_FIXED(i))
      Tcl_AppendResult(interp, "1 ", (char *)NULL);
    else
	    Tcl_AppendResult(interp, "0 ", (char *)NULL);
  }
  return;
}

void part_print_ext_force(Particle part, char *buffer, Tcl_Interp *interp)
{
  if(part.l.ext_flag & PARTICLE_EXT_FORCE) {
    Tcl_PrintDouble(interp, part.l.ext_force[0], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	  Tcl_PrintDouble(interp, part.l.ext_force[1], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	  Tcl_PrintDouble(interp, part.l.ext_force[2], buffer);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  else {
    Tcl_AppendResult(interp, "0.0 0.0 0.0 ", (char *)NULL);
  }
  return;
}
#endif

/** append particle data in ASCII form to the Tcl result.
    @param part_num the particle which data is appended
    @param interp   the Tcl interpreter to which result to add to */
int printParticleToResult(Tcl_Interp *interp, int part_num)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle part;
  IntList *bl = &(part.bl);

  if (get_particle_data(part_num, &part) == TCL_ERROR)
    return (TCL_ERROR);

  sprintf(buffer, "%d", part.p.identity);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_AppendResult(interp, " pos ", (char *)NULL);
  part_print_position(part, buffer, interp);
  Tcl_AppendResult(interp, " type ", (char *)NULL);
  sprintf(buffer, "%d", part.p.type);
  Tcl_AppendResult(interp, buffer, " molecule ", (char *)NULL);
  sprintf(buffer, "%d", part.p.mol_id);
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, buffer, " q ", (char *)NULL);
  Tcl_PrintDouble(interp, part.p.q, buffer);
#endif
  Tcl_AppendResult(interp, buffer, " v ", (char *)NULL);
  part_print_v(part, buffer, interp);
  Tcl_AppendResult(interp, " f ", (char *)NULL);
  part_print_f(part, buffer, interp);

#ifdef ROTATION
  /* print information about rotation */
  Tcl_AppendResult(interp, " quat ", (char *)NULL);
  part_print_quat(part, buffer, interp);
             
  Tcl_AppendResult(interp, " omega ", (char *)NULL);
  part_print_omega(part, buffer, interp);
               
  Tcl_AppendResult(interp, " torque ", (char *)NULL);
  part_print_torque(part, buffer, interp);   
#endif      

#ifdef EXTERNAL_FORCES
  /* print external force information. */
  if (part.l.ext_flag & PARTICLE_EXT_FORCE) {
    Tcl_AppendResult(interp, " ext_force ", (char *)NULL);
    part_print_ext_force(part, buffer, interp);  
  }

  /* print fix information. */
  if (part.l.ext_flag & COORDS_FIX_MASK) {
    Tcl_AppendResult(interp, " fix ", (char *)NULL);
    part_print_fix( part, buffer, interp);
  }
#endif

  /* print bonding structure */
  if(bl->n > 0) {
    part_print_bonding_structure(part, buffer, interp, bl);
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
  IntList *bl = &(part.bl);
    
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
      part_print_position(part, buffer, interp);   
    else if (ARG0_IS_S("force"))
      part_print_f(part, buffer, interp);   
    else if (ARG0_IS_S("folded_position"))
      part_print_folded_position(part, buffer, interp);   
    else if (ARG0_IS_S("type")) {
      sprintf(buffer, "%d", part.p.type);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else if (ARG0_IS_S("molecule_id")) {
      sprintf(buffer, "%d", part.p.mol_id);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#ifdef ELECTROSTATICS
    else if (ARG0_IS_S("q")) {
      Tcl_PrintDouble(interp, part.p.q, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#endif
    else if (ARG0_IS_S("v"))
      part_print_v(part, buffer, interp);   

#ifdef ROTATION
    else if (ARG0_IS_S("quat"))
      part_print_quat(part, buffer, interp);   
    else if (ARG0_IS_S("omega"))
      part_print_omega(part, buffer, interp);   
    else if (ARG0_IS_S("torque"))
      part_print_torque(part, buffer, interp);   
#endif         

#ifdef EXTERNAL_FORCES
    else if (ARG0_IS_S("ext_force"))
      part_print_ext_force(part,buffer,interp);
    else if (ARG0_IS_S("fix"))
      part_print_fix(part, buffer, interp);
#endif
    else if (ARG0_IS_S("bonds"))
      part_print_bonding_structure(part, buffer, interp, bl);
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
  if (! ARG0_IS_I(mid))
    return TCL_ERROR;
      
  if (mid < 0) {
    Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
    return TCL_ERROR;	  
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
  
  /* check type_num */
  if (argc <= 1 || ! ARG0_IS_I(type_num)) {
    if (delete == 0)
      return TCL_ERROR;
    bond = NULL;
   
    // since there is not even a type...
    n_partners = -1;
  } else {
    if(type_num < 0 || type_num >= n_bonded_ia) {
      Tcl_AppendResult(interp, "invalid bonded interaction type_num"
		       "(Set bonded interaction parameters first)", (char *) NULL);
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
      if(bond[j] > max_seen_particle || particle_node[bond[j]] == -1) {
	Tcl_AppendResult(interp, "partner atom %d (identity %d) not known"
			 ,j+1,bond[j],
			 "(Set all partner atoms first)", (char *) NULL);
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

int part_parse_cmd(Tcl_Interp *interp, int argc, char **argv,
		   int part_num)
{
  int change = 0, err = TCL_OK;

  while (argc > 0) {
    if (ARG0_IS_S("print"))
      return part_parse_print(interp, argc-1, argv+1, part_num);

    else if (ARG0_IS_S("delete"))
      err = part_cmd_delete(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("pos"))
      err = part_parse_pos(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("type"))
      err = part_parse_type(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("molecule_id"))
      err = part_parse_mol_id(interp, argc-1, argv+1, part_num, &change);

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

    else {
      Tcl_AppendResult(interp, "unknown particle parameter \"",
		       argv[0],"\"", (char *)NULL);

      return TCL_ERROR;
    }
  
    argc -= (change + 1);
    argv += (change + 1);

    if (err == TCL_ERROR)
      return TCL_ERROR;
  }

  return TCL_OK;
}

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;

  /* if no further arguments are given, print out all stored particles */
  if (argc == 1)
    return part_print_all(interp);

  /* only one argument is given */
  if (argc == 2) {

    if (ARG1_IS_S("deleteall")) {
      remove_all_particles();

      return TCL_OK;
    }
    
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
  if(delete != 0 || bond==NULL) { 
    delete = 1; bond=NULL; }
  else {
    if(bond[0] < 0 ) {
      fprintf(stderr, "ERROR: Failed upon changing bonds of particle %d due to invalid bonded interaction type_num %d (must be >0)!\n", part,bond[0]);
      return TCL_ERROR; }
    if(bond[0] >= n_bonded_ia) {
      fprintf(stderr, "ERROR: Failed upon changing bonds of particle %d due to currently unknown bonded interaction %d!\n", part,bond[0]);
      fprintf(stderr, "       Specify its parameters first (using 'inter <bond_type_number> { FENE | Angle | ... } <parameters> ...') "
	      "before adding bonds of that type!\n");
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
    fprintf(stderr, "%d: could not find cell of particle %d, exiting\n",
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

void local_place_particle(int part, double p[3])
{
  Cell *cell;
  double pp[3];
  int i[3], rl;
  Particle *pt = (part <= max_seen_particle) ? local_particles[part] : NULL;

  i[0] = 0;
  i[1] = 0;
  i[2] = 0;
  pp[0] = p[0];
  pp[1] = p[1];
  pp[2] = p[2];
  fold_position(pp, i);
  
  if (!pt) {
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

  PART_TRACE(fprintf(stderr, "%d: local_place_particle: got particle id=%d @ %f %f %f\n",
		     this_node, part, p[0], p[1], p[2]));

  memcpy(pt->r.p, pp, 3*sizeof(double));
  memcpy(pt->l.i, i, 3*sizeof(int));
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
		 bl->n - i);
	  realloc_intlist(bl, bl->n);
	}
	else
          i += 1 + partners;
      }
      if (i != bl->n) {
	fprintf(stderr, "%d: bond information corrupt for particle %d, exiting...\n",
		this_node, part[p].p.identity);
	errexit();
      }
    }
  }
}

void send_particles(ParticleList *particles, int node)
{
  int pc;
  IntList local_bi;

  PART_TRACE(fprintf(stderr, "%d: send_particles %d to %d\n", this_node, particles->n, node));

  MPI_Send(&particles->n, 1, MPI_INT, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);
  MPI_Send(particles->part, particles->n*sizeof(Particle),
	   MPI_BYTE, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);

  init_intlist(&local_bi);
  for (pc = 0; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    realloc_intlist(&local_bi, local_bi.n + p->bl.n);
    memcpy(&local_bi.e[local_bi.n], p->bl.e, p->bl.n*sizeof(int));
    local_bi.n += p->bl.n;
  }

  PART_TRACE(fprintf(stderr, "%d: send_particles sending %d bond ints\n", this_node, local_bi.n));
  if (local_bi.n > 0) {
    MPI_Send(local_bi.e, local_bi.n*sizeof(int),
	     MPI_BYTE, node, REQ_SNDRCV_PART, MPI_COMM_WORLD);
    realloc_intlist(&local_bi, 0);
  }

  /* remove particles from this nodes local list */
  for (pc = 0; pc < particles->n; pc++) {
    local_particles[particles->part[pc].p.identity] = NULL;
    free_particle(&particles->part[pc]);
  }

  realloc_particlelist(particles, particles->n = 0);
}

void recv_particles(ParticleList *particles, int node)
{
  int transfer, read, pc;
  MPI_Status status;
  IntList local_bi;

  PART_TRACE(fprintf(stderr, "%d: recv_particles from %d\n", this_node, node));

  MPI_Recv(&transfer, 1, MPI_INT, node,
	   REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);

  PART_TRACE(fprintf(stderr, "%d: recv_particles get %d\n", this_node, transfer));

  realloc_particlelist(particles, particles->n + transfer);
  MPI_Recv(&particles->part[particles->n], transfer*sizeof(Particle), MPI_BYTE, node,
	   REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);
  particles->n += transfer;

  local_bi.n = 0;
  for (pc = particles->n - transfer; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    local_bi.n += p->bl.n;

    PART_TRACE(fprintf(stderr, "%d: recv_particles got particle %d\n", this_node, p->p.identity));
#ifdef ADDITIONAL_CHECKS
    if (local_particles[p->p.identity] != NULL) {
      fprintf(stderr, "%d: transmitted particle %d is already here...\n", this_node, p->p.identity);
      errexit();
    }
#endif
  }

  update_local_particles(particles);

  PART_TRACE(fprintf(stderr, "%d: recv_particles expecting %d bond ints\n", this_node, local_bi.n));
  if (local_bi.n > 0) {
    alloc_intlist(&local_bi, local_bi.n);
    MPI_Recv(local_bi.e, local_bi.n*sizeof(int), MPI_BYTE, node,
	     REQ_SNDRCV_PART, MPI_COMM_WORLD, &status);
  }
  read = 0;
  for (pc = particles->n - transfer; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    if (p->bl.n > 0) {
      alloc_intlist(&p->bl, p->bl.n);
      memcpy(p->bl.e, &local_bi.e[read], p->bl.n*sizeof(int));
      read += p->bl.n;
    }
    else
      p->bl.e = NULL;
  }
  if (local_bi.n > 0)
    realloc_intlist(&local_bi, 0);
}

