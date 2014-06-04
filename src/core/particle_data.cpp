/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file particle_data.cpp
    This file contains everything related to particle storage. If you want to add a new
    property to the particles, it is probably a good idea to modify \ref Particle to give
    scripts access to that property. You always have to modify two positions: first the
    print section, where you should add your new data at the end, and second the read
    section where you have to find a nice and short name for your property to appear in
    the Tcl code. Then you just parse your part out of argc and argv.

    The corresponding header file is particle_data.hpp.
*/
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include "utils.hpp"
#include "particle_data.hpp"
#include "global.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "rotation.hpp"
#include "virtual_sites.hpp"

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
// List of particles for grandcanonical simulations
TypeOfIndex Type; 
IndexOfType Index; 
TypeList *type_array;
int number_of_type_lists;
int GC_init;
int Type_array_init = 0;

int max_seen_particle = -1;
int n_part = 0;
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

#ifdef SHANCHEN
  int ii;
  for(ii=0;ii<2*LB_COMPONENTS;ii++){ 
    part->p.solvation[ii]=0.0;
  }
  for(ii=0;ii<LB_COMPONENTS;ii++){ 
    part->r.composition[ii]=0.0;
  }
#endif

#ifdef ROTATIONAL_INERTIA
  part->p.rinertia[0] = 1.0;
  part->p.rinertia[1] = 1.0;
  part->p.rinertia[2] = 1.0;
#endif
#ifdef ROTATION_PER_PARTICLE
  part->p.rotation =1;
#endif


#ifdef ELECTROSTATICS
  part->p.q        = 0.0;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  part->p.mu_E[0]   = 0.0;
  part->p.mu_E[1]   = 0.0;
  part->p.mu_E[2]   = 0.0;
#endif

#ifdef CATALYTIC_REACTIONS
  part->p.catalyzer_count = 0;
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

  part->r.quatu[0]  = 0.0;
  part->r.quatu[1]  = 0.0;
  part->r.quatu[2]  = 1.0;
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

  #ifdef GHMC

    /* Last Saved ParticlePosition */
    part->l.r_ls.p[0]     = 0.0;
    part->l.r_ls.p[1]     = 0.0;
    part->l.r_ls.p[2]     = 0.0;

  #ifdef BOND_CONSTRAINT
    part->l.r_ls.p_old[0] = 0.0;
    part->l.r_ls.p_old[1] = 0.0;
    part->l.r_ls.p_old[2] = 0.0;
  #endif

  #ifdef ROTATION
    part->l.r_ls.quat[0]  = 1.0;
    part->l.r_ls.quat[1]  = 0.0;
    part->l.r_ls.quat[2]  = 0.0;
    part->l.r_ls.quat[3]  = 0.0;

    part->l.r_ls.quatu[0]  = 0.0;
    part->l.r_ls.quatu[1]  = 0.0;
    part->l.r_ls.quatu[2]  = 1.0;
  #endif

  #ifdef DIPOLES
    part->l.r_ls.dip[0]    = 0.0;
    part->l.r_ls.dip[1]    = 0.0;
    part->l.r_ls.dip[2]    = 0.0;
    //part->l.p_ls.dipm      = 0.0;
  #endif

    /* Last Saved ParticleMomentum */
    part->l.m_ls.v[0]     = 0.0;
    part->l.m_ls.v[1]     = 0.0;
    part->l.m_ls.v[2]     = 0.0;
  #ifdef ROTATION
    part->l.m_ls.omega[0] = 0.0;
    part->l.m_ls.omega[1] = 0.0;
    part->l.m_ls.omega[2] = 0.0;
  #endif

#endif

#ifdef EXTERNAL_FORCES
  part->l.ext_flag   = 0;
  part->l.ext_force[0] = 0.0;
  part->l.ext_force[1] = 0.0;
  part->l.ext_force[2] = 0.0;
  #ifdef ROTATION
    part->l.ext_torque[0] = 0.0;
    part->l.ext_torque[1] = 0.0;
    part->l.ext_torque[2] = 0.0;
  #endif
#endif

  init_intlist(&(part->bl));
#ifdef EXCLUSIONS
  init_intlist(&(part->el));
#endif

#ifdef VIRTUAL_SITES
  part->p.isVirtual      = 0;
#endif

#ifdef VIRTUAL_SITES_RELATIVE
  part->p.vs_relative_to_particle_id      = 0;
  part->p.vs_relative_distance =0;
#endif

#ifdef GHOST_FLAG
  part->l.ghost        = 0;
#endif

#ifdef LANGEVIN_PER_PARTICLE
  part->p.T = -1.0;
  part->p.gamma = -1.0;
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

int updatePartCfg(int bonds_flag)
{
  int j;

  if(partCfg)
    return 1;

  partCfg = (Particle*)malloc(n_part*sizeof(Particle));
  if (bonds_flag != WITH_BONDS)
    mpi_get_particles(partCfg, NULL);
  else
    mpi_get_particles(partCfg,&partCfg_bl);

  for(j=0; j<n_part; j++)
    unfold_position(partCfg[j].r.p,partCfg[j].l.i);
  partCfgSorted = 0;
#ifdef VIRTUAL_SITES

  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return 0;
  }
  if (!updatePartCfg(bonds_flag)) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not update positions of virtual sites in partcfg } ");
    return 0;
  }
#endif
return 1;
}

int sortPartCfg()
{
  int i;
  Particle *sorted;

  if (!partCfg)
    updatePartCfg(WITHOUT_BONDS);

  if (partCfgSorted)
    return 1;

  if (n_part != max_seen_particle + 1)
    return 0;

  sorted = (Particle*)malloc(n_part*sizeof(Particle));
  for(i = 0; i < n_part; i++)
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

int get_particle_data(int part, Particle *data)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;

  pnode = particle_node[part];
  if (pnode == -1)
    return ES_ERROR;
  mpi_recv_part(pnode, part, data);
  return ES_OK;
}

int place_particle(int part, double p[3])
{
  int i;
  int pnode, retcode = ES_PART_OK;

  if (part < 0)
    return ES_PART_ERROR;

  if (!particle_node)
    build_particle_node();

  pnode = (part <= max_seen_particle) ? particle_node[part] : -1;
  if (pnode == -1) {
    /* new particle, node by spatial position */
    pnode = cell_structure.position_to_node(p);

    /* master node specific stuff */
    realloc_particle_node(part);
    particle_node[part] = pnode;

    /* fill up possible gap */
    for (i = max_seen_particle + 1; i < part; i++)
      particle_node[i] = -1;

    retcode = ES_PART_CREATED;

    mpi_place_new_particle(pnode, part, p);

  } else {
    mpi_place_particle(pnode, part, p);
  }

  return retcode;
}

int set_particle_v(int part, double v[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_v(pnode, part, v);
  return ES_OK;
}

int set_particle_f(int part, double F[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_f(pnode, part, F);
  return ES_OK;
}

#ifdef SHANCHEN 
int set_particle_solvation(int part, double * solvation)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_solvation(pnode, part, solvation);
  return ES_OK;
}

#endif


#ifdef MASS
int set_particle_mass(int part, double mass)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_mass(pnode, part, mass);
  return ES_OK;
}
#endif

#ifdef ROTATIONAL_INERTIA
int set_particle_rotational_inertia(int part, double rinertia[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_rotational_inertia(pnode, part, rinertia);
  return ES_OK;
}
#endif


#ifdef ROTATION_PER_PARTICLE
int set_particle_rotation(int part, int rot)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_rotation(pnode, part, rot);
  return ES_OK;
}
#endif

#ifdef DIPOLES
int set_particle_dipm(int part, double dipm)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_dipm(pnode, part, dipm);
  return ES_OK;
}

int set_particle_dip(int part, double dip[3])
{
  int pnode;
  
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_dip(pnode, part, dip);

  return ES_OK;
}

#endif

#ifdef VIRTUAL_SITES
int set_particle_virtual(int part, int isVirtual)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_virtual(pnode, part, isVirtual); 
  return ES_OK;
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
int set_particle_vs_relative(int part, int vs_relative_to, double vs_distance)
{
  // Find out, on what node the particle is
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  
  // Send the stuff
  mpi_send_vs_relative(pnode, part, vs_relative_to, vs_distance);
  return ES_OK;
}
#endif

int set_particle_q(int part, double q)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_q(pnode, part, q);
  return ES_OK;
}

#ifdef LB_ELECTROHYDRODYNAMICS
int set_particle_mu_E(int part, double mu_E[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_mu_E(pnode, part, mu_E);
  return ES_OK;
}
#endif

int set_particle_type(int part, int type)
{

  int pnode;
  make_particle_type_exist(type);

  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;

// check if the particle exists already and the type is changed, then remove it from the list which contains it
  Particle *cur_par = (Particle *) malloc( sizeof(Particle) );
  if ( Type_array_init ) {
	  if ( cur_par != (Particle *) 0 ) {
		  if ( get_particle_data(part, cur_par) != ES_ERROR ) {
			  int prev_type = cur_par->p.type;
			  if ( prev_type != type ) {
			  // particle existed before so delete it from the list
			  remove_id_type_array(part, prev_type);
			  }
		  }
	  }
	  free(cur_par);
  }

  mpi_send_type(pnode, part, type);

  if ( Type_array_init ) { 
	  if ( add_particle_to_list(part, type) ==  ES_ERROR ){
		  //Tcl_AppendResult(interp, "gc particle add failed", (char *) NULL);
		  return ES_ERROR;
	  }
  }

  return ES_OK;
}

int set_particle_mol_id(int part, int mid)
{
  int pnode;

  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_mol_id(pnode, part, mid);
  return ES_OK;
}

#ifdef ROTATION
int set_particle_quat(int part, double quat[4])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_quat(pnode, part, quat);
  return ES_OK;
}

int set_particle_omega_lab(int part, double omega_lab[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;

  /* Internal functions require the body coordinates
     so we need to convert to these from the lab frame */

  double A[9];
  double omega[3];
  Particle particle;

  get_particle_data(part, &particle);
  define_rotation_matrix(&particle, A);

  omega[0] = A[0 + 3*0]*omega_lab[0] + A[0 + 3*1]*omega_lab[1] + A[0 + 3*2]*omega_lab[2];
  omega[1] = A[1 + 3*0]*omega_lab[0] + A[1 + 3*1]*omega_lab[1] + A[1 + 3*2]*omega_lab[2];
  omega[2] = A[2 + 3*0]*omega_lab[0] + A[2 + 3*1]*omega_lab[1] + A[2 + 3*2]*omega_lab[2];

  mpi_send_omega(pnode, part, omega);
  return ES_OK;
}

int set_particle_omega_body(int part, double omega[3])
{
  /* Nothing to be done but pass, since the coordinates
     are already in the proper frame */

  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_omega(pnode, part, omega);
  return ES_OK;
}

int set_particle_torque_lab(int part, double torque_lab[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;

  /* Internal functions require the body coordinates
     so we need to convert to these from the lab frame */

  double A[9];
  double torque[3];
  Particle particle;

  get_particle_data(part, &particle);
  define_rotation_matrix(&particle, A);

  torque[0] = A[0 + 3*0]*torque_lab[0] + A[0 + 3*1]*torque_lab[1] + A[0 + 3*2]*torque_lab[2];
  torque[1] = A[1 + 3*0]*torque_lab[0] + A[1 + 3*1]*torque_lab[1] + A[1 + 3*2]*torque_lab[2];
  torque[2] = A[2 + 3*0]*torque_lab[0] + A[2 + 3*1]*torque_lab[1] + A[2 + 3*2]*torque_lab[2];

  mpi_send_torque(pnode, part, torque);
  return ES_OK;
}

int set_particle_torque_body(int part, double torque[3])
{
  /* Nothing to be done but pass, since the coordinates
     are already in the proper frame */

  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_torque(pnode, part, torque);
  return ES_OK;
}

#endif

#ifdef LANGEVIN_PER_PARTICLE
int set_particle_temperature(int part, double T)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
    
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
    
  mpi_set_particle_temperature(pnode, part, T);
  return ES_OK;
}

int set_particle_gamma(int part, double gamma)
{
  int pnode;
  
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
    
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
    
  mpi_set_particle_gamma(pnode, part, gamma);
  return ES_OK;
}
#endif

#ifdef EXTERNAL_FORCES
  #ifdef ROTATION
    int set_particle_ext_torque(int part, int flag, double torque[3])
    {
      int pnode;
      if (!particle_node)
        build_particle_node();

      if (part < 0 || part > max_seen_particle)
        return ES_ERROR;
      pnode = particle_node[part];

      if (pnode == -1)
        return ES_ERROR;

      mpi_send_ext_torque(pnode, part, flag, PARTICLE_EXT_TORQUE, torque);
        return ES_OK;
    }
  #endif

int set_particle_ext_force(int part, int flag, double force[3])
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;

  mpi_send_ext_force(pnode, part, flag, PARTICLE_EXT_FORCE, force);
    return ES_OK;
}

int set_particle_fix(int part,  int flag)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_ext_force(pnode, part, flag, COORDS_FIX_MASK, NULL);
  return ES_OK;
}

#endif

int change_particle_bond(int part, int *bond, int _delete)
{
  int pnode;
  if (!particle_node)
    build_particle_node();

  if (part < 0 || part > max_seen_particle)
    return ES_ERROR;
  pnode = particle_node[part];

  if (pnode == -1)
    return ES_ERROR;
  if(_delete != 0 || bond == NULL)
    _delete = 1;

  if (bond != NULL) {
    if (bond[0] < 0 || bond[0] >= n_bonded_ia) {
      char *errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{048 invalid/unknown bonded interaction type %d}", bond[0]);
      return ES_ERROR;
    }
  }
  return mpi_send_bond(pnode, part, bond, _delete);
}

void remove_all_particles()
{
  mpi_remove_particle(-1, -1);
  realloc_particle_node(0);
}

int remove_particle(int part)
{
  int pnode;

  Particle *cur_par = (Particle *) malloc (sizeof(Particle));
  if (get_particle_data(part, cur_par) == ES_ERROR )
	  return ES_ERROR;
  int type = cur_par->p.type;
  free(cur_par);
  if (remove_id_type_array(part, type) == ES_ERROR )
	  return ES_ERROR;

  if (!particle_node)
    build_particle_node();

  if (part > max_seen_particle)
    return ES_ERROR;

  pnode = particle_node[part];
  if (pnode == -1)
    return ES_ERROR;

  particle_node[part] = -1;

  mpi_remove_particle(pnode, part);

  if (part == max_seen_particle) {
    while (max_seen_particle >= 0 && particle_node[max_seen_particle] == -1)
      max_seen_particle--;
    mpi_bcast_parameter(FIELD_MAXPART);
  }
  return ES_OK;
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

void local_place_particle(int part, double p[3], int _new)
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
  
  if (_new) {
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
  n_part = 0;
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

  n_part++;

  if (part > max_seen_particle) {
    realloc_local_particles(part);
    /* fill up possible gap. Part itself is ESSENTIAL!!!  */
    for (i = max_seen_particle + 1; i <= part; i++)
      local_particles[i] = NULL;
    max_seen_particle = part;
  }
}

int local_change_bond(int part, int *bond, int _delete)
{
  IntList *bl;
  Particle *p;
  int bond_size;
  int i;

  p = local_particles[part];
  if (_delete)
    return try_delete_bond(p, bond);

  bond_size = bonded_ia_params[bond[0]].num + 1;
  bl = &(p->bl);
  realloc_intlist(bl, bl->n + bond_size);
  for(i = 0; i < bond_size; i++)
    bl->e[bl->n++] = bond[i];
  return ES_OK;
}

int try_delete_bond(Particle *part, int *bond)
{
  IntList *bl = &part->bl;
  int i, j, type, partners;

  if (!bond) {
    realloc_intlist(bl, bl->n = 0);
    return ES_OK;
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
	return ES_OK;
      }
      i += 1 + partners;
    }
  }
  return ES_ERROR;
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
void local_change_exclusion(int part1, int part2, int _delete)
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
    if (_delete)
      try_delete_exclusion(part, part2);
    else
      try_add_exclusion(part, part2);
  }

  /* part2, if here */
  part = local_particles[part2];
  if (part) {
    if (_delete)
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

  MPI_Send(&particles->n, 1, MPI_INT, node, REQ_SNDRCV_PART, comm_cart);
  MPI_Send(particles->part, particles->n*sizeof(Particle),
	   MPI_BYTE, node, REQ_SNDRCV_PART, comm_cart);

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
	     MPI_BYTE, node, REQ_SNDRCV_PART, comm_cart);
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
  IntList local_dyn;

  PART_TRACE(fprintf(stderr, "%d: recv_particles from %d\n", this_node, node));

  MPI_Recv(&transfer, 1, MPI_INT, node,
	   REQ_SNDRCV_PART, comm_cart, MPI_STATUS_IGNORE);

  PART_TRACE(fprintf(stderr, "%d: recv_particles get %d\n", this_node, transfer));

  realloc_particlelist(particles, particles->n + transfer);
  MPI_Recv(&particles->part[particles->n], transfer*sizeof(Particle), MPI_BYTE, node,
	   REQ_SNDRCV_PART, comm_cart, MPI_STATUS_IGNORE);
  particles->n += transfer;

  init_intlist(&local_dyn);
  for (pc = particles->n - transfer; pc < particles->n; pc++) {
    Particle *p = &particles->part[pc];
    local_dyn.n += p->bl.n;
#ifdef EXCLUSIONS
    local_dyn.n += p->el.n;
#endif

    PART_TRACE(fprintf(stderr, "%d: recv_particles got particle %d\n", this_node, p->p.identity));
    if (local_particles[p->p.identity] != NULL) {
      fprintf(stderr, "%d: transmitted particle %d is already here...\n", this_node, p->p.identity);
      errexit();
    }
  }

  update_local_particles(particles);

  PART_TRACE(fprintf(stderr, "%d: recv_particles expecting %d bond ints\n", this_node, local_dyn.n));
  if (local_dyn.n > 0) {
    alloc_intlist(&local_dyn, local_dyn.n);
    MPI_Recv(local_dyn.e, local_dyn.n*sizeof(int), MPI_BYTE, node,
	     REQ_SNDRCV_PART, comm_cart, MPI_STATUS_IGNORE);
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

void add_partner(IntList *il, int i, int j, int distance)
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

#ifdef EXCLUSIONS

int change_exclusion(int part1, int part2, int _delete)
{
  if (!particle_node)
    build_particle_node();

  if (part1 < 0 || part1 > max_seen_particle ||
      part2 < 0 || part2 > max_seen_particle ||
      part1 == part2 ||
      particle_node[part1] == -1 ||
      particle_node[part2] == -1)
    return ES_ERROR;

  mpi_send_exclusion(part1, part2, _delete);
  return ES_OK;
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
  partners    = (IntList*)malloc((max_seen_particle + 1)*sizeof(IntList));
  for (p = 0; p <= max_seen_particle; p++)
    init_intlist(&partners[p]);

  /* determine initial connectivity */
  for (p = 0; p < n_part; p++) {
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

#endif


int init_gc(void){
	if ( type_array == (TypeList *) NULL) {
	//stores the number of currently available type_list's 
	number_of_type_lists=10;

	Type.max_entry = 0;
	Index.max_entry = 0;

	type_array = (TypeList *) malloc(sizeof(TypeList) * number_of_type_lists);
	if ( type_array == (TypeList *) 0 )
		return ES_ERROR;

	GC_init = 1;
	Type_array_init = 0;
	}
	return ES_OK;
}

int init_type_array(int type){
	if (init_gc() == ES_ERROR)
		return ES_ERROR;

	updatePartCfg(WITHOUT_BONDS);

	if (!partCfg)
		return ES_ERROR;

	int type_index=-1;
	type_index= (Type.max_entry++);
	if ( type_index == number_of_type_lists ) {
		reallocate_global_type_list(number_of_type_lists*2);
	}

	Type.index = (int *) realloc ( (void *) Type.index, sizeof(int)*Type.max_entry);

	//reallocate the array that holds the particle type and points to the type index used for the type_list
	
	if ( type >= Index.max_entry ) { 
		Index.type= (int * ) realloc( (void *) Index.type, (type+1)*sizeof(int));
		Index.max_entry = type + 1;
	}
	for (int i=0; i<Type.max_entry; i++) 
		Index.type[i]=-1;

	if (Type.index == (int *) 0 || Index.type == (int *) 0)
		return ES_ERROR;

	//allocates a list for ids for as many entries as there are particles right now
	if ( !(partCfg) ||  type < 0 ) { 
		return ES_ERROR;
	}
	Type.index[type_index]=type;
	//fill in array type_index_of_type
	for (int i=0; i<Type.max_entry; i++ ) {
		Index.type[Type.index[i]] = i;
	}

	int t_c = 0; //index
	type_array[Index.type[type]].id_list = (int *) malloc (sizeof (int) * n_part);
	for (int i=0; i<n_part; i++) {
		if ( partCfg[i].p.type==type ) 
			type_array[Index.type[type]].id_list[t_c++]=partCfg[i].p.identity;
	}
	int max_size=n_part;
	if ( t_c != 0 ) { 
		while ( t_c < (double) max_size/4.0) { 
			max_size= floor( (double ) max_size/2.0);
		}
		// now the array is shrinked to at least 4 times the highest entry
		type_array[Index.type[type]].id_list= (int *) realloc( (void *) type_array[Index.type[type]].id_list, sizeof(int)*2*max_size);
		type_array[Index.type[type]].max_entry = t_c;
		type_array[Index.type[type]].cur_size = max_size*2;
	} else {
		//no particles of the given type were found, so leave array size fixed at a reasonable start entry 64 ints in this case
		type_array[Index.type[type]].id_list= (int *) realloc( (void *) type_array[Index.type[type]].id_list, sizeof(int)*64);
		type_array[Index.type[type]].max_entry = t_c;
		type_array[Index.type[type]].cur_size = 64;
	}
	//fill remaining entries with -1
	for (int i=type_array[Index.type[type]].max_entry; i<type_array[Index.type[type]].cur_size; i++) {
		type_array[Index.type[type]].id_list[i] = -1;
	}
	Type_array_init = 1;
	return ES_OK;
}

int reallocate_type_array(int type){
	type_array[Index.type[type]].id_list = (int *) realloc ( (void *) type_array[Index.type[type]].id_list, sizeof (int) * type_array[Index.type[type]].cur_size * 2);	
	if (type_array[Index.type[type]].id_list == (int *) 0) {
		return ES_ERROR;
	}
	type_array[Index.type[type]].cur_size = type_array[Index.type[type]].cur_size*2;
	return ES_OK;
}

int remove_id_type_array(int part_id, int type){

	int l_err = 1;
	for ( int j = 0; j<Type.max_entry; j++){
		if ( Type.index[j] == type ) {
			l_err = 0;
			break;
		}
	}
	if ( l_err ) { 
		//there is no list which contains this type
		return ES_OK;
	}
	int in_type = Index.type[type];
	int temp_id = -1;
	int max = type_array[in_type].max_entry;
	for (int i=0; i<max; i++){
		if ( type_array[in_type].id_list[i] == part_id ) {
			temp_id = i;
			break;
		}
	}
	if ( temp_id == -1 ){
		//particle is not in the list
		return ES_OK;
	}
	if ( temp_id == max - 1 ) { 
		type_array[in_type].id_list[temp_id] = -1;
	} else {
		int temp=type_array[in_type].id_list[max-1]; 
		type_array[in_type].id_list[max-1] = -1;
		type_array[in_type].id_list[temp_id]=temp;
	}
	type_array[in_type].max_entry--;
	return ES_OK;
}

int update_particle_array(int type) {
	updatePartCfg(WITHOUT_BONDS);
	if (!partCfg) 
		return ES_ERROR;

	int t_c = 0;
	for (int i=0; i<n_part; i++) {
		if (partCfg[i].p.type == type ) {
			type_array[Index.type[type]].id_list[t_c++] = partCfg[i].p.identity;
		}	
		if ( t_c > (double) type_array[Index.type[type]].cur_size/2.0 ) {
			if ( reallocate_type_array(type) == ES_ERROR ) 
				return ES_ERROR;
		}
	}
	type_array[Index.type[type]].max_entry = t_c;
	for ( int i = t_c; i<type_array[Index.type[type]].cur_size; i++)
		type_array[Index.type[type]].id_list[i] = -1;

	return ES_OK;
}

int reallocate_global_type_list(int size){
	if (size <= 0 ) 
		return ES_ERROR;
	type_array = (TypeList *) realloc( (void *) type_array, sizeof(TypeList)*size);
	number_of_type_lists=size;
	if ( type_array == (TypeList *) 0)
		return ES_ERROR;

	return ES_OK;
}

int find_particle_type(int type, int *id){ 
	int l_err=1;
	// type i not indexed, so no list for this particle exists
	for (int i = 0; i<Type.max_entry; i++){
		if ( Type.index[i] == type ){
			l_err=0;
			break;
		}
	}
	if ( l_err ) {
		return ES_ERROR;
	}
	if ( type_array[Index.type[type]].max_entry == 0 ){ 
		return ES_ERROR;
	}
	int rand_index = i_random ( type_array[Index.type[type]].max_entry );
	*id = type_array[Index.type[type]].id_list[rand_index];

	return ES_OK;
}

int find_particle_type_id(int type, int *id, int *in_id ){
	int l_err=1;
	// type i not indexed, so no list for this particle exists
	for (int i = 0; i<Type.max_entry; i++){
		if ( Type.index[i] == type ){
			l_err=0;
			break;
		}
	}
	if ( l_err ) {
		return ES_ERROR;
	}
	if ( type_array[Index.type[type]].max_entry == 0 )
		return ES_ERROR;

	int rand_index = i_random ( type_array[Index.type[type]].max_entry );
	if (id == (int *) 0 && in_id == (int *) 0) 
		return ES_ERROR;
	else {
		*in_id = rand_index;
		*id = type_array[Index.type[type]].id_list[*in_id];
		return ES_OK;
	}
}

int delete_particle_of_type(int type) { 
	int *p_id, *index_id;
	p_id=(int *) malloc (sizeof(int));
	index_id=(int *) malloc (sizeof(int));
	if (find_particle_type_id(type, p_id, index_id) == ES_ERROR )
		return ES_ERROR;

	int in_type = Index.type[type];
	// maximal possible index id
	int max = type_array[in_type].max_entry - 1;
	if ( max < 0 ) 
		return ES_ERROR;

	if ( remove_particle(*p_id) == ES_ERROR ) {
		// takes also care of removing the index from the array
		return ES_ERROR;
	}
	return ES_OK;
}

int add_particle_to_list(int part_id, int type){
	int l_err=1;
	int already_in = 0;
//	int already_in_other_list = 0;
	// type i not indexed, so no list for this particle exists
	for (int i = 0; i<Type.max_entry; i++){
		if ( Type.index[i] == type ){
			l_err=0;
			break;
		}
	}
	if ( l_err ) {
		return NOT_INDEXED;
	}

	int in_type = Index.type[type];
	int max = type_array[in_type].max_entry;
	for ( int i=0; i<max; i++ ) {
		if ( type_array[in_type].id_list[i] == part_id ) {
			already_in = 1;
			break;
		}
	}
	if ( already_in ) { 
		return ES_OK;
	}

	if ( max >= (double ) type_array[in_type].cur_size/2.0 )
		if (reallocate_type_array(type)== ES_ERROR)
			return ES_ERROR;

	//add particle id to list:
	type_array[in_type].id_list[max]=part_id;
	type_array[in_type].max_entry++;
	return ES_OK;
}

int gc_status(int type){
	int l_err=1;
	// type i not indexed, so no list for this particle exists
	for (int i = 0; i<Type.max_entry; i++){
		if ( Type.index[i] == type ){
			l_err=0;
			break;
		}
	}
	if ( l_err ) {
		return ES_ERROR;
	}
	int in_type = Index.type[type];
	for (int i = 0; i<type_array[in_type].max_entry; i++) {
		printf("%d\n", type_array[in_type].id_list[i]);
	}
	return ES_OK;
}

int free_particle_lists(void){
	if (type_array == (TypeList *) 0 || Type.index == (int *) 0){
		return ES_OK;
	}
	for (int i=0; i<Type.max_entry; i++){
		free(type_array[i].id_list);
	}
	free(type_array);
	free(Type.index);
	free(Index.type);
	return ES_OK;
}


int number_of_particles_with_type(int type, int *number){
	int indexed=0;
	if ( type_array == (TypeList *) 0 ) 
		init_type_array(type);

	for ( int i = 0; i<Type.max_entry; i++) {
		if ( type == Type.index[i] ){
			indexed=1;
			break;
		}
	}
	if ( indexed ) {
		*number = type_array[Index.type[type]].max_entry;
		return ES_OK;
	}
	return NOT_INDEXED;
}






// The following functions are used by the python interface to obtain 
// properties of a particle, which are only compiled in in some configurations
// This is needed, because cython does not support conditional compilation 
// within a ctypedef definition


#ifdef ROTATION
void pointer_to_omega_body(Particle* p, double*&  res)
{
  res=p->m.omega;
}

void pointer_to_torque_lab(Particle* p, double*& res)
{
  res=p->f.torque;
}

void pointer_to_quat(Particle* p, double*& res)
{
  res=p->r.quat;
}

void pointer_to_quatu(Particle* p, double*& res)
{
  res=p->r.quatu;
}
#endif

#ifdef ELECTROSTATICS
void pointer_to_q(Particle* p, double*& res)
{
 res =&(p->p.q);
}
#endif

#ifdef VIRTUAL_SITES
void pointer_to_virtual(Particle* p, int*&  res)
{
  res=&(p->p.isVirtual);
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_relative(Particle* p, int*& res1,double*& res2)
{
  res1=&(p->p.vs_relative_to_particle_id);
  res2=&(p->p.vs_relative_distance);
}
#endif


#ifdef MASS
void pointer_to_mass(Particle* p, double*&  res)
{
  res=&(p->p.mass);
}
#endif


#ifdef DIPOLES
void pointer_to_dip(Particle* p, double*& res)
{
res=p->r.dip;
}

void pointer_to_dipm(Particle* p, double*& res)
{
res=&(p->p.dipm);
}
#endif

