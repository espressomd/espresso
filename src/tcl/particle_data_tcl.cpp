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
#include "parser.hpp"
#include "rotation.hpp"
#include "virtual_sites.hpp"


/* Add a link of size size consisting of (link[l], p) */
void add_link(IntList *il, IntList *link, int l, int p, int size)
{
  int i;
  realloc_intlist(il, il->n + size);
  for ( i = 0; i < size-1; i++ ) {
    il->e[il->n++] = link->e[l+i];
  } 
  il->e[il->n++] = p;
}

#ifdef ROTATIONAL_INERTIA
void tclcommand_part_print_rotational_inertia(Particle *part, char *buffer, Tcl_Interp *interp)
  {double rinertia[3];

  rinertia[0]=part->p.rinertia[0];
  rinertia[1]=part->p.rinertia[1];
  rinertia[2]=part->p.rinertia[2];

  Tcl_PrintDouble(interp, rinertia[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, rinertia[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, rinertia[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

#ifdef ROTATION

/* tclcommand_part_print_*_body_frame: function to have the possiblility of
   printing the angular velocities and torques in the body frame" */ 

void tclcommand_part_print_omega_body_frame(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->m.omega[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.omega[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.omega[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void tclcommand_part_print_torque_body_frame(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->f.torque[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.torque[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.torque[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

/* tclcommand_part_print_*_lab_frame: function to have the possiblility of
   printing the angular velocities and torques in the laboratory frame" */ 

void tclcommand_part_print_omega_lab_frame(Particle *part, char *buffer, Tcl_Interp *interp)
{
  double omega[3];
//in Espresso angular velocities are in body-fixed frames. We should convert they to the space-fixed coordinates.
  convert_omega_body_to_space(part, omega);

  Tcl_PrintDouble(interp, omega[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, omega[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, omega[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void tclcommand_part_print_torque_lab_frame(Particle *part, char *buffer, Tcl_Interp *interp)
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



void tclcommand_part_print_quat(Particle *part, char *buffer, Tcl_Interp *interp)
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

/* TODO: This function is not used anywhere. To be removed?  */
void tclcommand_part_print_quatu(Particle *part, char *buffer, Tcl_Interp *interp)
{
  Tcl_PrintDouble(interp, part->r.quatu[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.quatu[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->r.quatu[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

#ifdef DIPOLES
void tclcommand_part_print_dip(Particle *part, char *buffer, Tcl_Interp *interp)
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
void tclcommand_part_print_virtual(Particle *part, char *buffer, Tcl_Interp *interp)
{
  sprintf(buffer,"%i", part->p.isVirtual);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

#ifdef ROTATION_PER_PARTICLE
void tclcommand_part_print_rotation(Particle *part, char *buffer, Tcl_Interp *interp)
{
  sprintf(buffer,"%i", part->p.rotation);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif

void tclcommand_part_print_v(Particle *part, char *buffer, Tcl_Interp *interp)
{
  /* unscale velocities ! */
  Tcl_PrintDouble(interp, part->m.v[0]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.v[1]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->m.v[2]/time_step, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

#ifdef LB_ELECTROHYDRODYNAMICS
void tclcommand_part_print_mu_E(Particle *part, char *buffer, Tcl_Interp *interp)
{
  /* unscale velocities ! */
  Tcl_PrintDouble(interp, part->p.mu_E[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->p.mu_E[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->p.mu_E[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}
#endif


#ifdef SHANCHEN
void tclcommand_part_print_solvation(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int ii;
  for(ii=0;ii<LB_COMPONENTS;++ii){
     Tcl_PrintDouble(interp, part->p.solvation[ii], buffer);
     Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
}


void tclcommand_part_print_composition(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int ii;
  for(ii=0;ii<LB_COMPONENTS;++ii){
     Tcl_PrintDouble(interp, part->r.composition[ii], buffer);
     Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
}

#endif 

void tclcommand_part_print_f(Particle *part, char *buffer, Tcl_Interp *interp)
{
  /* unscale forces ! */
  Tcl_PrintDouble(interp, part->f.f[0]*PMASS(*part)/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.f[1]*PMASS(*part)/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part->f.f[2]*PMASS(*part)/(0.5*time_step*time_step), buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

void tclcommand_part_print_position(Particle *part, char *buffer, Tcl_Interp *interp)
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

void tclcommand_part_print_folded_position(Particle *part, char *buffer, Tcl_Interp *interp)
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

void tclcommand_part_print_bonding_structure(Particle *part, char *buffer, Tcl_Interp *interp)
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


 
/** Return all bond partners of a particle including bonds that are not stored at the particle itself up to a certain distance in numbers of bonds. Return a tcl list to the interpreter c*/
void tclcommand_part_print_bond_partners(Particle *part, char *buffer, Tcl_Interp *interp, int distance)
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
  partners    = (IntList*)malloc((max_seen_particle + 1)*sizeof(IntList));
  for (p = 0; p <= max_seen_particle; p++) init_intlist(&partners[p]);
  updatePartCfg(WITH_BONDS);

  /* determine initial connectivity */
  for (p = 0; p < n_part; p++) {
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
  links    = (IntList*)malloc((distance+1)*sizeof(IntList));
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
void tclcommand_part_print_exclusions(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int i;
  for (i = 0; i < part->el.n; i++) {
    sprintf(buffer, "%d ", part->el.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
}
#endif

#ifdef EXTERNAL_FORCES
void tclcommand_part_print_fix(Particle *part, char *buffer, Tcl_Interp *interp)
{
  int i;
  for (i = 0; i < 3; i++) {
    if (part->l.ext_flag & COORD_FIXED(i))
      Tcl_AppendResult(interp, "1 ", (char *)NULL);
    else
	    Tcl_AppendResult(interp, "0 ", (char *)NULL);
  }
}

void tclcommand_part_print_ext_force(Particle *part, char *buffer, Tcl_Interp *interp)
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

#ifdef ROTATION
void tclcommand_part_print_ext_torque(Particle *part, char *buffer, Tcl_Interp *interp)
{
  if(part->l.ext_flag & PARTICLE_EXT_TORQUE) {
    Tcl_PrintDouble(interp, part->l.ext_torque[0], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part->l.ext_torque[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part->l.ext_torque[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  else {
    Tcl_AppendResult(interp, "0.0 0.0 0.0 ", (char *)NULL);
  }
}
#endif
#endif

/** append particle data in ASCII form to the Tcl result.
    @param part_num the particle which data is appended
    @param interp   the Tcl interpreter to which result to add to */
int tclprint_to_result_Particle(Tcl_Interp *interp, int part_num)
{
  /* print to the buffer AT MOST a double or an integer, plus
     at most 8 letters (like "{ ")  */
  char buffer[8 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle part;

  if (get_particle_data(part_num, &part) == TCL_ERROR)
    return (TCL_ERROR);

  sprintf(buffer, "%d", part.p.identity);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_AppendResult(interp, " pos ", (char *)NULL);
  tclcommand_part_print_position(&part, buffer, interp);
  sprintf(buffer, "%d", part.p.type);
  Tcl_AppendResult(interp, " type ", buffer, (char *)NULL);

  if (part.p.mol_id > -1) {
    sprintf(buffer, "%d", part.p.mol_id);
    Tcl_AppendResult(interp, " molecule ", buffer, (char *)NULL);
  }
#ifdef MASS
  Tcl_PrintDouble(interp, part.p.mass, buffer);
  Tcl_AppendResult(interp, " mass ", buffer, (char *)NULL);
#endif
#ifdef ELECTROSTATICS
  Tcl_PrintDouble(interp, part.p.q, buffer);
  Tcl_AppendResult(interp, " q ", buffer, (char *)NULL);
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  Tcl_AppendResult(interp, " mu_E ", (char *)NULL);
  tclcommand_part_print_mu_E(&part, buffer, interp);
#endif
  Tcl_AppendResult(interp, " v ", (char *)NULL);
  tclcommand_part_print_v(&part, buffer, interp);
  Tcl_AppendResult(interp, " f ", (char *)NULL);
  tclcommand_part_print_f(&part, buffer, interp);

#ifdef ROTATION
  /* print information about rotation */
  Tcl_AppendResult(interp, " quat ", (char *)NULL);
  tclcommand_part_print_quat(&part, buffer, interp);

  Tcl_AppendResult(interp, " omega_lab ", (char *)NULL);
  tclcommand_part_print_omega_lab_frame(&part, buffer, interp);

  Tcl_AppendResult(interp, " omega_body ", (char *)NULL);
  tclcommand_part_print_omega_body_frame(&part, buffer, interp);

  Tcl_AppendResult(interp, " torque_lab ", (char *)NULL);
  tclcommand_part_print_torque_lab_frame(&part, buffer, interp);

  Tcl_AppendResult(interp, " torque_body ", (char *)NULL);
  tclcommand_part_print_torque_body_frame(&part, buffer, interp);
#endif

#ifdef ROTATION_PER_PARTICLE
  Tcl_AppendResult(interp, " rotation ", (char *)NULL);
  tclcommand_part_print_rotation(&part, buffer, interp);
#endif

#ifdef ROTATIONAL_INERTIA
  /* print information about rotational inertia */
  Tcl_AppendResult(interp, " rinertia ", (char *)NULL);
  tclcommand_part_print_rotational_inertia(&part, buffer, interp);
#endif

#ifdef DIPOLES
#ifndef ROTATION
  /* print full information about dipoles, no quaternions to keep
     the rest */
  Tcl_AppendResult(interp, " dip ", (char *)NULL);
  tclcommand_part_print_dip(&part, buffer, interp);
#else
  /* quaternions are set, just put the scalar dipole moment, otherwise
     reading back fails. */
  Tcl_PrintDouble(interp, part.p.dipm, buffer);
  Tcl_AppendResult(interp, " dipm ", buffer, (char *)NULL);
#endif

#endif

#ifdef VIRTUAL_SITES
  /* print information about isVirtual */
  Tcl_AppendResult(interp, " virtual ", (char *)NULL);
  tclcommand_part_print_virtual(&part, buffer, interp);
#endif

#ifdef VIRTUAL_SITES_RELATIVE
  // print the particle attributes used by the "relative" implementation of virtual sites
  Tcl_AppendResult(interp, " vs_relative ", (char *)NULL);
  sprintf(buffer, "%d %f", part.p.vs_relative_to_particle_id, part.p.vs_relative_distance);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
#endif

#ifdef EXCLUSIONS
  if (part.el.n > 0) {
    Tcl_AppendResult(interp, " exclude ", (char *)NULL);
    tclcommand_part_print_exclusions(&part, buffer, interp);
  }
#endif

#ifdef SHANCHEN
  Tcl_AppendResult(interp, " solvation ", (char *)NULL);
  tclcommand_part_print_solvation(&part, buffer, interp);

  Tcl_AppendResult(interp, " composition ", (char *)NULL);
  tclcommand_part_print_composition(&part, buffer, interp);
#endif

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
  if (part.l.ext_flag & PARTICLE_EXT_TORQUE) {
    Tcl_AppendResult(interp, " ext_torque ", (char *)NULL);
    tclcommand_part_print_ext_torque(&part, buffer, interp);
  }
#endif

  /* print external force information. */
  if (part.l.ext_flag & PARTICLE_EXT_FORCE) {
    Tcl_AppendResult(interp, " ext_force ", (char *)NULL);
    tclcommand_part_print_ext_force(&part, buffer, interp);
  }

  /* print fix information. */
  if (part.l.ext_flag & COORDS_FIX_MASK) {
    Tcl_AppendResult(interp, " fix ", (char *)NULL);
    tclcommand_part_print_fix(&part, buffer, interp);
  }
#endif

  /* print bonding structure */
  if (part.bl.n > 0) {
    Tcl_AppendResult(interp, " bond ", (char *)NULL);
    tclcommand_part_print_bonding_structure(&part, buffer, interp);
  }

  free_particle(&part);
  return (TCL_OK);
}

int tclcommand_part_print_all(Tcl_Interp *interp)
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

      tclprint_to_result_Particle(interp, i);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
  }

  return TCL_OK;
}

int tclcommand_part_parse_print(Tcl_Interp *interp, int argc, char **argv,
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
      tclcommand_part_print_position(&part, buffer, interp);
    else if (ARG0_IS_S("force"))
      tclcommand_part_print_f(&part, buffer, interp);
    else if (ARG0_IS_S("folded_position"))
      tclcommand_part_print_folded_position(&part, buffer, interp);
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
#ifdef SHANCHEN
    else if (ARG0_IS_S("solvation")) {
      tclcommand_part_print_solvation(&part, buffer, interp);
    }
    else if (ARG0_IS_S("composition")) {
      tclcommand_part_print_composition(&part, buffer, interp);
    }
#endif
#ifdef ELECTROSTATICS
    else if (ARG0_IS_S("q")) {
      Tcl_PrintDouble(interp, part.p.q, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
    else if (ARG0_IS_S("mu_E")) {
      tclcommand_part_print_mu_E(&part, buffer, interp);
    }
#endif
    else if (ARG0_IS_S("v"))
      tclcommand_part_print_v(&part, buffer, interp);

#ifdef ROTATION
   else if (ARG_IS_S_EXACT(0,"quatu"))
      tclcommand_part_print_quatu(&part, buffer, interp);
    else if (ARG0_IS_S("quat"))
      tclcommand_part_print_quat(&part, buffer, interp);
    else if (ARG0_IS_S("omega") || ARG0_IS_S_EXACT("omega_lab"))
      tclcommand_part_print_omega_lab_frame(&part, buffer, interp);
    else if (ARG0_IS_S_EXACT("omega_body"))
      tclcommand_part_print_omega_body_frame(&part, buffer, interp);
    else if ( ARG0_IS_S("torque") || ARG0_IS_S_EXACT("torque_lab") )
      tclcommand_part_print_torque_lab_frame(&part, buffer, interp);
    else if ( ARG0_IS_S_EXACT("torque_body") || ARG0_IS_S("tbf") )
      tclcommand_part_print_torque_body_frame(&part, buffer, interp);
#endif
#ifdef ROTATION_PER_PARTICLE
    else if (ARG0_IS_S("rotation"))
      tclcommand_part_print_rotation(&part, buffer, interp);
#endif

#ifdef ROTATIONAL_INERTIA
    else if (ARG0_IS_S("rinertia"))
      tclcommand_part_print_rotational_inertia(&part, buffer, interp);
#endif

#ifdef DIPOLES
    else if (ARG0_IS_S("dip"))
      tclcommand_part_print_dip(&part, buffer, interp);
    else if (ARG0_IS_S("dipm")) {
      Tcl_PrintDouble(interp, part.p.dipm, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
#endif

#ifdef VIRTUAL_SITES
    else if (ARG0_IS_S("virtual"))
      tclcommand_part_print_virtual(&part, buffer, interp);
#endif
#ifdef VIRTUAL_SITES_RELATIVE    
     else if (ARG0_IS_S("vs_relative")) {
       sprintf(buffer, "%d %f", part.p.vs_relative_to_particle_id, part.p.vs_relative_distance);
       Tcl_AppendResult(interp, buffer, (char *)NULL);
     }
#endif

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
    else if (ARG0_IS_S("ext_torque"))
      tclcommand_part_print_ext_torque(&part,buffer,interp);
#endif

    else if (ARG0_IS_S("ext_force"))
      tclcommand_part_print_ext_force(&part,buffer,interp);
    else if (ARG0_IS_S("fix"))
      tclcommand_part_print_fix(&part, buffer, interp);
#endif

#ifdef EXCLUSIONS
    else if (ARG0_IS_S("exclusions"))
      tclcommand_part_print_exclusions(&part, buffer, interp);
#endif

    else if (ARG0_IS_S("bonds"))
      tclcommand_part_print_bonding_structure(&part, buffer, interp);
    else if (ARG0_IS_S("connections")) {
      int distance = 1;
      if (argc ==2) {
	if (!ARG1_IS_I(distance)) {
	  free_particle(&part);
	  return TCL_ERROR;
	}
	argc--; argv++;
      }
      tclcommand_part_print_bond_partners(&part, buffer, interp, distance);
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


int tclcommand_part_parse_delete(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_pos(Tcl_Interp *interp, int argc, char **argv,
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

  if (place_particle(part_num, pos) == ES_PART_ERROR) {
    Tcl_AppendResult(interp, "particle could not be set", (char *) NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}


#ifdef SHANCHEN
int tclcommand_part_parse_solvation(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    /* For each fluid component we need 2 constants, one for particle-fluid and one for fluid-particl interaction */
    double solvation[2*LB_COMPONENTS];
    int ii;
    *change = 2*LB_COMPONENTS;
    if (argc < 2*LB_COMPONENTS) {
      Tcl_AppendResult(interp, "solvation requires \"", 2*LB_COMPONENTS, "\"  arguments", (char *) NULL);
      return TCL_ERROR;
    }

    /* set mass */
    for(ii=0;ii<2*LB_COMPONENTS;++ii){
       if (! ARG_IS_D(ii,solvation[ii]))
         return TCL_ERROR;
    }

    if (set_particle_solvation(part_num, solvation) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}


#endif


#ifdef MASS
int tclcommand_part_parse_mass(Tcl_Interp *interp, int argc, char **argv,
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

#ifdef ROTATION_PER_PARTICLE
int tclcommand_part_parse_rotation(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    int rot;

    *change = 1;

    if (argc < 1) {
      Tcl_AppendResult(interp, "rotation requires 1 argument", (char *) NULL);
      return TCL_ERROR;
    }

    /* set rotation flag */
    if (! ARG0_IS_I(rot))
      return TCL_ERROR;

    if (set_particle_rotation(part_num, rot) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}
#endif


#ifdef  ROTATIONAL_INERTIA
int tclcommand_part_parse_rotational_inertia(Tcl_Interp *interp, int argc, char **argv,
				  int part_num, int * change)
{
  double rinertia[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "rotational inertia requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }

  /* set rotational inertia */
  if (! ARG_IS_D(0, rinertia[0]) || ! ARG_IS_D(1, rinertia[1]) || ! ARG_IS_D(2, rinertia[2]))
    return TCL_ERROR;

  if (set_particle_rotational_inertia(part_num, rinertia) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}
#endif

#ifdef DIPOLES
int tclcommand_part_parse_dipm(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_dip(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double dip[3];
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

  /* convenience error message, dipm is not used otherwise. */
  dipm = dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2];
  if (dipm < ROUND_ERROR_PREC) {
    Tcl_AppendResult(interp, "cannot set dipole with zero length", (char *)NULL);
    return TCL_ERROR;
  }

  if (set_particle_dip(part_num, dip) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

#endif

#ifdef VIRTUAL_SITES
int tclcommand_part_parse_virtual(Tcl_Interp *interp, int argc, char **argv,
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

    if (set_particle_virtual(part_num, isVirtual) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}
#endif


#ifdef VIRTUAL_SITES_RELATIVE
int part_parse_vs_relative(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    // See particle_data.hpp for explanation of the quantities
    int vs_relative_to;
    double vs_distance;

    // We consume two arguments after the vs_relative:
    *change = 2;

    // Validate input
    if (argc < 2) {
      Tcl_AppendResult(interp, "vs_relative needs the id of the particle to which the virtual site is related and the distnace it should have from that particle as arguments.", (char *) NULL);
      return TCL_ERROR;
    }

    // Get parameters from tcl
    if (! ARG0_IS_I(vs_relative_to))
      return TCL_ERROR;
    
    if (! ARG1_IS_D(vs_distance))
      return TCL_ERROR;

    if (set_particle_vs_relative(part_num, vs_relative_to, vs_distance) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}


// Set's up the particle, so that it's relative distance to a spicific other 
// particle P, and the relative orientation between the director of particle P
// and the vector pointing from P to this particle stays constant.
// It sets the virtual, vs_relative_distance, ps_relative_to_particle_id and quaternion attributes
// of the this particle.
int part_parse_vs_relate_to(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    // We consume one argument:
    *change = 1;

    // Validate input
    if (argc < 1) {
      Tcl_AppendResult(interp, "vs_auto_relate_to needs the id of the particle to which this particle should be related.", (char *) NULL);
      return TCL_ERROR;
    }

    int relate_to;

    // Get parameters from tcl
    if (! ARG0_IS_I(relate_to))
      return TCL_ERROR;
    
 return vs_relate_to(part_num, relate_to);
}    


#endif




int tclcommand_part_parse_q(Tcl_Interp *interp, int argc, char **argv,
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

#ifdef LB_ELECTROHYDRODYNAMICS
int tclcommand_part_parse_mu_E(Tcl_Interp *interp, int argc, char **argv,
		 int part_num, int * change)
{
    double mu_E[3];

    *change = 3;

    if (argc < 3) {
      Tcl_AppendResult(interp, "mu_E requires 3 arguments", (char *) NULL);
      return TCL_ERROR;
    }

    /* set mobility */
    if (! ARG_IS_D(0, mu_E[0]))
      return TCL_ERROR;
    if (! ARG_IS_D(1, mu_E[1]))
      return TCL_ERROR;
    if (! ARG_IS_D(2, mu_E[2]))
      return TCL_ERROR;

    if (set_particle_mu_E(part_num, mu_E) == TCL_ERROR) {
      Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

      return TCL_ERROR;
    }

    return TCL_OK;
}
#endif

int tclcommand_part_parse_v(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_f(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_type(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_mol_id(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_quat(Tcl_Interp *interp, int argc, char **argv,
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

/* Internal omega (body frame) gets set from values in the body frame */
int tclcommand_part_parse_omega_body(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double omega[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "omega_body requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set angular velocity */
  if (! ARG_IS_D(0, omega[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, omega[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, omega[2]))
    return TCL_ERROR;

   if (set_particle_omega_body(part_num, omega) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

/* Internal omega (body frame) gets set from values in the lab frame */
int tclcommand_part_parse_omega_lab(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double omega[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "omega_lab requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set angular velocity */
  if (! ARG_IS_D(0, omega[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, omega[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, omega[2]))
    return TCL_ERROR;

   if (set_particle_omega_lab(part_num, omega) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

/* Internal torque (body frame) gets set from values in the body frame */
int tclcommand_part_parse_torque_body(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double torque[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "torque_body requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set torque */
  if (! ARG_IS_D(0, torque[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, torque[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, torque[2]))
    return TCL_ERROR;

  if (set_particle_torque_body(part_num, torque) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}

/* Internal torque (body frame) gets set from values in the lab frame */
int tclcommand_part_parse_torque_lab(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double torque[3];

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "torque_lab requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set torque */
  if (! ARG_IS_D(0, torque[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, torque[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, torque[2]))
    return TCL_ERROR;

  if (set_particle_torque_lab(part_num, torque) == TCL_ERROR) {
   Tcl_AppendResult(interp, "set particle position first", (char *)NULL);

    return TCL_ERROR;
  }

  return TCL_OK;
}
#endif

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
int tclcommand_part_parse_ext_torque(Tcl_Interp *interp, int argc, char **argv,
    int part_num, int * change)
{
  double ext_t[3];
  int ext_flag;

  *change = 3;

  if (argc < 3) {
    Tcl_AppendResult(interp, "ext_torque requires 3 arguments", (char *) NULL);
    return TCL_ERROR;
  }
  /* set external torque */
  if (! ARG_IS_D(0, ext_t[0]))
    return TCL_ERROR;

  if (! ARG_IS_D(1, ext_t[1]))
    return TCL_ERROR;

  if (! ARG_IS_D(2, ext_t[2]))
    return TCL_ERROR;

  if (ext_t[0] == 0 && ext_t[1] == 0 && ext_t[2] == 0)
    ext_flag = 0;
  else
    ext_flag = PARTICLE_EXT_TORQUE;

  if (set_particle_ext_torque(part_num, ext_flag, ext_t) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}
#endif

int tclcommand_part_parse_ext_force(Tcl_Interp *interp, int argc, char **argv,
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

  if (set_particle_ext_force(part_num, ext_flag, ext_f) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tclcommand_part_parse_fix(Tcl_Interp *interp, int argc, char **argv,
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

int tclcommand_part_parse_unfix(Tcl_Interp *interp, int argc, char **argv,
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

#ifdef LANGEVIN_PER_PARTICLE
int part_parse_temp(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double T;

  *change = 1;

  if (argc < 1) {
    Tcl_AppendResult(interp, "temp requires 1 argument", (char *) NULL);
    return TCL_ERROR;
  }
  /* set temperature */
  if (! ARG_IS_D(0, T))
    return TCL_ERROR;

  if (set_particle_temperature(part_num, T) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int part_parse_gamma(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  double gamma;

  *change = 1;

  if (argc < 1) {
    Tcl_AppendResult(interp, "gamma requires 1 argument", (char *) NULL);
    return TCL_ERROR;
  }
  /* set temperature scaling factor */
  if (! ARG_IS_D(0, gamma))
    return TCL_ERROR;

  if (set_particle_gamma(part_num, gamma) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

#endif

int part_parse_gc(Tcl_Interp *interp, int argc, char **argv){

  char buffer[100 + TCL_DOUBLE_SPACE + 3*TCL_INTEGER_SPACE];
  int type;

	if (argc < 1) {
	  Tcl_AppendResult(interp, "gc requires at least particle type as argument\nUsage: part gc {<part_type> | {find | delete} <part_type> }\nThe call with only a part_type will init the array lists for the given type\nfind will return a randomly chosen particle id that has the given type\ndelete consequently removes a random particle of given type", (char *) NULL);
	  return TCL_ERROR;
	} else if (ARG0_IS_S("delete")) {
		//delete particle of certain type
		argc--;
		argv++;
		if ( argc < 1 ){
		  Tcl_AppendResult(interp, "gc delete needs particle type as argument", (char *) NULL);
		  return TCL_ERROR;
		}

    if( !ARG_IS_I(0,type) )
    { 
      sprintf(buffer, "\nWrong type for the <type> parameter");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

		if (type < 0 ) {
		  Tcl_AppendResult(interp, "no negative types", (char *) NULL);
		  return TCL_ERROR;
		}
		if ( Type_array_init ) {
			if (delete_particle_of_type(type) == ES_ERROR ) {
			  Tcl_AppendResult(interp, "no particles with type left", (char *) NULL);
			  return TCL_ERROR;
			}
		} else {
			Tcl_AppendResult(interp, "particle lists not initialized", (char *) NULL);
			return TCL_ERROR;
		}
	} else if (ARG0_IS_S("find")) {
		argc--;
		argv++;
		if ( argc < 1 ){ 
		  Tcl_AppendResult(interp, "gc find needs a particle type as argument", (char *) NULL);
		  return TCL_ERROR;
		}

    if( !ARG_IS_I(0,type) )
    { 
      sprintf(buffer, "\nWrong type for the <type> parameter");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

		if (type < 0 ) {
		  Tcl_AppendResult(interp, "no negative types", (char *) NULL);
		  return TCL_ERROR;
		}
		int id;
		if (find_particle_type(type, &id) == ES_ERROR ){
			Tcl_AppendResult(interp, "-1", (char *) NULL);
			return TCL_OK;
		}
		char buffer[32+TCL_INTEGER_SPACE];
		sprintf(buffer, "%d", id);
		Tcl_AppendResult(interp, buffer, (char *) NULL);
	
	} else if ( ARG0_IS_S("status") ) {
		argc--;
		argv++;
		if ( argc < 1 ) {
			Tcl_AppendResult(interp, "gc status need type as argument", (char *) NULL);
			return TCL_ERROR;
		}

    if( !ARG_IS_I(0,type) )
    { 
      sprintf(buffer, "\nWrong type for the <type> parameter");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

		if ( type < 0 ) {
		  Tcl_AppendResult(interp, "no negative types", (char *) NULL);
		  return TCL_ERROR;
		}
		if ( type_array!=(TypeList *) 0 && type_array[Index.type[type]].max_entry!= 0 ) {
			int indexed=0;
			for ( int i=0; i<Type.max_entry; i++) {
				if (type==Type.index[i]) {
					indexed=1;
					break;
				}
			}
			if ( indexed ) {
				char buffer[32 + TCL_INTEGER_SPACE];
				Tcl_AppendResult(interp, "{ ", (char *) NULL);
				for (int i=0; i<type_array[Index.type[type]].max_entry; i++ ) {
					sprintf(buffer, "%d ", type_array[Index.type[type]].id_list[i]);
					Tcl_AppendResult(interp, buffer, (char *) NULL);
				}
				Tcl_AppendResult(interp, " }", (char *) NULL);
			}
		}
		else {
			Tcl_AppendResult(interp, "no list for particle", (char *) NULL);
			return TCL_ERROR;
		}
		return TCL_OK;
	} else if ( ARG0_IS_S("number") ) {
		argc--;
		argv++;
		if ( argc < 1 ) {
			Tcl_AppendResult(interp, "number expects type as argument", (char *) NULL);
			return TCL_ERROR;
		}

    if( !ARG_IS_I(0,type) )
    { 
      sprintf(buffer, "\nWrong type for the <type> parameter");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

		int number;
		if ( type < 0 ) {
		  Tcl_AppendResult(interp, "no negative types", (char *) NULL);
		  return TCL_ERROR;
		}
		if ( number_of_particles_with_type(type, &number) == NOT_INDEXED ) {
			Tcl_AppendResult(interp, "no list for particle", (char *) NULL);
			return TCL_OK;
		}
		char buffer[32 + TCL_INTEGER_SPACE];
		sprintf(buffer, "%d", number);
		Tcl_AppendResult(interp, buffer, (char *) NULL);
		return TCL_OK;
	} else if ( argc == 1 ) {
		// initialize particle array for given type

    if( !ARG_IS_I(0,type) )
    { 
      sprintf(buffer, "\nWrong type for the <type> parameter");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

		if (type < 0 ) {
		  Tcl_AppendResult(interp, "no negative types", (char *) NULL);
		  return TCL_ERROR;
		}
		if (init_type_array(type) == ES_ERROR ) {
			Tcl_AppendResult(interp, "gc init failed", (char *) NULL);
			return TCL_ERROR;
		}

	}
	return TCL_OK;
}

int tclcommand_part_parse_bond(Tcl_Interp *interp, int argc, char **argv,
		    int part_num, int *change)
{
  bool deleteIt = false;
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
    deleteIt = true;
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
      if (change_particle_bond(part_num, bond, deleteIt) != ES_OK) {
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
  if (argc <= 0) {
    if (!deleteIt) {
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
      printf("\nIn particle_data.c\n \n");
      printf("bonded_ia_params[%d].num=%d   ",type_num,bonded_ia_params[type_num].num);
      printf("bonded_ia_params[%d].type=%d  ",type_num,bonded_ia_params[type_num].type);
      printf("n_partners=%d  ",n_partners);
      printf("argc=%d\n",argc);
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
  if (change_particle_bond(part_num, bond, deleteIt) != TCL_OK) {
    Tcl_AppendResult(interp, "bond to delete did not exist", (char *)NULL);
    free(bond);
    return TCL_ERROR;
  }
  free(bond);

  *change += (1 + n_partners);

  return TCL_OK;
}

#ifdef EXCLUSIONS
int tclcommand_part_parse_exclusion(Tcl_Interp *interp, int argc, char **argv,
			 int part_num, int * change)
{
  bool deleteIt = false;
  int partner;

  *change = 0;

  /* check number of arguments */
  if (argc < 1) {
    Tcl_AppendResult(interp, "exclusion requires at least 1 arguments: "
		     "[delete] { <partner 1> <partner 2> ... }", (char *) NULL);
    return TCL_ERROR;
  }

  /* parse away delete eventually */
  if (ARG0_IS_S("delete")) {
    deleteIt = true;
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
    if (change_exclusion(part_num, partner, deleteIt) != TCL_OK) {
      if (deleteIt)
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

int tclcommand_part_parse_cmd(Tcl_Interp *interp, int argc, char **argv,
		   int part_num)
{
  int change = 0, err = TCL_OK;
#ifdef DIPOLES
  int dipm_set = 0;
#endif
#if defined(ROTATION) || defined(DIPOLES)
  int quat_set = 0, dip_set = 0;
#endif

#ifdef ADDITIONAL_CHECKS
  if (!particle_node)
    build_particle_node();
  mpi_bcast_event(CHECK_PARTICLES);
#endif

  if (ARG0_IS_S("print"))
    return tclcommand_part_parse_print(interp, argc-1, argv+1, part_num);

  // various setters
  while (argc > 0) {
    if (ARG0_IS_S("delete"))
      err = tclcommand_part_parse_delete(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("pos"))
      err = tclcommand_part_parse_pos(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("type"))
      err = tclcommand_part_parse_type(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("molecule_id"))
      err = tclcommand_part_parse_mol_id(interp, argc-1, argv+1, part_num, &change);
#ifdef MASS
    else if (ARG0_IS_S("mass"))
      err = tclcommand_part_parse_mass(interp, argc-1, argv+1, part_num, &change);
#endif
#ifdef SHANCHEN 
    else if (ARG0_IS_S("solvation"))
      err = tclcommand_part_parse_solvation(interp, argc-1, argv+1, part_num, &change);
#endif
    else if (ARG0_IS_S("q"))
      err = tclcommand_part_parse_q(interp, argc-1, argv+1, part_num, &change);

#ifdef LB_ELECTROHYDRODYNAMICS
    else if (ARG0_IS_S("mu_E"))
      err = tclcommand_part_parse_mu_E(interp, argc-1, argv+1, part_num, &change);
#endif

    else if (ARG0_IS_S("v"))
      err = tclcommand_part_parse_v(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("f"))
      err = tclcommand_part_parse_f(interp, argc-1, argv+1, part_num, &change);

#ifdef ROTATION

    else if (ARG0_IS_S("quat")) {
      if (dip_set) {
	      Tcl_AppendResult(interp, "(vector) dipole and orientation can not be set at the same time", (char *)NULL);	
        return TCL_ERROR;
      }
      err = tclcommand_part_parse_quat(interp, argc-1, argv+1, part_num, &change);
      quat_set = 1;
    }

    /* Unfortunately a somewhat complex routine is required to make it backwards compatible */
    else if ( ARG0_IS_S("omega") || ARG0_IS_S("omega_body") || ARG0_IS_S("omega_lab") ) 
    {
      if (ARG0_IS_S_EXACT("omega_body"))
      {
        err = tclcommand_part_parse_omega_body(interp, argc-1, argv+1, part_num, &change);
      }
      else if (ARG0_IS_S_EXACT("omega_lab"))
      {
        err = tclcommand_part_parse_omega_lab(interp, argc-1, argv+1, part_num, &change);
      }
      else 
      {
        err = tclcommand_part_parse_omega_body(interp, argc-1, argv+1, part_num, &change);
      }
    }

    /* Unfortunately a somewhat complex routine is required to make it backwards compatible */
    else if ( ARG0_IS_S("torque") || ARG0_IS_S("torque_body") || ARG0_IS_S("torque_lab") ) 
    {
      if (ARG0_IS_S_EXACT("torque_body"))
      {
        err = tclcommand_part_parse_torque_body(interp, argc-1, argv+1, part_num, &change);
      }
      else if (ARG0_IS_S_EXACT("torque_lab"))
      {
        err = tclcommand_part_parse_torque_lab(interp, argc-1, argv+1, part_num, &change);
      }
      else 
      {
        err = tclcommand_part_parse_torque_body(interp, argc-1, argv+1, part_num, &change);
      }
    }

#endif

#ifdef ROTATIONAL_INERTIA
    else if (ARG0_IS_S("rinertia"))
      err = tclcommand_part_parse_rotational_inertia(interp, argc-1, argv+1, part_num, &change);
#endif

#ifdef ROTATION_PER_PARTICLE
    else if (ARG0_IS_S("rotation"))
      err = tclcommand_part_parse_rotation(interp, argc-1, argv+1, part_num, &change);
#endif


#ifdef DIPOLES
    else if (ARG0_IS_S("dip")) {
      if (quat_set) {
	      Tcl_AppendResult(interp, "(vector) dipole and orientation can not be set at the same time", (char *)NULL);	
        return TCL_ERROR;
      }
      if (dipm_set) {
	      Tcl_AppendResult(interp, "(vector) dipole and scalar dipole moment can not be set at the same time", (char *)NULL);	
        return TCL_ERROR;
      }
      err = tclcommand_part_parse_dip(interp, argc-1, argv+1, part_num, &change);
      dip_set = 1;
    }

    else if (ARG0_IS_S("dipm")) {
      if (dip_set) {
	      Tcl_AppendResult(interp, "(vector) dipole and scalar dipole moment can not be set at the same time", (char *)NULL);	
        return TCL_ERROR;
      }
      err = tclcommand_part_parse_dipm(interp, argc-1, argv+1, part_num, &change);
      dipm_set = 1;
    }

#endif

#ifdef VIRTUAL_SITES

    else if (ARG0_IS_S("virtual"))
      err = tclcommand_part_parse_virtual(interp, argc-1, argv+1, part_num, &change);

#endif

#ifdef VIRTUAL_SITES_RELATIVE

    else if (ARG0_IS_S("vs_relative"))
      err = part_parse_vs_relative(interp, argc-1, argv+1, part_num, &change);
    else if (ARG0_IS_S("vs_auto_relate_to"))
      err = part_parse_vs_relate_to(interp, argc-1, argv+1, part_num, &change);
#endif

#ifdef EXTERNAL_FORCES

    else if (ARG0_IS_S("unfix"))
      err = tclcommand_part_parse_unfix(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("fix"))
      err = tclcommand_part_parse_fix(interp, argc-1, argv+1, part_num, &change);

    else if (ARG0_IS_S("ext_force"))
      err = tclcommand_part_parse_ext_force(interp, argc-1, argv+1, part_num, &change);

  #ifdef ROTATION
    else if (ARG0_IS_S("ext_torque"))
      err = tclcommand_part_parse_ext_torque(interp, argc-1, argv+1, part_num, &change);
  #endif
#endif

    else if (ARG0_IS_S("bond"))
      err = tclcommand_part_parse_bond(interp, argc-1, argv+1, part_num, &change);

#ifdef EXCLUSIONS
    else if (ARG0_IS_S("exclude"))
      err = tclcommand_part_parse_exclusion(interp, argc-1, argv+1, part_num, &change);
#endif

#ifdef LANGEVIN_PER_PARTICLE
    else if (ARG0_IS_S("temp"))
      err = part_parse_temp(interp, argc-1, argv+1, part_num, &change);
	else if (ARG0_IS_S("gamma"))
	  err = part_parse_gamma(interp, argc-1, argv+1, part_num, &change);
#endif
	else if (ARG0_IS_S("gc")) { 
	  argc--;
	  argv++;
	  err = part_parse_gc(interp, argc, argv);
	}

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

  return gather_runtime_errors(interp, err);
}


int tclcommand_part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;

  /* if no further arguments are given, print out all stored particles */
  if (argc == 1)
    return tclcommand_part_print_all(interp);
  
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

  else if ( ARG1_IS_S("gc")) {
	 argc-=2;
	 argv+=2;
	 if (part_parse_gc(interp, argc, argv) == TCL_ERROR)
		 return TCL_ERROR;
	 return TCL_OK;
  }
  

  /* only one argument is given */
  if (argc == 2) {

    if (ARG1_IS_I(part_num)) {
      if (tclprint_to_result_Particle(interp, part_num) == TCL_ERROR)
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

  return tclcommand_part_parse_cmd(interp, argc-2, argv+2, part_num);
}
