/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

#ifndef VIRTUAL_SITES_H
#define VIRTUAL_SITES_H

#include "particle_data.h"
#include <tcl.h>

/** \file virtual_sites.h
 *  This file contains routine to handle virtual sites
 *  Virtual sites are like particles, but they will be not integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move no-virtual particles
 *  - update virtual sites
 */

#ifdef VIRTUAL_SITES
void update_mol_vel_pos();
void update_mol_vel();
void update_mol_pos();
void update_mol_pos_particle(Particle *);
void update_mol_vel_particle(Particle *);

void distribute_mol_force();

Particle *get_mol_com_particle(Particle *calling_p);

double get_mol_dist(Particle *p1,Particle *p2);
double get_mol_dist_partcfg(Particle *p1,Particle *p2);
MDINLINE int ifParticleIsVirtual(Particle *p){
   if (p->p.isVirtual == 0) {
      return 0;
   }
   else{
      return 1;
   }
}

int update_mol_pos_cfg();
int parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv);
int parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv);
int parse_and_check_mol_pos(Tcl_Interp *interp,int argc, char **argv);
int parse_and_print_dipole_mol(Tcl_Interp *interp,int argc, char **argv);
int parse_pressure_profile_cross_section(Tcl_Interp *interp,int argc, char **argv);

#endif

#endif
