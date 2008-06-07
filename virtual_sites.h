// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

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
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:junghans@mpip-mainz.mpg.de">Ben</a>
 */

#ifdef VIRTUAL_SITES
void update_mol_vel_pos();
void update_mol_vel();
void update_mol_pos();
void distribute_mol_force();

double get_mol_dist(Particle *p1,Particle *p2);
double get_mol_dist_per(Particle *p1,Particle *p2);
double get_mol_dist_cfg(Particle *p1,Particle *p2);

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
#endif

#endif
