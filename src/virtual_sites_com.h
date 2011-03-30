/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#ifndef _VIRTUAL_SITES_COM_H
#define _VIRTUAL_SITES_COM_H

#include "config.h"

#ifdef VIRTUAL_SITES_COM

// The following three functions have to be provided by all implementations
// of virtual sites
// Update the vel/pos of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_pos_particle(Particle *);
void update_mol_vel_particle(Particle *);

// Distribute forces that have accumulated on virtual particles to the 
// associated real particles
void distribute_mol_force();


// Gets the (first) virtual particle of the same molecule as the given (real)
// particle
Particle *get_mol_com_particle(Particle *calling_p);

// Gets the (first) virtual particles in the molecules of the given particles
// and returns their distance
double get_mol_dist(Particle *p1,Particle *p2);

//Distance between molecules in the partcfg data structure
double get_mol_dist_partcfg(Particle *p1,Particle *p2);

// Analyze the pressure on the molecule level
int tclcommand_analyze_parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv);
// Analyze kinetic energy of the molecules
int tclcommand_analyze_parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv);
// Sanity checks the positions of virtual sites
int tclcommand_analyze_parse_and_print_check_mol(Tcl_Interp *interp,int argc, char **argv);
// Analze dipole moment on melecular basis
int tclcommand_analyze_parse_and_print_dipole_mol(Tcl_Interp *interp,int argc, char **argv);
#endif

#endif
