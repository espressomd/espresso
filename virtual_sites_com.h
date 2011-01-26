

// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

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


// Unknown
// Analyze the pressure on the molecule level
int parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv);
// Analyze kinetic energy of the molecules
int parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv);
// Sanity checks the positions of virtual sites
int parse_and_check_mol_pos(Tcl_Interp *interp,int argc, char **argv);
// Analze dipole moment on melecular basis
int parse_and_print_dipole_mol(Tcl_Interp *interp,int argc, char **argv);
#endif

