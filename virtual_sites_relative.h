
// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#ifdef VIRTUAL_SITES_RELATIVE

// The following three functions have to be provided by all implementations
// of virtual sites
// Update the vel/pos of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_pos_particle(Particle *);
void update_mol_vel_particle(Particle *);

// Distribute forces that have accumulated on virtual particles to the 
// associated real particles
void distribute_mol_force();

// Setup the virtual_sites_relative properties of a particle so that the given virtaul particle will follow the given real particle
int vs_relate_to(int part_num, int relate_to);

#endif

