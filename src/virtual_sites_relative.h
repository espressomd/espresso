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

