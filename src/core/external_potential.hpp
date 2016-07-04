/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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
#ifndef _EXTERNAL_POTENTIAL_HPP
#define _EXTERNAL_POTENTIAL_HPP

#include "lattice.hpp"

#define MAX_FILENAME_SIZE 256

#define EXTERNAL_POTENTIAL_TYPE_TABULATED 0
#define EXTERNAL_POTENTIAL_TYPE_ROD 1

void external_potential_pre_init();

typedef struct {
  char filename[MAX_FILENAME_SIZE];
  Lattice potential;
} ExternalPotentialTabulated;

typedef struct {
  int dummy;
  int type;
  double* scale;
  int n_particle_types;
  ExternalPotentialTabulated tabulated;
  double energy;
} ExternalPotential;

extern ExternalPotential* external_potentials;
extern int n_external_potentials;

int external_potential_tabulated_init(int number, char* filename, int n_particle_types, double* scale);

void external_potential_init_energies();

int generate_external_potential(ExternalPotential** externalPotential);

int external_potential_init(int number, char* filename, int n_particle_types, double* scale);

int external_potential_tabulated_read_potential_file(int number);

void add_external_potential_forces(Particle* p);
void add_external_potential_energy(Particle* p);

int write_local_lattice_to_file(const char* filename_prefix, Lattice* lattice); 


#endif

