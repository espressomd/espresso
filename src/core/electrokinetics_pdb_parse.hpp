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
/* vim: set ts=8 sts=2 sw=2 et: */
#ifndef _ELECTROKINETICS_PDB_PARSE_HPP
#define _ELECTROKINETICS_PDB_PARSE_HPP

#include "electrokinetics.hpp"

#ifdef ELECTROKINETICS

extern float* pdb_charge_lattice;
extern int* pdb_boundary_lattice;

/* Returns 0/1 if reading the files was successful/unsuccessful */
int pdb_parse(char* pdb_filename, char* itp_filename, double scale);

int print_charge_field(char* filename);

int print_boundary_lattice(char* filename);

#else
/* that is tested for in a number of places, make sure that pdb
   appears disabled if not compiled in.
 */
#define pdb_boundary_lattice 0

#endif

#endif
