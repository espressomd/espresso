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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>

#include "electrokinetics_pdb_parse.hpp"

#include "PdbParser.hpp"

#ifdef EK_BOUNDARIES

/* Replacements for bool variables */
const int pdb_SUCCESS = 0;
const int pdb_ERROR = 1;

float* pdb_charge_lattice = NULL;
int* pdb_boundary_lattice = NULL;

typedef struct {
  float max_x;
  float max_y;
  float max_z;
  float min_x;
  float min_y;
  float min_z;
  float center[3];
} bounding_box;

/* BEGIN CODE */

void galloc(void** ptr, size_t size) {
  if (!*ptr) {
    if (size > 0) {
      *ptr = (void*) Utils::malloc(size);
    }
    else {
      printf("You cannot malloc to size 0\n");
    }
  }
  else {
    if (size > 0) {
      *ptr = (void*) Utils::realloc(*ptr, size);
    }
    else {
      free(*ptr);
    }
  }
}

unsigned int pdb_rhoindex_cartesian2linear(unsigned int x, unsigned int y, unsigned int z) {
  return z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;
}

int print_charge_field(char* filename) {
  FILE* fp;
  if ((fp = fopen(filename,"w")) == NULL) return pdb_ERROR;
  
  if( fp == NULL ) {
    return 1;
  }
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
charge_density\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS charge_density float 1\n\
LOOKUP_TABLE default\n",
  ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
  ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
  ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
  ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z);
  
  for( unsigned int i = 0; i < (ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z); i++ ) {
    fprintf( fp, "%e ", pdb_charge_lattice[i] );
  }
  
  fclose( fp );
  return pdb_SUCCESS;
}

int print_boundary_lattice(char* filename) {
  FILE* fp;
  if ((fp = fopen(filename,"w")) == NULL) return pdb_ERROR;
  
  if( fp == NULL ) {
    return 1;
  }
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
boundary_flag\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS boundary_flag float 1\n\
LOOKUP_TABLE default\n",
  ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
  ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
  ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
  ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z);

  for( unsigned int i = 0; i < (ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z); i++ ) {
    fprintf( fp, "%d ", pdb_boundary_lattice[i] );
  }
  
  fclose( fp );
  return pdb_SUCCESS;
}

int populate_lattice(PdbParser::PdbParser &parser, double scale) {
  /*
   * This routine will populate the lattice using the
   * values read from the pdb and itp files.
   * WARNING: It contains much logic and interpolation stuff!
   */
#ifdef DEBUG
  printf("pdb_n_particles=%u, itp_n_particles=%u, itp_n_parameters=%u\n",atom_data->pdb_n_particles,atom_data->itp_n_particles,atom_data->itp_n_parameters);
#endif
  // TODO: Check if bounding box fits into simbox
  const PdbParser::BoundingBox bbox = parser.calc_bounding_box();
  float center[3];
  
  center[0] = ( bbox.urx + bbox.llx )/2 / scale;
  center[1] = ( bbox.ury + bbox.lly )/2 / scale;
  center[2] = ( bbox.urz + bbox.llz )/2 / scale;

  // calculate the shift of the bounding box
  float shift[3];
  shift[0] = ek_parameters.agrid / 2.0 * ek_parameters.dim_x - center[0];
  shift[1] = ek_parameters.agrid / 2.0 * ek_parameters.dim_y - center[1];
  shift[2] = ek_parameters.agrid / 2.0 * ek_parameters.dim_z - center[2];

#ifdef DEBUG
  printf("agrid=%f, dim_x=%d, dim_y=%d, dim_z=%d\n",ek_parameters.agrid, ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z);
  printf("shift=[%f; %f; %f]\n",shift[0], shift[1], shift[2]);
#endif

  // joining the array
  int lowernode[3];
  float cellpos[3];
  float gridpos;
  float a_x_shifted, a_y_shifted, a_z_shifted;
  float a_x_scaled, a_y_scaled, a_z_scaled;

  for (std::vector<PdbParser::pdb_atom>::const_iterator a = parser.pdb_atoms.begin(); a != parser.pdb_atoms.end(); ++a) {
    PdbParser::itp_atom b = parser.itp_atoms[a->i];
    PdbParser::itp_atomtype c = parser.itp_atomtypes[b.type];
#ifdef DEBUG
    printf("i=%d x=%f y=%f z=%f type=%s charge=%f sigma=%f epsilon=%f\n",a->i,a->x,a->y,a->z,b.type,b.charge,c.sigma,c.epsilon);
#endif

    a_x_scaled = a->x / scale;
    a_y_scaled = a->y / scale;
    a_z_scaled = a->z / scale;

    // Interpolate the charge to the lattice
    gridpos      = (a_x_scaled + shift[0]) / ek_parameters.agrid - 0.5f;
    lowernode[0] = (int) floorf( gridpos );
    cellpos[0]   = gridpos - lowernode[0];
                                                
    gridpos      = (a_y_scaled + shift[1]) / ek_parameters.agrid - 0.5f;
    lowernode[1] = (int) floorf( gridpos );
    cellpos[1]   = gridpos - lowernode[1];
                                                
    gridpos      = (a_z_scaled + shift[2]) / ek_parameters.agrid - 0.5f;
    lowernode[2] = (int) floorf( gridpos );
    cellpos[2]   = gridpos - lowernode[2];
                                                
    lowernode[0] = (lowernode[0] + ek_parameters.dim_x) % ek_parameters.dim_x;
    lowernode[1] = (lowernode[1] + ek_parameters.dim_y) % ek_parameters.dim_y;
    lowernode[2] = (lowernode[2] + ek_parameters.dim_z) % ek_parameters.dim_z;

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],lowernode[1],lowernode[2] )]
      += b.charge * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,lowernode[1],lowernode[2] )]
      += b.charge * cellpos[0] * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],( lowernode[1] + 1 ) % ek_parameters.dim_y,lowernode[2] )]
      += b.charge * ( 1 - cellpos[0] ) * cellpos[1] * ( 1 - cellpos[2] );

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],lowernode[1],( lowernode[2] + 1 ) % ek_parameters.dim_z )]
      += b.charge * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * cellpos[2];

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,( lowernode[1] + 1 ) % ek_parameters.dim_y,lowernode[2] )]
      += b.charge * cellpos[0] * cellpos[1] * ( 1 - cellpos[2] );

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,lowernode[1],( lowernode[2] + 1 ) % ek_parameters.dim_z )]
      += b.charge * cellpos[0] * ( 1 - cellpos[1] ) * cellpos[2];

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],( lowernode[1] + 1 ) % ek_parameters.dim_y,( lowernode[2] + 1 ) % ek_parameters.dim_z )]
      += b.charge * ( 1 - cellpos[0] ) * cellpos[1] * cellpos[2];

    pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,( lowernode[1] + 1 ) % ek_parameters.dim_y,( lowernode[2] + 1 ) % ek_parameters.dim_z )]
      += b.charge * cellpos[0] * cellpos[1] * cellpos[2];
    // Interpolate lennard-jones parameters to boundary
    float r = pow(2,1./6.)*c.sigma * 10 / scale;

    a_x_shifted = (a_x_scaled + shift[0]) / ek_parameters.agrid - 0.5f;
    a_y_shifted = (a_y_scaled + shift[1]) / ek_parameters.agrid - 0.5f;
    a_z_shifted = (a_z_scaled + shift[2]) / ek_parameters.agrid - 0.5f;

    for (float z = a_z_scaled - r; z <= a_z_scaled + r + ek_parameters.agrid; z += ek_parameters.agrid) {
      for (float y = a_y_scaled - r; y <= a_y_scaled + r + ek_parameters.agrid; y += ek_parameters.agrid) {
	for (float x = a_x_scaled - r; x <= a_x_scaled + r + ek_parameters.agrid; x += ek_parameters.agrid) {
	  gridpos      = (x + shift[0]) / ek_parameters.agrid - 0.5f;
	  lowernode[0] = (int) floorf( gridpos );

	  gridpos      = (y + shift[1]) / ek_parameters.agrid - 0.5f;
	  lowernode[1] = (int) floorf( gridpos );

	  gridpos      = (z + shift[2]) / ek_parameters.agrid - 0.5f;
	  lowernode[2] = (int) floorf( gridpos );

	  lowernode[0] = (lowernode[0] + ek_parameters.dim_x) % ek_parameters.dim_x;
	  lowernode[1] = (lowernode[1] + ek_parameters.dim_y) % ek_parameters.dim_y;
	  lowernode[2] = (lowernode[2] + ek_parameters.dim_z) % ek_parameters.dim_z;
#ifdef DEBUG
	  printf("shifted: %f %f %f\n", a_x_shifted, a_y_shifted, a_z_shifted);
	  printf("lowernode: %d %d %d\n", lowernode[0], lowernode[1], lowernode[2]);
	  printf("distance: %f %f %f\n", lowernode[0] - a_x_shifted, lowernode[1] - a_y_shifted, lowernode[2] - a_z_shifted);
	  printf("distance: %f <= %f\n\n", pow(lowernode[0] - a_x_shifted,2) + pow(lowernode[1] - a_y_shifted,2) + pow(lowernode[2] - a_z_shifted,2), pow(r/ek_parameters.agrid,2));
#endif
	  if ( pow(lowernode[0] - a_x_shifted,2) + pow(lowernode[1] - a_y_shifted,2) + pow(lowernode[2] - a_z_shifted,2) <= pow(r/ek_parameters.agrid,2) ) {
	    pdb_boundary_lattice[ek_parameters.dim_y*ek_parameters.dim_x*lowernode[2] + ek_parameters.dim_x*lowernode[1] + lowernode[0]] = 1;
	  }
	}
      }
    }

  }

  return pdb_SUCCESS;
}

int pdb_parse(char* pdb_filename, char* itp_filename, double scale) {
  /*
   * This is the main parsing routine, which is visible to the outside
   * through the header electrokinetics_pdb_parse.h. It doesn't contain any logic and just
   * deploys the input to the subroutines.
   */

  /* BEGIN DEPLOY */
  galloc( (void**) &pdb_charge_lattice, ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z * sizeof(float));
  galloc( (void**) &pdb_boundary_lattice, ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z * sizeof(int));
  for ( unsigned int i = 0; i < ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z; i++ ) {
    pdb_charge_lattice[i] = 0.0;
    pdb_boundary_lattice[i] = 0;
  }

  PdbParser::PdbParser parser;
  if(!parser.parse_file(pdb_filename, itp_filename))
    return pdb_ERROR;

  return populate_lattice(parser, scale);
}

#endif
