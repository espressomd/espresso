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
/** \file galilei.cpp
 *
 */

#include "galilei.hpp"
#include "utils.hpp"
#include "initialize.hpp"
#include "forces.hpp"

galilei_struct gal;

/* Stop the particle motion by setting the 
   velocity of each particle to zero */
void local_kill_particle_motion( int omega ) {
  int c, np, i;
  Particle *part;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      part[i].m.v[0] = 0.0;
      part[i].m.v[1] = 0.0;
      part[i].m.v[2] = 0.0;

      if( omega != 0 ) { 
#ifdef ROTATION
      part[i].m.omega[0] = 0.0;
      part[i].m.omega[1] = 0.0;
      part[i].m.omega[2] = 0.0;
#endif
      }
    }
  }
}

/* Set all the forces acting on the particles
   to zero */
void local_kill_particle_forces( int torque ) {
  int c, np, i;
  Particle *part;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      part[i].f.f[0] = 0.0;
      part[i].f.f[1] = 0.0;
      part[i].f.f[2] = 0.0;

      if( torque != 0 ) { 
#ifdef ROTATION
      part[i].f.torque[0] = 0.0;
      part[i].f.torque[1] = 0.0;
      part[i].f.torque[2] = 0.0;
#endif
      }
    }
  }
}

/* Calculate the CMS of the system */
void local_system_CMS( double *sdata ) {
  int c, np, i;
  Particle *part;
  Cell *cell;
  double x = 0.0, 
         y = 0.0, 
         z = 0.0;
  double ppos[3];
  int img[3];

#ifdef MASS

  double mass = 0.0;
  double M;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      M = part[i].p.mass;
      mass += M;

      memcpy(ppos, part[i].r.p, 3*sizeof(double));
      memcpy(img, part[i].l.i, 3*sizeof(int));
      unfold_position(ppos, img);

      x += M*ppos[0];
      y += M*ppos[1];
      z += M*ppos[2];
    }
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = mass;

#else

  int npart = 0;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      npart++;

      memcpy(ppos, part[i].r.p, 3*sizeof(double));
      memcpy(img, part[i].l.i, 3*sizeof(int));
      unfold_position(ppos, img);

      x += ppos[0];
      y += ppos[1];
      z += ppos[2];
    }
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = (double)npart;

#endif
}

/* Calculate the CMS velocity of the system */
void local_system_CMS_velocity( double *sdata ) {
  int c, np, i;
  Particle *part;
  Cell *cell;
  double x = 0.0, 
         y = 0.0, 
         z = 0.0;

#ifdef MASS

  double mass = 0.0;
  double M;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      M = part[i].p.mass;
      mass += M;

      x += M*part[i].m.v[0];
      y += M*part[i].m.v[1];
      z += M*part[i].m.v[2];
    }
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = mass;

#else

  int npart = 0;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      npart++;

      x += part[i].m.v[0];
      y += part[i].m.v[1];
      z += part[i].m.v[2];
    }
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = (double)npart;

#endif
}

/* Remove the CMS velocity */
void local_galilei_transform( double *sdata ) {
  int c, np, i;
  Particle *part;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      part[i].m.v[0] -= sdata[0];
      part[i].m.v[1] -= sdata[1];
      part[i].m.v[2] -= sdata[2];
    }
  }
}
