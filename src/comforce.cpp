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
/** \file comforce.cpp
 *
 *  Implementation of \ref comforce.hpp
 */
#include "utils.hpp"
#include "comforce.hpp"
#include "cells.hpp"
#include "communication.hpp"

#ifdef COMFORCE
int comforce_set_params(int part_type_a, int part_type_b,
			int flag, int dir, double force, double fratio)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  
  if (!data)
    return 1;

  if (n_nodes > 1)
    return 2;

  data->COMFORCE_flag   = flag;
  data->COMFORCE_dir    = dir;
  data->COMFORCE_force  = force;
  data->COMFORCE_fratio = fratio;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return 0;
}

void calc_comforce()
{
  IA_parameters *ia_params;
  double com0[3], com1[3], MofImatrix[9], diff[3];
  double vect0[3], vect1[3], eva[3], eve[3], fvect[3];
  Particle *p;
  int np;
  Cell *cell;
  
  for (int t0=0; t0<n_particle_types-1; t0++) {
    for (int t1=t0+1; t1<n_particle_types; t1++) {
      ia_params = get_ia_param(t0,t1);
      if (ia_params->COMFORCE_flag == 1) {
        centermass(t0,com0);
        centermass(t1,com1);
        for (int i = 0; i < 3; i++) {
          diff[i]=com1[i]-com0[i];
        }
        momentofinertiamatrix(t0, MofImatrix);
        calc_eigenvalues_3x3(MofImatrix, eva);
        /* perpendicular force */
        if(ia_params->COMFORCE_dir == 1) {
          calc_eigenvector_3x3(MofImatrix,eva[0],eve);
          /*By doing two vector products find radial axis along the target system */
          vector_product(eve,diff,vect0);
          vector_product(vect0,eve,vect1);
          
          /* normalize vect1, return is fvect */
          unit_vector(vect1,fvect);
        } else {
          /* parallel force */
          calc_eigenvector_3x3(MofImatrix,eva[0],fvect);
        }
        
        /* orient it along the com vector */
        if (scalar(fvect,diff) < 0.) {
          for (int i = 0; i < 3; i++) {
            fvect[i] = -fvect[i];
          }
        }
        
        /* Now apply the force */
        for (int c = 0; c < local_cells.n; c++) {
          cell = local_cells.cell[c];
          p  = cell->part;
          np = cell->n;
          for(int i = 0; i < np; i++) {
            if(p[i].p.type==t0) {
      	      for(int j = 0; j < 3; j++) {
                p[i].f.f[j] -= ia_params->COMFORCE_fratio * ia_params->COMFORCE_force * fvect[j];
              }
            }
            if(p[i].p.type==t1) {
      	      for (int j = 0; j < 3; j++) {
                p[i].f.f[j] +=  ia_params->COMFORCE_force * fvect[j];
              }
            }
          }
        }
        /*end of force application */
      }
    }
  }
  
}

#endif
