/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group,
  
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
/** \file lb-boundaries.c
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.h.
 *
 */
#include "utils.h"
#include "constraint.h"
#include "lb-boundaries.h"
#include "lb.h"
#include "interaction_data.h"
#include "communication.h"

#if defined (LB_BOUNDARIES) || defined (LB_BOUNDARIES_GPU)

int n_lb_boundaries = 0;
LB_Boundary *lb_boundaries = NULL;

void lbboundary_mindist_position(double pos[3], double* mindist, double distvec[3], int* no) {
  double vec[3] = {1e100, 1e100, 1e100};
  double dist=1e100;
  *mindist = 1e100;
  int n;

  Particle* p1=0;
  
  for(n=0;n<n_lb_boundaries;n++) {
    switch(lb_boundaries[n].type) {
      case CONSTRAINT_WAL: 
	      calculate_wall_dist(p1, pos, (Particle*) NULL, &lb_boundaries[n].c.wal, &dist, vec); 
        break;
        
      case CONSTRAINT_SPH:
	      calculate_sphere_dist(p1, pos, (Particle*) NULL, &lb_boundaries[n].c.sph, &dist, vec); 
        break;
        
      case CONSTRAINT_CYL: 
	      calculate_cylinder_dist(p1, pos, (Particle*) NULL, &lb_boundaries[n].c.cyl, &dist, vec); 
        break;
        
      case CONSTRAINT_RHOMBOID: 
	      calculate_rhomboid_dist(p1, pos, (Particle*) NULL, &lb_boundaries[n].c.rhomboid, &dist, vec); 
        break;
        
      case CONSTRAINT_PORE: 
	      calculate_pore_dist(p1, pos, (Particle*) NULL, &lb_boundaries[n].c.pore, &dist, vec); 
        break;
    }
    
    if (dist<*mindist) {
      *no=n;
      *mindist=dist;
      distvec[0] = vec[0];
      distvec[1] = vec[1];
      distvec[2] = vec[2];
    } 
  }
}


/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {

  int n, x, y, z;
  char *errtxt;
  double pos[3], dist, dist_tmp=0.0, dist_vec[3];
  
  if (lattice_switch & LATTICE_LB_GPU) {
#if defined (LB_GPU) && defined (LB_BOUNDARIES_GPU)
    int number_of_boundnodes = 0;
    int *host_boundary_node_list= (int*)malloc(sizeof(int));
    int *host_boundary_index_list= (int*)malloc(sizeof(int));
    size_t size_of_index;
    int boundary_number = -1; // the number the boundary will actually belong to.

    for(z=0; z<lbpar_gpu.dim_z; z++) {
      for(y=0; y<lbpar_gpu.dim_y; y++) {
        for (x=0; x<lbpar_gpu.dim_x; x++) {	    
          pos[0] = (x+0.5)*lbpar_gpu.agrid;
          pos[1] = (y+0.5)*lbpar_gpu.agrid;
          pos[2] = (z+0.5)*lbpar_gpu.agrid;
             
          dist = 1e99;

          for (n=0;n<n_lb_boundaries;n++) {
            switch (lb_boundaries[n].type) {
              case LB_BOUNDARY_WAL:
                calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.wal, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_SPH:
                calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.sph, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_CYL:
                calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.cyl, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_RHOMBOID:
                calculate_rhomboid_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.rhomboid, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_POR:
                calculate_pore_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.pore, &dist_tmp, dist_vec);
                break;
                
              default:
                errtxt = runtime_error(128);
                ERROR_SPRINTF(errtxt, "{109 lbboundary type %d not implemented in lb_init_boundaries()\n", lb_boundaries[n].type);
            }
            
            if (dist > dist_tmp || n == 0) {
              dist = dist_tmp;
              boundary_number = n;
            }
          }
          
          if (dist <= 0 && boundary_number > 0 && n_lb_boundaries > 0) {
            size_of_index = (number_of_boundnodes+1)*sizeof(int);
            host_boundary_node_list = realloc(host_boundary_node_list, size_of_index);
            host_boundary_index_list = realloc(host_boundary_index_list, size_of_index);
            host_boundary_node_list[number_of_boundnodes] = x + lbpar_gpu.dim_x*y + lbpar_gpu.dim_x*lbpar_gpu.dim_y*z;
            host_boundary_index_list[number_of_boundnodes] = boundary_number; 
            //printf("boundindex %i: \n", host_boundindex[number_of_boundnodes]);  
            number_of_boundnodes++;  
          }
        }
      }
    }

    /**call of cuda fkt*/
    float* boundary_velocity = malloc(3*n_lb_boundaries*sizeof(float));
    for (n=0; n<n_lb_boundaries; n++) {
      boundary_velocity[3*n+0]=lb_boundaries[n].velocity[0];
      boundary_velocity[3*n+1]=lb_boundaries[n].velocity[1];
      boundary_velocity[3*n+2]=lb_boundaries[n].velocity[2];
    }
    if (n_lb_boundaries)
      lb_init_boundaries_GPU(n_lb_boundaries, number_of_boundnodes, host_boundary_node_list, host_boundary_index_list, boundary_velocity);
    free(boundary_velocity);
    free(host_boundary_node_list);
    free(host_boundary_index_list);
#endif
  } else {
#if defined (LB) && defined (LB_BOUNDARIES)   
    int node_domain_position[3], offset[3];
    int the_boundary=-1;
    map_node_array(this_node, node_domain_position);

    offset[0] = node_domain_position[0]*lblattice.grid[0];
    offset[1] = node_domain_position[1]*lblattice.grid[1];
    offset[2] = node_domain_position[2]*lblattice.grid[2];
    
    for (n=0;n<lblattice.halo_grid_volume;n++) {
      lbfields[n].boundary = 0;
    }
    
    if (lblattice.halo_grid_volume==0)
      return;
    
    for (z=0; z<lblattice.grid[2]+2; z++) {
      for (y=0; y<lblattice.grid[1]+2; y++) {
        for (x=0; x<lblattice.grid[0]+2; x++) {	    
          pos[0] = (offset[0]+(x-0.5))*lblattice.agrid;
          pos[1] = (offset[1]+(y-0.5))*lblattice.agrid;
          pos[2] = (offset[2]+(z-0.5))*lblattice.agrid;
          
          dist = 1e99;

          for (n=0;n<n_lb_boundaries;n++) {
            switch (lb_boundaries[n].type) {
              case LB_BOUNDARY_WAL:
                calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.wal, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_SPH:
                calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.sph, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_CYL:
                calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.cyl, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_RHOMBOID:
                calculate_rhomboid_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.rhomboid, &dist_tmp, dist_vec);
                break;
                
              case LB_BOUNDARY_POR:
                calculate_pore_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.pore, &dist_tmp, dist_vec);
                break;
                
              default:
                errtxt = runtime_error(128);
                ERROR_SPRINTF(errtxt, "{109 lbboundary type %d not implemented in lb_init_boundaries()\n", lb_boundaries[n].type);
            }
            
            if (dist_tmp<dist) {
              dist = dist_tmp;
              the_boundary = n;
            }
          }       
          
    	    if (dist <= 0 && n_lb_boundaries > 0) {
     	      lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = the_boundary+1;
     	      //printf("boundindex %i: \n", get_linear_index(x,y,z,lblattice.halo_grid));   
          }
          else {
            lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 0;
          }
        }
      }
    } 
#endif
  }
}

int lbboundary_get_force(int no, double* f) {
#ifdef LB_BOUNDARIES

  double* forces=malloc(3*n_lb_boundaries*sizeof(double));
  
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_BOUNDARIES_GPU
    lb_gpu_get_boundary_forces(forces);
#else 
    return ES_ERROR;
#endif
  } else { 
    mpi_gather_stats(8, forces, NULL, NULL, NULL);
  }
  
  f[0]=forces[3*no+0]/lbpar.tau/lbpar.tau*lbpar.agrid;
  f[1]=forces[3*no+1]/lbpar.tau/lbpar.tau*lbpar.agrid;
  f[2]=forces[3*no+2]/lbpar.tau/lbpar.tau*lbpar.agrid;
  
  free(forces);
#endif
  return 0;
}

#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */
