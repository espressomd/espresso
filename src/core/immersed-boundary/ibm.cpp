/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file ibm.cpp
 *
 * Immersed Boundary methods 
 *
 * includes the fluid membrane interaction functions that distributes the forces on the fluid
 * and sets the velocities of the virtual particles by the interpolation
 *
 * ideally this file should encapsulate all logic related to the immersed_boundary (i.e the pseudo-code of the algorithm)
 * such that, changes in implementation can be carried out here only !.
 *
 */

#include "lb.hpp"
#include "lb-boundaries.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "halo.hpp"
#include "communication.hpp"
#include "ibm.hpp"

#ifdef IMMERSED_BOUNDARY

void lb_ibm_coupling() { 
  
  int i, c, np;
  Cell *cell ;
  Particle *p ;
   
    /* exchange halo regions (for fluid-particle coupling) */
    // halo_communication(&update_halo_comm, (char*)**lbfluid);
      
   
    /* communicate the random numbers */
    // ghost_communicator(&cell_structure.ghost_lbcoupling_comm) ;

#ifdef STRETCHING_FORCE_IMMERSED_BOUNDARY
    ghost_communicator(&cell_structure.ghost_stretching_force_ibm_comm);
#endif

    /* local cells */
    for (c=0;c<local_cells.n;c++) {
    cell = local_cells.cell[c] ;
    p = cell->part ;
    np = cell->n ;

    for (i=0;i<np;i++) {
    if(ifParticleIsVirtual(&p[i])) { 
    couple_trace_to_fluid(&p[i]);
  }
  }}
     for (c=0;c<ghost_cells.n;c++) {
      cell = ghost_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {
        /* for ghost particles we have to check if they lie
        * in the range of the local lattice nodes */
        if (p[i].r.p[0] >= my_left[0]-0.5*lblattice.agrid[0] 
            && p[i].r.p[0] < my_right[0]+0.5*lblattice.agrid[0]
            && p[i].r.p[1] >= my_left[1]-0.5*lblattice.agrid[1] 
            && p[i].r.p[1] < my_right[1]+0.5*lblattice.agrid[1]
            && p[i].r.p[2] >= my_left[2]-0.5*lblattice.agrid[2] 
            && p[i].r.p[2] < my_right[2]+0.5*lblattice.agrid[2]) {

          //Triangles are not subject to viscous coupling, but interact with the fluid via elastic forces
          if(ifParticleIsVirtual(&p[i])) {
            couple_trace_to_fluid(&p[i]);
  } }}}
}

void couple_trace_to_fluid(Particle *p) {
  double *local_f, delta_j[3];
  double force[3], interpolated_u[3];
  index_t node_index[8];
        double delta[6];
  int k,x,y,z;
  
  lblattice.map_position_to_lattice(p->r.p,node_index,delta);
  
  for(k=0; k<3; k++) {
    force[k]=p->f.f[k];
  }
   
   
  //Distribute force among adjacent nodes, just as in viscous coupling
  //if ( p->p.identity == 0 ) printf("unscaled fx = %f\n", force[0]);
  delta_j[0] = force[0]*time_step*lbpar.tau/lbpar.agrid;
  delta_j[1] = force[1]*time_step*lbpar.tau/lbpar.agrid;
  delta_j[2] = force[2]*time_step*lbpar.tau/lbpar.agrid;

  
  // DEBUG
  /*if ( p->p.identity == 0)
  {
    double p_temp[3], v_int[3];
    p_temp[0] = p->r.p[0];
    p_temp[1] = p->r.p[1];
    p_temp[2] = p->r.p[2];
    lb_lbfluid_get_interpolated_velocity_global(p_temp,v_int);
    printf("Before adding force: vx = %e\n", v_int[0]);
  }*/
  // End DEBUG
    
  for (z=0;z<2;z++) {
       for (y=0;y<2;y++) {
            for (x=0;x<2;x++) {
  
        local_f = lbfields[node_index[(z*2+y)*2+x]].force;

        local_f[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[0];
        local_f[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[1];
        local_f[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[2];

            }
        }
    }
  
  // DEBUG
  /*if ( p->p.identity == 0)
  {
    double p_temp[3], v_int[3];
    p_temp[0] = p->r.p[0];
    p_temp[1] = p->r.p[1];
    p_temp[2] = p->r.p[2];
    lb_lbfluid_get_interpolated_velocity_global(p_temp,v_int);
    printf("After adding force: vx = %e\n", v_int[0]);
  }*/
  // End DEBUG  
}

int lb_lbfluid_get_interpolated_velocity_ibm(double* p, double* v, int id) {
  index_t node_index[8], index;
  double delta[6];
  double local_rho, local_j[3], interpolated_u[3];
  double modes[19];
  int x,y,z;
  double *f;

  double lbboundary_mindist, distvec[3];
  double pos[3];

#ifdef LB_BOUNDARIES
  int boundary_no;
  int boundary_flag=-1; // 0 if more than agrid/2 away from the boundary, 1 if 0<dist<agrid/2, 2 if dist <0 

  lbboundary_mindist_position(p, &lbboundary_mindist, distvec, &boundary_no);
  if (lbboundary_mindist>lbpar.agrid/2) {
    boundary_flag=0;
    pos[0]=p[0];
    pos[1]=p[1];
    pos[2]=p[2];
  } else if (lbboundary_mindist > 0 ) {
    boundary_flag=1;
    pos[0]=p[0] - distvec[0]+ distvec[0]/lbboundary_mindist*lbpar.agrid/2.;
    pos[1]=p[1] - distvec[1]+ distvec[1]/lbboundary_mindist*lbpar.agrid/2.;
    pos[2]=p[2] - distvec[2]+ distvec[2]/lbboundary_mindist*lbpar.agrid/2.;
  } else {
    boundary_flag=2;
    v[0]= lb_boundaries[boundary_no].velocity[0]*lbpar.agrid/lbpar.tau;
    v[1]= lb_boundaries[boundary_no].velocity[1]*lbpar.agrid/lbpar.tau;
    v[2]= lb_boundaries[boundary_no].velocity[2]*lbpar.agrid/lbpar.tau;
    return 0; // we can return without interpolating
  }
#else
  pos[0]=p[0];
  pos[1]=p[1];
  pos[2]=p[2];
#endif
  
  /* determine elementary lattice cell surrounding the particle 
     and the relative position of the particle in this cell */ 
  lblattice.map_position_to_lattice(pos,node_index,delta);

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  interpolated_u[0] = interpolated_u[1] = interpolated_u[2] = 0.0 ;

  for (z=0;z<2;z++) {
    for (y=0;y<2;y++) {
      for (x=0;x<2;x++) {
          
        index = node_index[(z*2+y)*2+x];
  f = lbfields[index].force;

        
#ifdef LB_BOUNDARIES
        if (lbfields[index].boundary) {
          local_rho=lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;
          local_j[0] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
          local_j[1] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
          local_j[2] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
        } else {
          lb_calc_modes(index, modes);
          local_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid + modes[0];
          local_j[0] = modes[1] + f[0];
          local_j[1] = modes[2] + f[1];
          local_j[2] = modes[3] + f[2];
        }
#else 
        lb_calc_modes(index, modes);
        local_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid + modes[0];
        local_j[0] = modes[1] + f[0];
        local_j[1] = modes[2] + f[1];
        local_j[2] = modes[3] + f[2];
#endif
  
        interpolated_u[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[0]/(local_rho);
        interpolated_u[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[1]/(local_rho);   
        interpolated_u[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[2]/(local_rho) ;

      }
    }
  }
#ifdef LB_BOUNDARIES
  if (boundary_flag==1) {
    v[0]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[0]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[0];
    v[1]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[1]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[1];
    v[2]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[2]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[2];
  } else {
    v[0] = interpolated_u[0];
    v[1] = interpolated_u[1];
    v[2] = interpolated_u[2];
  }
#else 
    v[0] = interpolated_u[0];
    v[1] = interpolated_u[1];
    v[2] = interpolated_u[2];
#endif
  v[0] *= lbpar.agrid/lbpar.tau;
  v[1] *= lbpar.agrid/lbpar.tau;
  v[2] *= lbpar.agrid/lbpar.tau;
  return 0;
  
}

void force_density_conversion_ibm() { 

  int i = 0;
  for (i=0; i<lblattice.halo_grid_volume; ++i) {

#ifdef EXTERNAL_FORCES
    // unit conversion: force density
    lbfields[i].force[0] = lbpar.ext_force[0]*pow(lbpar.agrid,4)*lbpar.tau*lbpar.tau;
    lbfields[i].force[1] = lbpar.ext_force[1]*pow(lbpar.agrid,4)*lbpar.tau*lbpar.tau;
    lbfields[i].force[2] = lbpar.ext_force[2]*pow(lbpar.agrid,4)*lbpar.tau*lbpar.tau;
#else
    lbfields[i].force[0] = 0.0;
    lbfields[i].force[1] = 0.0;
    lbfields[i].force[2] = 0.0;
    lbfields[i].has_force = 0;
#endif
  }
}

#endif
