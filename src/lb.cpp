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
/** \file lb.cpp
 *
 * Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 * Includes fluctuating LB and coupling to MD particles via frictional 
 * momentum transfer.
 *
 */

#include <mpi.h>
#include <cstdio>
#include "utils.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "interaction_data.hpp"
#include "thermostat.hpp"
#include "lattice.hpp"
#include "halo.hpp"
#include "lb-d3q19.hpp"
#include "lb-boundaries.hpp"
#include "lb.hpp"

#include "cuda_interface.hpp"

// global variable holding the number of fluid components (see global.cpp)
int lb_components = LB_COMPONENTS; 

#ifdef LB

#ifdef ADDITIONAL_CHECKS
static void lb_check_halo_regions();
#endif

/** Flag indicating momentum exchange between particles and fluid */
int transfer_momentum = 0;

/** Struct holding the Lattice Boltzmann parameters */
// LB_Parameters lbpar = { .rho={0.0}, .viscosity={0.0}, .bulk_viscosity={-1.0}, .agrid=-1.0, .tau=-1.0, .friction={0.0}, .ext_force={ 0.0, 0.0, 0.0},.rho_lb_units={0.},.gamma_odd={0.}, .gamma_even={0.} };
LB_Parameters lbpar = {
  // rho
  {0.0},
  // viscosity
  {0.0},
  // bulk_viscosity
  {-1.0},
  // agrid
  -1.0,
  // tau
  -1.0,
  // friction
  {0.0},
  // ext_force
  { 0.0, 0.0, 0.0},
  // rho_lb_units
  {0.},
  // gamma_odd
  {0.},
  // gamma_even
  {0.},
  // resend_halo
  0
};

/** The DnQm model to be used. */
LB_Model lbmodel = { 19, d3q19_lattice, d3q19_coefficients, d3q19_w, NULL, 1./3. };
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * ! MAKE SURE THAT D3Q19 is #undefined WHEN USING OTHER MODELS !
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* doesn't work yet */
#ifndef D3Q19
#error The implementation only works for D3Q19 so far!
#endif

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif

/** The underlying lattice structure */
Lattice lblattice;

/** Pointer to the velocity populations of the fluid nodes */
double **lbfluid[2] = { NULL, NULL };

/** Pointer to the hydrodynamic fields of the fluid nodes */
LB_FluidNode *lbfields = NULL;

/** Communicator for halo exchange between processors */
HaloCommunicator update_halo_comm = { 0, NULL };

/** \name Derived parameters */
/*@{*/
/** Flag indicating whether fluctuations are present. */
static int fluct;

/** relaxation rate of shear modes */
double gamma_shear = 0.0;
/** relaxation rate of bulk modes */
double gamma_bulk = 0.0;
/** relaxation of the odd kinetic modes */
static double gamma_odd  = 0.0;
/** relaxation of the even kinetic modes */
static double gamma_even = 0.0;
/** amplitudes of the fluctuations of the modes */
static double lb_phi[19];
/** amplitude of the fluctuations in the viscous coupling */
static double lb_coupl_pref = 0.0;
/** amplitude of the fluctuations in the viscous coupling with gaussian random numbers */
static double lb_coupl_pref2 = 0.0;
/*@}*/

/** The number of velocities of the LB model.
 * This variable is used for convenience instead of having to type lbmodel.n_veloc everywhere. */
/* TODO: Remove convenience variables */
static int n_veloc;

/** Lattice spacing.
 * This variable is used for convenience instead of having to type lbpar.agrid everywhere. */
static double agrid;

/** Lattice Boltzmann time step
 * This variable is used for convenience instead of having to type lbpar.tau everywhere. */
static double tau;

/** measures the MD time since the last fluid update */
static double fluidstep=0.0;

#ifdef ADDITIONAL_CHECKS
/** counts the random numbers drawn for fluctuating LB and the coupling */
static int rancounter=0;
/** counts the occurences of negative populations due to fluctuations */
static int failcounter=0;
#endif

/***********************************************************************/
#endif

#if defined (LB) || defined (LB_GPU)

/* *********************** C Interface part *************************************/
/* ******************************************************************************/

#ifdef SHANCHEN

int lb_lbfluid_set_shanchen_coupling(double *p_coupling) {
#ifdef LB_GPU
  int ii,jj,n=0;
  switch(LB_COMPONENTS){
        case 1: 
           lbpar_gpu.coupling[0] = (float)p_coupling[0];
           lbpar_gpu.coupling[1] = (float)p_coupling[1];
        break;
        default:
           for(ii=0;ii<LB_COMPONENTS;ii++){
             for(jj=ii;jj<LB_COMPONENTS;jj++){
                  lbpar_gpu.coupling[LB_COMPONENTS*ii+jj] = (float)p_coupling[n]; 
                  lbpar_gpu.coupling[LB_COMPONENTS*jj+ii] = (float)p_coupling[n]; 
                  n++;
             }
           }
        break;

  }
  on_lb_params_change_gpu(LBPAR_COUPLING);
#endif
#ifdef LB
    #error not implemented
#endif
  return 0;
}

int lb_lbfluid_set_mobility(double *p_mobility) {
  int ii;
  for(ii=0;ii<LB_COMPONENTS-1;ii++){
  if ( p_mobility[ii] <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.mobility[ii] = (float)p_mobility[ii];
    on_lb_params_change_gpu(LBPAR_MOBILITY);
#endif
  } else {
#ifdef LB
    #error not implemented
#endif
    }
  }
  return 0;
}

int affinity_set_params(int part_type_a, int part_type_b, double * affinity)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  data->affinity_on=0;
  if (!data) return ES_ERROR;
  for(int ii=0;ii<LB_COMPONENTS;ii++) { 
       if(affinity[ii]<0 || affinity[ii]>1) { return ES_ERROR; }
       data->affinity[ii]= affinity[ii];
       if (data->affinity[ii]>0) data->affinity_on=1;
  }
  
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif   // SHANCHEN

int lb_lbfluid_set_density(double *p_dens) { //
  
  for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( p_dens[ii] <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.rho[ii] = (float)p_dens[ii];
    on_lb_params_change_gpu(LBPAR_DENSITY);
#endif
  } else {
#ifdef LB
    lbpar.rho[ii] = p_dens[ii];
    mpi_bcast_lb_params(LBPAR_DENSITY);
#endif
    }
  }
  return 0;
}
int lb_lbfluid_set_visc(double * p_visc){// 
  
  for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( p_visc[ii] <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.viscosity[ii] = (float)p_visc[ii];
    on_lb_params_change_gpu(LBPAR_VISCOSITY);
#endif
  } else {
#ifdef LB
    lbpar.viscosity[ii] = p_visc[ii];
    mpi_bcast_lb_params(LBPAR_VISCOSITY);
#endif
    }
  }
  return 0;
}
int lb_lbfluid_set_bulk_visc(double *p_bulk_visc){ //
  
 for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( p_bulk_visc[ii] <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.bulk_viscosity[ii] = (float)p_bulk_visc[ii];
    on_lb_params_change_gpu(LBPAR_BULKVISC);
#endif
  } else {
#ifdef LB
    lbpar.bulk_viscosity[ii] = p_bulk_visc[ii];
    mpi_bcast_lb_params(LBPAR_BULKVISC);
#endif
    }
 }
 return 0;
}

int lb_lbfluid_set_gamma_odd(double *p_gamma_odd){//
  
  for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( fabs(p_gamma_odd[ii]) > 1 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_odd[ii] = (float)p_gamma_odd[ii];
    on_lb_params_change_gpu(0);
#endif
  } else {
#ifdef LB
    lbpar.gamma_odd[ii] = gamma_odd = p_gamma_odd[ii];
    mpi_bcast_lb_params(0);
#endif
    }
  }
  return 0;
}

int lb_lbfluid_set_gamma_even(double *p_gamma_even){//

  for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( fabs(p_gamma_even[ii]) > 1 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_even[ii] = (float)p_gamma_even[ii];
    on_lb_params_change_gpu(0);
#endif
  } else {
#ifdef LB
    lbpar.gamma_even[ii] = gamma_even = p_gamma_even[ii];
    mpi_bcast_lb_params(0);
#endif
    }
  }
  return 0;
}

int lb_lbfluid_set_friction(double * p_friction){//

  
  for(int ii=0;ii<LB_COMPONENTS;ii++){
  if ( p_friction[ii] <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.friction[ii] = (float)p_friction[ii];
    on_lb_params_change_gpu(LBPAR_FRICTION);
#endif
  } else {
#ifdef LB
    lbpar.friction[ii] = p_friction[ii];
    mpi_bcast_lb_params(LBPAR_FRICTION);
#endif
    }
  }
  return 0;
}

int lb_lbfluid_set_couple_flag(int couple_flag) { 
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if ( couple_flag != LB_COUPLE_TWO_POINT && couple_flag != LB_COUPLE_THREE_POINT ) {
      return -1;
    }
    lbpar_gpu.lb_couple_switch = couple_flag;
#endif
  } else {
#ifdef LB
    /* Only the two point nearest neighbor coupling is present in the case of the cpu, 
       so just throw an error if something else is tried */
    if ( couple_flag != LB_COUPLE_TWO_POINT ) {
      return -1;
    }   
#endif
  }
  return 0;
}

int lb_lbfluid_set_agrid(double p_agrid){
  if ( p_agrid <= 0) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.agrid = (float)p_agrid;

    lbpar_gpu.dim_x = (unsigned int)rint(box_l[0]/p_agrid);
    lbpar_gpu.dim_y = (unsigned int)rint(box_l[1]/p_agrid);
    lbpar_gpu.dim_z = (unsigned int)rint(box_l[2]/p_agrid);
    unsigned int tmp[3];
    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;
    /* sanity checks */
    int dir;
    for (dir=0;dir<3;dir++) {
    /* check if box_l is compatible with lattice spacing */
      if (fabs(box_l[dir]-tmp[dir]*p_agrid) > ROUND_ERROR_PREC) {
        char *errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{097 Lattice spacing p_agrid=%f is incompatible with box_l[%i]=%f, factor=%d err= %g} ",p_agrid,dir,box_l[dir],tmp[dir],fabs(box_l[dir]-tmp[dir]*p_agrid));
      }
    }
    lbpar_gpu.number_of_nodes = lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
    on_lb_params_change_gpu(LBPAR_AGRID);
#endif
  } else {
#ifdef LB
    lbpar.agrid = agrid = p_agrid;
    mpi_bcast_lb_params(LBPAR_AGRID);
#endif
    }
  return 0;
}



int lb_lbfluid_set_tau(double p_tau){
  if ( p_tau <= 0 ) {
    return -1;
  }
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.tau = (float)p_tau;
    on_lb_params_change_gpu(0);
#endif
  } else {
#ifdef LB
    lbpar.tau = p_tau;
    mpi_bcast_lb_params(0);
#endif
    }
  return 0;
}
#ifdef SHANCHEN
int lb_lbfluid_set_remove_momentum(void){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.remove_momentum = 1;
    on_lb_params_change_gpu(0);
#endif
  } else {
#ifdef LB
     return -1;
#endif
    }
  return 0;
}
#endif

int lb_lbfluid_set_ext_force(int component, double p_fx, double p_fy, double p_fz) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    /* external force density is stored in MD units */
    lbpar_gpu.ext_force[3*component+0] = (float)p_fx;
    lbpar_gpu.ext_force[3*component+1] = (float)p_fy;
    lbpar_gpu.ext_force[3*component+2] = (float)p_fz;
    lbpar_gpu.external_force = 1;
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
#endif
  } else { 
#ifdef LB 
    lbpar.ext_force[0] = p_fx;
    lbpar.ext_force[1] = p_fy;
    lbpar.ext_force[2] = p_fz;
    mpi_bcast_lb_params(LBPAR_EXTFORCE);
#endif
    }
  return 0;
}

int lb_lbfluid_get_density(double* p_dens){ 
// TODO: implement all lb_lbfluid_get_*() correctly for all components!!
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_dens = lbpar_gpu.rho[0];  
#endif
  } else {
#ifdef LB
    *p_dens = lbpar.rho[0];
#endif
    }
  return 0;
}
int lb_lbfluid_get_visc(double* p_visc){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_visc = lbpar_gpu.viscosity[0];
#endif
  } else {
#ifdef LB
    *p_visc = lbpar.viscosity[0];
#endif
  }
  return 0;
}

int lb_lbfluid_get_bulk_visc(double* p_bulk_visc){ 
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_bulk_visc = lbpar_gpu.bulk_viscosity[0];
#endif
  } else {
#ifdef LB
    *p_bulk_visc = lbpar.bulk_viscosity[0];
#endif
  }
  return 0;
}

int lb_lbfluid_get_gamma_odd(double* p_gamma_odd){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_gamma_odd = lbpar_gpu.gamma_odd[0];
#endif
  } else {
#ifdef LB
    *p_gamma_odd = lbpar.gamma_odd[0];
#endif
  }  
return 0;
}

int lb_lbfluid_get_gamma_even(double* p_gamma_even){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_gamma_even = lbpar_gpu.gamma_even[0];
#endif
  } else {
#ifdef LB
    *p_gamma_even = lbpar.gamma_even[0];
#endif
  }
  return 0;
}

int lb_lbfluid_get_agrid(double* p_agrid){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_agrid = lbpar_gpu.agrid;
#endif
  } else {
#ifdef LB
    *p_agrid = lbpar.agrid;
#endif
    }
  return 0;
}
int lb_lbfluid_get_ext_force(double* p_fx, double* p_fy, double* p_fz){
#ifdef SHANCHEN
  exit(printf("Not implemented yet (%s:%d) ",__FILE__,__LINE__));
#endif
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_fx = lbpar_gpu.ext_force[0];
    *p_fy = lbpar_gpu.ext_force[1];
    *p_fz = lbpar_gpu.ext_force[2];
#endif
  } else {
#ifdef LB
    *p_fx = lbpar.ext_force[0];
    *p_fy = lbpar.ext_force[1];
    *p_fz = lbpar.ext_force[2];
#endif
  }
  return 0;
}

int lb_lbfluid_print_vtk_boundary(char* filename) {
  FILE* fp = fopen(filename, "w");
  
  if(fp == NULL)
    return 1;
  
  if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
    unsigned int* bound_array;
    bound_array = (unsigned int*) malloc(lbpar_gpu.number_of_nodes*sizeof(unsigned int));
    lb_get_boundary_flags_GPU(bound_array);
    
    int j;	
    /** print of the calculated phys values */
    fprintf(fp, "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %u\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n", 
            lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z, 
            lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, 
            lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid, 
            lbpar_gpu.number_of_nodes);
    for(j=0; j<int(lbpar_gpu.number_of_nodes); ++j){
      /** print of the calculated phys values */
      fprintf(fp, "%d \n", bound_array[j]);
    }     
    free(bound_array);
#endif
  } else {	
#ifdef LB
    int pos[3];
    int boundary;
    int gridsize[3];
    
    gridsize[0] = box_l[0] / lbpar.agrid;
    gridsize[1] = box_l[1] / lbpar.agrid;
    gridsize[2] = box_l[2] / lbpar.agrid;
    
    fprintf(fp, "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n", 
            gridsize[0], gridsize[1], gridsize[2], 
            lblattice.agrid[0]*0.5, lblattice.agrid[1]*0.5, lblattice.agrid[2]*0.5, 
            lblattice.agrid[0], lblattice.agrid[1], lblattice.agrid[2], 
            gridsize[0]*gridsize[1]*gridsize[2]);
    
    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++)
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++)
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          lb_lbnode_get_boundary(pos, &boundary);
          fprintf(fp, "%d \n", boundary);
        }
#endif
				}
	  fclose(fp);				
	return 0;
}

int lb_lbfluid_print_vtk_density(char** filename) {
  int ii;
  
  for(ii=0;ii<LB_COMPONENTS;++ii) {
    FILE* fp = fopen(filename[ii], "w");
    
    if(fp == NULL) {
      perror("lb_lbfluid_print_vtk_density");
      return 1;
    }
    
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU

     int j;	
     size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
     host_values = (LB_rho_v_pi_gpu*)malloc(size_of_values);
     lb_get_values_GPU(host_values);
     
     fprintf(fp, "# vtk DataFile Version 2.0\nlbfluid_gpu\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %u\nSCALARS density float 1\nLOOKUP_TABLE default\n", 
             lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z, 
             lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, 
             lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.number_of_nodes);
     
     for(j=0; j<int(lbpar_gpu.number_of_nodes); ++j){
         /** print the calculated phys values */
         fprintf(fp, "%f\n", host_values[j].rho[ii]);
     }
     
     free(host_values);

#endif
    }
    else {
#ifdef LB
      exit(printf("Not implemented yet (%s:%d) ",__FILE__,__LINE__));
#endif // LB
    }
    
    fclose(fp);
  }
  
  return 0;
}

int lb_lbfluid_print_vtk_velocity(char* filename) {

  FILE* fp = fopen(filename, "w");

  if(fp == NULL)
    return 1;
  
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu*)malloc(size_of_values);
    lb_get_values_GPU(host_values);
    fprintf(fp, "# vtk DataFile Version 2.0\nlbfluid_gpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %u\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n", 
            lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z, 
            lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, lbpar_gpu.agrid*0.5, 
            lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid, 
            lbpar_gpu.number_of_nodes);
    for(int j=0; j<int(lbpar_gpu.number_of_nodes); ++j){
      /** print of the calculated phys values */
      fprintf(fp, "%f %f %f\n", host_values[j].v[0], host_values[j].v[1], host_values[j].v[2]);
    }
    free(host_values);
#endif
  } else {
#ifdef LB
    int pos[3];
    double u[3];
    int gridsize[3];
    
    gridsize[0] = box_l[0] / lbpar.agrid;
    gridsize[1] = box_l[1] / lbpar.agrid;
    gridsize[2] = box_l[2] / lbpar.agrid;
    
    fprintf(fp, "# vtk DataFile Version 2.0\nlbfluid_cpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n", 
            gridsize[0], gridsize[1], gridsize[2], 
            lblattice.agrid[0]*0.5, lblattice.agrid[1]*0.5, lblattice.agrid[2]*0.5, 
            lblattice.agrid[0], lblattice.agrid[1], lblattice.agrid[2], 
            gridsize[0]*gridsize[1]*gridsize[2]);
	
    for(pos[2] = 0; pos[2] < gridsize[2]; pos[2]++)
      for(pos[1] = 0; pos[1] < gridsize[1]; pos[1]++)
        for(pos[0] = 0; pos[0] < gridsize[0]; pos[0]++)
          {
            lb_lbnode_get_u(pos, u);
            fprintf(fp, "%f %f %f\n", u[0], u[1], u[2]);
          }
#endif     
         }
  fclose(fp);	
			
	return 0;
}

int lb_lbfluid_print_boundary(char* filename) {
	  FILE* fp = fopen(filename, "w");
	
	  if(fp == NULL)
	  	return 1;

  if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
    unsigned int* bound_array; 
    bound_array = (unsigned int*) malloc(lbpar_gpu.number_of_nodes*sizeof(unsigned int));
    lb_get_boundary_flags_GPU(bound_array);

    int xyz[3];
    int j;	
    for(j=0; j<int(lbpar_gpu.number_of_nodes); ++j){
      xyz[0] = j%lbpar_gpu.dim_x;
      int k = j/lbpar_gpu.dim_x;
      xyz[1] = k%lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      /** print of the calculated phys values */
      fprintf(fp, "%f %f %f %d\n", (xyz[0]+0.5)*lbpar_gpu.agrid, (xyz[1]+0.5)*lbpar_gpu.agrid, (xyz[2]+0.5)*lbpar_gpu.agrid, bound_array[j]);
    }
    free(bound_array);
#endif  
  } else {
#ifdef LB
    int pos[3];
    int boundary;
    int gridsize[3];
	
    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];	

    for(pos[2] = 0; pos[2] < gridsize[2]; pos[2]++)
      for(pos[1] = 0; pos[1] < gridsize[1]; pos[1]++)
        for(pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          lb_lbnode_get_boundary(pos, &boundary);
          boundary = ( boundary != 0 ? 1 : 0 );
          fprintf(fp, "%f %f %f %d\n", 
                  (pos[0]+0.5)*lblattice.agrid[0], 
                  (pos[1]+0.5)*lblattice.agrid[1], 
                  (pos[2]+0.5)*lblattice.agrid[2], 
                  boundary);
        }
#endif
  }
	
  fclose(fp);
  return 0;
}

int lb_lbfluid_print_velocity(char* filename) {
    FILE* fp = fopen(filename, "w");

	   if(fp == NULL)
		   return 1;

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
#ifdef SHANCHEN
    printf("TODO:adapt for SHANCHEN (%s:%d)\n",__FILE__,__LINE__);
    exit(1);
#endif
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu*)malloc(size_of_values);
    lb_get_values_GPU(host_values);
    int xyz[3];
    int j;	
    for(j=0; j<int(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j%lbpar_gpu.dim_x;
      int k = j/lbpar_gpu.dim_x;
      xyz[1] = k%lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;     
            /** print of the calculated phys values */
      fprintf(fp, "%f %f %f %f %f %f\n", (xyz[0]+0.5)*lbpar_gpu.agrid, (xyz[1]+0.5)*lbpar_gpu.agrid, (xyz[2]+0.5)*lbpar_gpu.agrid, host_values[j].v[0], host_values[j].v[1], host_values[j].v[2]);
    }
    free(host_values);
#endif
  } else {
#ifdef LB
    int pos[3];
    double u[3];
    int gridsize[3];
    
    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];
 	
    for(pos[2] = 0; pos[2] < gridsize[2]; pos[2]++)
      for(pos[1] = 0; pos[1] < gridsize[1]; pos[1]++)
        for(pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
#ifdef SHANCHEN
          printf("SHANCHEN not implemented for the CPU LB\n",__FILE__,__LINE__);
          exit(1);
#endif
          lb_lbnode_get_u(pos, u);
          fprintf(fp, "%f %f %f %f %f %f\n", 
                  (pos[0]+0.5)*lblattice.agrid[0], 
                  (pos[1]+0.5)*lblattice.agrid[1], 
                  (pos[2]+0.5)*lblattice.agrid[2], 
                  u[0], u[1], u[2]);
        }
#endif		
    }
    
    fclose(fp);		
    return 0;
}

int lb_lbfluid_save_checkpoint(char* filename, int binary) {
  if(lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
     FILE* cpfile;
     cpfile=fopen(filename, "w");
     if (!cpfile) {
        return ES_ERROR;
     }
     float* host_checkpoint_vd = (float *) malloc(lbpar_gpu.number_of_nodes * 19 * sizeof(float));
     unsigned int* host_checkpoint_seed = (unsigned int *) malloc(lbpar_gpu.number_of_nodes * sizeof(unsigned int));
     unsigned int* host_checkpoint_boundary = (unsigned int *) malloc(lbpar_gpu.number_of_nodes * sizeof(unsigned int));
     float* host_checkpoint_force = (float *) malloc(lbpar_gpu.number_of_nodes * 3 * sizeof(float));
     lb_save_checkpoint_GPU(host_checkpoint_vd, host_checkpoint_seed, host_checkpoint_boundary, host_checkpoint_force);
     for (int n=0; n<(19*int(lbpar_gpu.number_of_nodes)); n++) {
         fprintf(cpfile, "%.8E \n", host_checkpoint_vd[n]); 
     }
     for (int n=0; n<int(lbpar_gpu.number_of_nodes); n++) {
         fprintf(cpfile, "%u \n", host_checkpoint_seed[n]); 
     }
     for (int n=0; n<int(lbpar_gpu.number_of_nodes); n++) {
         fprintf(cpfile, "%u \n", host_checkpoint_boundary[n]); 
     }
     for (int n=0; n<(3*int(lbpar_gpu.number_of_nodes)); n++) {
         fprintf(cpfile, "%.8E \n", host_checkpoint_force[n]); 
     }
     fclose(cpfile);
     free(host_checkpoint_vd);
     free(host_checkpoint_seed);
     free(host_checkpoint_boundary);
     free(host_checkpoint_force);
     //fprintf(stderr, "LB checkpointing not implemented for GPU\n");
#endif
  }
  else if(lattice_switch & LATTICE_LB) {
#ifdef LB
		FILE* cpfile;
		cpfile=fopen(filename, "w");
		if (!cpfile) {
			return ES_ERROR;
		}
		double pop[19];
		int ind[3];
		
		int gridsize[3];

		gridsize[0] = box_l[0] / lbpar.agrid;
		gridsize[1] = box_l[1] / lbpar.agrid;
		gridsize[2] = box_l[2] / lbpar.agrid;

		for (int i=0; i < gridsize[0]; i++) {
			for (int j=0; j < gridsize[1]; j++) {
				for (int k=0; k < gridsize[2]; k++) {
					ind[0]=i;
					ind[1]=j;
					ind[2]=k;
					lb_lbnode_get_pop(ind, pop);
					if (!binary) {
						for (int n=0; n<19; n++) {
							fprintf(cpfile, "%.16e ", pop[n]); 
						}
						fprintf(cpfile, "\n"); 
					}
					else {
						fwrite(pop, sizeof(double), 19, cpfile);
					}
				} 
			}
		}
		fclose(cpfile);
#endif
	}
  return ES_OK;
}
int lb_lbfluid_load_checkpoint(char* filename, int binary) {
  if(lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    FILE* cpfile;
    cpfile=fopen(filename, "r");
    if (!cpfile) {
      return ES_ERROR;
    }
    float* host_checkpoint_vd = (float *) malloc(lbpar_gpu.number_of_nodes * 19 * sizeof(float));
    unsigned int* host_checkpoint_seed = (unsigned int *) malloc(lbpar_gpu.number_of_nodes * sizeof(unsigned int));
    unsigned int* host_checkpoint_boundary = (unsigned int *) malloc(lbpar_gpu.number_of_nodes * sizeof(unsigned int));
    float* host_checkpoint_force = (float *) malloc(lbpar_gpu.number_of_nodes * 3 * sizeof(float));

    if (!binary) {
      for (int n=0; n<(19*int(lbpar_gpu.number_of_nodes)); n++) {
          fscanf(cpfile, "%f", &host_checkpoint_vd[n]); 
      }
      for (int n=0; n<int(lbpar_gpu.number_of_nodes); n++) {
          fscanf(cpfile, "%u", &host_checkpoint_seed[n]); 
      }
      for (int n=0; n<int(lbpar_gpu.number_of_nodes); n++) {
          fscanf(cpfile, "%u", &host_checkpoint_boundary[n]); 
      }
      for (int n=0; n<(3*int(lbpar_gpu.number_of_nodes)); n++) {
         fscanf(cpfile, "%f", &host_checkpoint_force[n]); 
      }
     lb_load_checkpoint_GPU(host_checkpoint_vd, host_checkpoint_seed, host_checkpoint_boundary, host_checkpoint_force);
    }
    fclose(cpfile);
    free(host_checkpoint_vd);
    free(host_checkpoint_seed);
    free(host_checkpoint_boundary);
    free(host_checkpoint_force);
#endif
  }
  else if(lattice_switch & LATTICE_LB) {
#ifdef LB
  FILE* cpfile;
  cpfile=fopen(filename, "r");
  if (!cpfile) {
    return ES_ERROR;
  }
  double pop[19];
  int ind[3];
   
  int gridsize[3];
  lbpar.resend_halo=1;
  mpi_bcast_lb_params(0);
  gridsize[0] = box_l[0] / lbpar.agrid;
  gridsize[1] = box_l[1] / lbpar.agrid;
  gridsize[2] = box_l[2] / lbpar.agrid;

  for (int i=0; i < gridsize[0]; i++) {
    for (int j=0; j < gridsize[1]; j++) {
      for (int k=0; k < gridsize[2]; k++) {
        ind[0]=i;
        ind[1]=j;
        ind[2]=k;
        if (!binary) {
          if (fscanf(cpfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", &pop[0],&pop[1],&pop[2],&pop[3],&pop[4],&pop[5],&pop[6],&pop[7],&pop[8],&pop[9],&pop[10],&pop[11],&pop[12],&pop[13],&pop[14],&pop[15],&pop[16],&pop[17],&pop[18]) != 19) {
            return ES_ERROR;
          } 
        }
        else {
          if (fread(pop, sizeof(double), 19, cpfile) != 19)
            return ES_ERROR;
        }
        lb_lbnode_set_pop(ind, pop);
      } 
    }
  }
  fclose(cpfile);
//  lbpar.resend_halo=1;
//  mpi_bcast_lb_params(0);
#endif
  }
  else {
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{To load an LB checkpoint one needs to have already initialized the LB fluid with the same grid size.}");
    return ES_ERROR;
  }
  return ES_OK;
}

int lb_lbnode_get_rho(int* ind, double* p_rho){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    static LB_rho_v_pi_gpu *host_print_values=NULL;

    if(host_print_values==NULL) host_print_values = (LB_rho_v_pi_gpu *) malloc(sizeof(LB_rho_v_pi_gpu));
    lb_print_node_GPU(single_nodeindex, host_print_values);
    for(int ii=0;ii<LB_COMPONENTS;ii++) { 
       p_rho[ii] = (double)(host_print_values->rho[ii]);
    }
#endif // LB_GPU
  } else {
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];
    double rho; double j[3]; double pi[6];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
  
    mpi_recv_fluid(node,index,&rho,j,pi);
    // unit conversion
    rho *= 1/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    *p_rho = rho;
#endif
    }
  return 0;
}

int lb_lbnode_get_u(int* ind, double* p_u) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    static LB_rho_v_pi_gpu *host_print_values=NULL;
    if(host_print_values==NULL) host_print_values= (LB_rho_v_pi_gpu *) malloc(sizeof(LB_rho_v_pi_gpu));

    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, host_print_values);
     
    p_u[0] = (double)(host_print_values[0].v[0]);
    p_u[1] = (double)(host_print_values[0].v[1]);
    p_u[2] = (double)(host_print_values[0].v[2]);
#endif
  } else {  
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];
    double rho; double j[3]; double pi[6];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
  
    mpi_recv_fluid(node,index,&rho,j,pi);
    // unit conversion
    p_u[0] = j[0]/rho/tau/lbpar.agrid;
    p_u[1] = j[1]/rho/tau/lbpar.agrid;
    p_u[2] = j[2]/rho/tau/lbpar.agrid;
#endif
    }
  return 0;
}


/** calculates the fluid velocity at a given position of the 
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. This version of the function
 * can be called without the position needing to be on the local processor.
 * Note that this gives a slightly different version then the values used to
 * couple to MD beads when near a wall, see lb_lbfluid_get_interpolated_velocity.
 */

int lb_lbfluid_get_interpolated_velocity_global (double* p, double* v) {
  double local_v[3] = {0, 0, 0}, delta[6]; //velocity field, relative positions to surrounding nodes
  int 	 ind[3] = {0,0,0}, tmpind[3]; //node indices
  int		 x, y, z; //counters

  // convert the position into lower left grid point
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
	map_position_to_lattice_global(p, ind, delta, lbpar_gpu.agrid);
#endif
  } else {  
#ifdef LB
	map_position_to_lattice_global(p, ind, delta,  lbpar.agrid);
#endif
  }

  //set the initial velocity to zero in all directions
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
  for (z=0;z<2;z++) {
    for (y=0;y<2;y++) {
      for (x=0;x<2;x++) {
				//give the index of the neighbouring nodes
        tmpind[0] = ind[0]+x;
        tmpind[1] = ind[1]+y;
        tmpind[2] = ind[2]+z;

	if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
					if (tmpind[0] == int(lbpar_gpu.dim_x)) tmpind[0] =0;
					if (tmpind[1] == int(lbpar_gpu.dim_y)) tmpind[1] =0;
					if (tmpind[2] == int(lbpar_gpu.dim_z)) tmpind[2] =0;
#endif
	} else {  
#ifdef LB
	  if (tmpind[0] == box_l[0]/lbpar.agrid) tmpind[0] =0;
	  if (tmpind[1] == box_l[1]/lbpar.agrid) tmpind[1] =0;
	  if (tmpind[2] == box_l[2]/lbpar.agrid) tmpind[2] =0;
	  
#endif
	}

#ifdef SHANCHEN
        //printf (" %d %d %d %f %f %f\n", tmpind[0], tmpind[1],tmpind[2],v[0], v[1], v[2]);
        printf("TODO:adapt for SHANCHEN (%s:%d)\n",__FILE__,__LINE__);
        exit(1);
#endif

				lb_lbnode_get_u(tmpind, local_v);
				
	v[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_v[0];
        v[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_v[1];	  
        v[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_v[2];

      }
    }
  }

  return 0;
}

int lb_lbnode_get_pi(int* ind, double* p_pi) {

  double p0 = 0;

  lb_lbnode_get_pi_neq(ind, p_pi);

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    for(int ii=0;ii<LB_COMPONENTS;ii++){
         p0 += lbpar_gpu.rho[ii]*lbpar_gpu.agrid*lbpar_gpu.agrid/lbpar_gpu.tau/lbpar_gpu.tau/3.;
    }
#endif
  } else {  
#ifdef LB
    p0 = lbpar.rho[0]*lbpar.agrid*lbpar.agrid/lbpar.tau/lbpar.tau/3.;
#endif
  }

  p_pi[0] += p0;
  p_pi[2] += p0;
  p_pi[5] += p0;

  return 0;
}

int lb_lbnode_get_pi_neq(int* ind, double* p_pi) {

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    static LB_rho_v_pi_gpu *host_print_values=NULL;
    if(host_print_values==NULL) host_print_values= (LB_rho_v_pi_gpu *) malloc(sizeof(LB_rho_v_pi_gpu));	
    
    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, host_print_values);
    for (int i = 0; i<6; i++) {
      p_pi[i]=host_print_values->pi[i];
    }
    return 0;
#endif
  } else {  
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];
    double rho; double j[3]; double pi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
  
    mpi_recv_fluid(node,index,&rho,j,pi);
    // unit conversion // TODO: Check Unit Conversion!
    p_pi[0] = pi[0]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    p_pi[1] = pi[1]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    p_pi[2] = pi[2]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    p_pi[3] = pi[3]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    p_pi[4] = pi[4]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
    p_pi[5] = pi[5]/tau/tau/lbpar.agrid/lbpar.agrid/lbpar.agrid;
#endif
  }
  return 0;
}

int lb_lbnode_get_boundary(int* ind, int* p_boundary) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    unsigned int host_flag;
    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    lb_get_boundary_flag_GPU(single_nodeindex, &host_flag);
    p_boundary[0] = host_flag;
#endif
  } else {  
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
  
    mpi_recv_fluid_boundary_flag(node,index,p_boundary);
#endif
    }
  return 0;
}

#endif  //defined (LB) || defined (LB_GPU)

int lb_lbnode_get_pop(int* ind, double* p_pop) {
  if (lattice_switch & LATTICE_LB_GPU) {
    fprintf(stderr, "Not implemented for GPU\n");
  } else {
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
    mpi_recv_fluid_populations(node, index, p_pop);
#endif
    }
  return 0;
}

int lb_lbnode_set_rho(int* ind, double *p_rho){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_rho[LB_COMPONENTS];
    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    int i;
    for(i=0;i<LB_COMPONENTS;i++){
       host_rho[i]=(float)p_rho[i];
    }
    lb_set_node_rho_GPU(single_nodeindex, host_rho);
#endif
  } else {
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];
    double rho; double j[3]; double pi[6];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);

    mpi_recv_fluid(node,index,&rho,j,pi);
    rho  = (*p_rho)*agrid*agrid*agrid;
    mpi_send_fluid(node,index,rho,j,pi) ;

//  lb_calc_average_rho();
//  lb_reinit_parameters();
#endif
    }
  return 0;
}

int lb_lbnode_set_u(int* ind, double* u){
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_velocity[3];
    host_velocity[0] = (float)u[0]*lbpar_gpu.tau*lbpar_gpu.agrid;
    host_velocity[1] = (float)u[1]*lbpar_gpu.tau*lbpar_gpu.agrid;
    host_velocity[2] = (float)u[2]*lbpar_gpu.tau*lbpar_gpu.agrid;
    int single_nodeindex = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif
  } else {
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];
    double rho; double j[3]; double pi[6];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);

    /* transform to lattice units */

    mpi_recv_fluid(node,index,&rho,j,pi);
    j[0] = rho*u[0]*tau*agrid;
    j[1] = rho*u[1]*tau*agrid;
    j[2] = rho*u[2]*tau*agrid;
    mpi_send_fluid(node,index,rho,j,pi) ;
#endif
    }
  return 0;
}

int lb_lbnode_set_pi(int* ind, double* pi) {
  return -100;
}

int lb_lbnode_set_pi_neq(int* ind, double* pi_neq) {
  return -100;
}

int lb_lbnode_set_pop(int* ind, double* p_pop) {
  if (lattice_switch & LATTICE_LB_GPU) {

    printf("Not implemented in the LB GPU code!\n");

  } else {
#ifdef LB
    index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0]; ind_shifted[1] = ind[1]; ind_shifted[2] = ind[2];
    node = map_lattice_to_node(&lblattice,ind_shifted,grid);
    index = get_linear_index(ind_shifted[0],ind_shifted[1],ind_shifted[2],lblattice.halo_grid);
    mpi_send_fluid_populations(node, index, p_pop);
#endif
    }
  return 0;
}

int lb_lbnode_set_extforce(int* ind, double* f) {
  return -100;
}

#ifdef LB
/********************** The Main LB Part *************************************/

/* Halo communication for push scheme */
static void halo_push_communication() {

  index_t index;
  int x, y, z, count;
  int rnode, snode;
  double *buffer=NULL, *sbuf=NULL, *rbuf=NULL;
  MPI_Status status;

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];

  /***************
   * X direction *
   ***************/
  count = 5*lblattice.halo_grid[1]*lblattice.halo_grid[2];
  sbuf = (double*) malloc(count*sizeof(double));
  rbuf = (double*) malloc(count*sizeof(double));

  /* send to right, recv from left i = 1, 7, 9, 11, 13 */
  snode = node_neighbors[1];
  rnode = node_neighbors[0];

  buffer = sbuf;
  index = get_linear_index(lblattice.grid[0]+1,0,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (y=0; y<lblattice.halo_grid[1]; y++) {     

      buffer[0] = lbfluid[1][1][index];
      buffer[1] = lbfluid[1][7][index];
      buffer[2] = lbfluid[1][9][index];
      buffer[3] = lbfluid[1][11][index];
      buffer[4] = lbfluid[1][13][index];      
      buffer += 5;

      index += yperiod;
    }    
  }
  
  if (node_grid[0] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(1,0,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (y=0; y<lblattice.halo_grid[1]; y++) {

      lbfluid[1][1][index] = buffer[0];
      lbfluid[1][7][index] = buffer[1];
      lbfluid[1][9][index] = buffer[2];
      lbfluid[1][11][index] = buffer[3];
      lbfluid[1][13][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }    
  }

  /* send to left, recv from right i = 2, 8, 10, 12, 14 */
  snode = node_neighbors[0];
  rnode = node_neighbors[1];

  buffer = sbuf;
  index = get_linear_index(0,0,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (y=0; y<lblattice.halo_grid[1]; y++) {     

      buffer[0] = lbfluid[1][2][index];
      buffer[1] = lbfluid[1][8][index];
      buffer[2] = lbfluid[1][10][index];
      buffer[3] = lbfluid[1][12][index];
      buffer[4] = lbfluid[1][14][index];      
      buffer += 5;

      index += yperiod;
    }    
  }

  if (node_grid[0] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(lblattice.grid[0],0,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (y=0; y<lblattice.halo_grid[1]; y++) {     

      lbfluid[1][2][index] = buffer[0];
      lbfluid[1][8][index] = buffer[1];
      lbfluid[1][10][index] = buffer[2];
      lbfluid[1][12][index] = buffer[3];
      lbfluid[1][14][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }    
  }

  /***************
   * Y direction *
   ***************/
  count = 5*lblattice.halo_grid[0]*lblattice.halo_grid[2];
  sbuf = (double*) realloc(sbuf, count*sizeof(double));
  rbuf = (double*) realloc(rbuf, count*sizeof(double));

  /* send to right, recv from left i = 3, 7, 10, 15, 17 */
  snode = node_neighbors[3];
  rnode = node_neighbors[2];

  buffer = sbuf;
  index = get_linear_index(0,lblattice.grid[1]+1,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {

      buffer[0] = lbfluid[1][3][index];
      buffer[1] = lbfluid[1][7][index];
      buffer[2] = lbfluid[1][10][index];
      buffer[3] = lbfluid[1][15][index];
      buffer[4] = lbfluid[1][17][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  if (node_grid[1] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0,1,0,lblattice.halo_grid);  
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {

      lbfluid[1][3][index] = buffer[0];
      lbfluid[1][7][index] = buffer[1];
      lbfluid[1][10][index] = buffer[2];
      lbfluid[1][15][index] = buffer[3];
      lbfluid[1][17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /* send to left, recv from right i = 4, 8, 9, 16, 18 */
  snode = node_neighbors[2];
  rnode = node_neighbors[3];

  buffer = sbuf;
  index = get_linear_index(0,0,0,lblattice.halo_grid);
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {

      buffer[0] = lbfluid[1][4][index];
      buffer[1] = lbfluid[1][8][index];
      buffer[2] = lbfluid[1][9][index];
      buffer[3] = lbfluid[1][16][index];
      buffer[4] = lbfluid[1][18][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  if (node_grid[1] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0,lblattice.grid[1],0,lblattice.halo_grid); 
  for (z=0; z<lblattice.halo_grid[2]; z++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {
      
      lbfluid[1][4][index] = buffer[0];
      lbfluid[1][8][index] = buffer[1];
      lbfluid[1][9][index] = buffer[2];
      lbfluid[1][16][index] = buffer[3];
      lbfluid[1][18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /***************
   * Z direction *
   ***************/
  count = 5*lblattice.halo_grid[0]*lblattice.halo_grid[1];
  sbuf = (double*) realloc(sbuf, count*sizeof(double));
  rbuf = (double*) realloc(rbuf, count*sizeof(double));
  
  /* send to right, recv from left i = 5, 11, 14, 15, 18 */
  snode = node_neighbors[5];
  rnode = node_neighbors[4];

  buffer = sbuf;
  index = get_linear_index(0,0,lblattice.grid[2]+1,lblattice.halo_grid);
  for (y=0; y<lblattice.halo_grid[1]; y++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {
      
      buffer[0] = lbfluid[1][5][index];
      buffer[1] = lbfluid[1][11][index];
      buffer[2] = lbfluid[1][14][index];
      buffer[3] = lbfluid[1][15][index];
      buffer[4] = lbfluid[1][18][index];      
      buffer += 5;

      ++index;
    }
  }

  if (node_grid[2] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0,0,1,lblattice.halo_grid);  
  for (y=0; y<lblattice.halo_grid[1]; y++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {
      
      lbfluid[1][5][index] = buffer[0];
      lbfluid[1][11][index] = buffer[1];
      lbfluid[1][14][index] = buffer[2];
      lbfluid[1][15][index] = buffer[3];
      lbfluid[1][18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  /* send to left, recv from right i = 6, 12, 13, 16, 17 */
  snode = node_neighbors[4];
  rnode = node_neighbors[5];

  buffer = sbuf;
  index = get_linear_index(0,0,0,lblattice.halo_grid);
  for (y=0; y<lblattice.halo_grid[1]; y++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {
      
      buffer[0] = lbfluid[1][6][index];
      buffer[1] = lbfluid[1][12][index];
      buffer[2] = lbfluid[1][13][index];
      buffer[3] = lbfluid[1][16][index];
      buffer[4] = lbfluid[1][17][index];      
      buffer += 5;

      ++index;
    }
  }

  if (node_grid[2] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD,
		 rbuf, count, MPI_DOUBLE, rnode, REQ_HALO_SPREAD,
		 comm_cart, &status);
  } else {
    memcpy(rbuf,sbuf,count*sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0,0,lblattice.grid[2],lblattice.halo_grid);
  for (y=0; y<lblattice.halo_grid[1]; y++) {
    for (x=0; x<lblattice.halo_grid[0]; x++) {
      
      lbfluid[1][6][index] = buffer[0];
      lbfluid[1][12][index] = buffer[1];
      lbfluid[1][13][index] = buffer[2];
      lbfluid[1][16][index] = buffer[3];
      lbfluid[1][17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  free(rbuf);
  free(sbuf);
}

/***********************************************************************/

/** Performs basic sanity checks. */
int lb_sanity_checks() {
  char *errtext;
  int ret = 0;

  if (lbpar.agrid <= 0.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    ret = 1;
  }
  if (lbpar.tau <= 0.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    ret = 1;
  }
  if (lbpar.rho[0] <= 0.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    ret = 1;
  }
  if (lbpar.viscosity[0] <= 0.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    ret = 1;
  }
  if (dd.use_vList && skin>=lbpar.agrid/2.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{104 LB requires either no Verlet lists or that the skin of the verlet list to be less than half of lattice-Boltzmann grid spacing.} ");
    ret = 1;
  }
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{103 LB requires domain-decomposition cellsystem} ");
    ret = -1;
  } 
  else if (dd.use_vList && skin>=lbpar.agrid/2.0) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{104 LB requires either no Verlet lists or that the skin of the verlet list to be less than half of lattice-Boltzmann grid spacing.} ");
    ret = -1;
  }
  
  if (thermo_switch & ~THERMO_LB) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{122 LB must not be used with other thermostats} ");
    ret = 1;
  }
  return ret;
}

/***********************************************************************/

/** (Pre-)allocate memory for data structures */
void lb_pre_init() {
  n_veloc = lbmodel.n_veloc;
  lbfluid[0]    = (double**) malloc(2*lbmodel.n_veloc*sizeof(double *));
  lbfluid[0][0] = (double*) malloc(2*lblattice.halo_grid_volume*lbmodel.n_veloc*sizeof(double));
}

/** (Re-)allocate memory for the fluid and initialize pointers. */
static void lb_realloc_fluid() {
  int i;

  LB_TRACE(printf("reallocating fluid\n"));

  lbfluid[0]    = (double**) realloc(*lbfluid,2*lbmodel.n_veloc*sizeof(double *));
  lbfluid[0][0] = (double*) realloc(**lbfluid,2*lblattice.halo_grid_volume*lbmodel.n_veloc*sizeof(double));
  lbfluid[1]    = (double **)lbfluid[0] + lbmodel.n_veloc;
  lbfluid[1][0] = (double *)lbfluid[0][0] + lblattice.halo_grid_volume*lbmodel.n_veloc;

  for (i=0; i<lbmodel.n_veloc; ++i) {
    lbfluid[0][i] = lbfluid[0][0] + i*lblattice.halo_grid_volume;
    lbfluid[1][i] = lbfluid[1][0] + i*lblattice.halo_grid_volume;
  }

  lbfields = (LB_FluidNode*) realloc(lbfields,lblattice.halo_grid_volume*sizeof(*lbfields));

}

/** Sets up the structures for exchange of the halo regions.
 *  See also \ref halo.cpp */
static void lb_prepare_communication() {
    int i;
    HaloCommunicator comm = { 0, NULL };

    /* since the data layout is a structure of arrays, we have to
     * generate a communication for this structure: first we generate
     * the communication for one of the arrays (the 0-th velocity
     * population), then we replicate this communication for the other
     * velocity indices by constructing appropriate vector
     * datatypes */

    /* prepare the communication for a single velocity */
    prepare_halo_communication(&comm, &lblattice, FIELDTYPE_DOUBLE, MPI_DOUBLE);

    update_halo_comm.num = comm.num;
    update_halo_comm.halo_info = (HaloInfo*) realloc(update_halo_comm.halo_info,comm.num*sizeof(HaloInfo));

    /* replicate the halo structure */
    for (i=0; i<comm.num; i++) {
      HaloInfo *hinfo = &(update_halo_comm.halo_info[i]);

      hinfo->source_node = comm.halo_info[i].source_node;
      hinfo->dest_node   = comm.halo_info[i].dest_node;
      hinfo->s_offset    = comm.halo_info[i].s_offset;
      hinfo->r_offset    = comm.halo_info[i].r_offset;
      hinfo->type        = comm.halo_info[i].type;

      /* generate the vector datatype for the structure of lattices we
       * have to use hvector here because the extent of the subtypes
       * does not span the full lattice and hence we cannot get the
       * correct vskip out of them */

      MPI_Aint lower;
      MPI_Aint extent;
      MPI_Type_get_extent(MPI_DOUBLE, &lower, &extent);      
      MPI_Type_create_hvector(lbmodel.n_veloc, 1,
                              lblattice.halo_grid_volume*extent, 
                              comm.halo_info[i].datatype, &hinfo->datatype);
      MPI_Type_commit(&hinfo->datatype);
      
      halo_create_field_hvector(lbmodel.n_veloc,1,
                                lblattice.halo_grid_volume*sizeof(double),
                                comm.halo_info[i].fieldtype,&hinfo->fieldtype);
    }      

    release_halo_communication(&comm);    
}

/** (Re-)initializes the fluid. */
void lb_reinit_parameters() {
  int i;

  n_veloc = lbmodel.n_veloc;

  agrid   = lbpar.agrid;
  tau     = lbpar.tau;

  if (lbpar.viscosity[0] > 0.0) {
    /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    // unit conversion: viscosity
    gamma_shear = 1. - 2./(6.*lbpar.viscosity[0]*tau/(agrid*agrid)+1.);
  }

  if (lbpar.bulk_viscosity[0] > 0.0) {
    /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    // unit conversion: viscosity
    gamma_bulk = 1. - 2./(9.*lbpar.bulk_viscosity[0]*tau/(agrid*agrid)+1.);
  }
  
  gamma_odd = lbpar.gamma_odd[0];
  gamma_even = lbpar.gamma_even[0];

  double mu = 0.0;

  if (temperature > 0.0) {  /* fluctuating hydrodynamics ? */

    fluct = 1;

    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).
     * Note that the modes are not normalized as in the paper here! */
    mu = temperature/lbmodel.c_sound_sq*tau*tau/(agrid*agrid);
    //mu *= agrid*agrid*agrid;  // Marcello's conjecture
#ifdef D3Q19
    double (*e)[19] = d3q19_modebase;
#else
    double **e = lbmodel.e;
#endif
    for (i=0; i<3; i++) lb_phi[i] = 0.0;
    lb_phi[4] = sqrt(mu*e[19][4]*(1.-SQR(gamma_bulk)));
    for (i=5; i<10; i++) lb_phi[i] = sqrt(mu*e[19][i]*(1.-SQR(gamma_shear)));
    for (i=10; i<n_veloc; i++) lb_phi[i] = sqrt(mu*e[19][i]);
 
    /* lb_coupl_pref is stored in MD units (force)
     * Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * time_step comes from the discretization.
     */
    lb_coupl_pref = sqrt(12.*2.*lbpar.friction[0]*temperature/time_step);
    lb_coupl_pref2 = sqrt(2.*lbpar.friction[0]*temperature/time_step);

  } else {
    /* no fluctuations at zero temperature */
    fluct = 0;
    for (i=0;i<n_veloc;i++) lb_phi[i] = 0.0;
    lb_coupl_pref = 0.0;
    lb_coupl_pref2 = 0.0;
  }

  LB_TRACE(fprintf(stderr,"%d: gamma_shear=%f gamma_bulk=%f shear_fluct=%f bulk_fluct=%f mu=%f, bulkvisc=%f, visc=%f\n",this_node,gamma_shear,gamma_bulk,lb_phi[9],lb_phi[4],mu, lbpar.bulk_viscosity, lbpar.viscosity));

}


/** Resets the forces on the fluid nodes */
void lb_reinit_forces() {
  index_t index;
  for (index=0; index<lblattice.halo_grid_volume; index++) {

#ifdef EXTERNAL_FORCES
    // unit conversion: force density
      lbfields[index].force[0] = lbpar.ext_force[0]*pow(lbpar.agrid,4)*tau*tau;
      lbfields[index].force[1] = lbpar.ext_force[1]*pow(lbpar.agrid,4)*tau*tau;
      lbfields[index].force[2] = lbpar.ext_force[2]*pow(lbpar.agrid,4)*tau*tau;
#else
      lbfields[index].force[0] = 0.0;
      lbfields[index].force[1] = 0.0;
      lbfields[index].force[2] = 0.0;
      lbfields[index].has_force = 0;
#endif

  }
#ifdef LB_BOUNDARIES
  for (int i =0; i<n_lb_boundaries; i++) {
    lb_boundaries[i].force[0]=0.;
    lb_boundaries[i].force[1]=0.;
    lb_boundaries[i].force[2]=0.;
  }
#endif

}

/** (Re-)initializes the fluid according to the given value of rho. */
void lb_reinit_fluid() {
    index_t index;

    /* default values for fields in lattice units */
    /* here the conversion to lb units is performed */
    double rho = lbpar.rho[0]*agrid*agrid*agrid;
    double j[3] = { 0., 0., 0. };
//    double pi[6] = { rho*lbmodel.c_sound_sq, 0., rho*lbmodel.c_sound_sq, 0., 0., rho*lbmodel.c_sound_sq };
double pi[6] = { 0., 0., 0., 0., 0., 0. };

    LB_TRACE(fprintf(stderr, "Initialising the fluid with equilibrium populations\n"););

    for (index=0; index<lblattice.halo_grid_volume; index++) {

      lb_calc_n_from_rho_j_pi(index,rho,j,pi);

      lbfields[index].recalc_fields = 1;
#ifdef LB_BOUNDARIES
      lbfields[index].boundary = 0;
#endif

    }

    lbpar.resend_halo = 0;
#ifdef LB_BOUNDARIES
    lb_init_boundaries();
#endif
}

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init() {

  LB_TRACE(printf("Begin initialzing fluid on CPU\n"));

  if (lbpar.agrid <= 0.0) {
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set when initializing fluid} ");
  }
  if (check_runtime_errors()) return;

  double temp_agrid[3];
  double temp_offset[3];
  for (int i =0; i<3; i++) {
    temp_agrid[i]=lbpar.agrid;
    temp_offset[i]=0.5;
  }

  /* initialize the local lattice domain */
  init_lattice(&lblattice, temp_agrid, temp_offset, 1, 0);  
//  int init_lattice(Lattice *lattice, double* agrid, double* offset, int halo_size, size_t element_size);

  if (check_runtime_errors()) return;

  /* allocate memory for data structures */
  lb_realloc_fluid();

  /* prepare the halo communication */
  lb_prepare_communication();

  /* initialize derived parameters */
  lb_reinit_parameters();

  /* setup the initial particle velocity distribution */
  lb_reinit_fluid();

  /* setup the external forces */
  lb_reinit_forces();

  LB_TRACE(printf("Initialzing fluid on CPU successful\n"));
}

/** Release the fluid. */
void lb_release_fluid() {
  free(lbfluid[0][0]);
  free(lbfluid[0]);
  free(lbfields);
}

/** Release fluid and communication. */
void lb_release() {
  
  lb_release_fluid();

  release_halo_communication(&update_halo_comm);

}

/***********************************************************************/
/** \name Mapping between hydrodynamic fields and particle populations */
/***********************************************************************/
/*@{*/

void lb_calc_n_from_rho_j_pi(const index_t index, const double rho, const double *j, double *pi) {

  int i;
  double local_rho, local_j[3], local_pi[6], trace;
  const double avg_rho = lbpar.rho[0]*agrid*agrid*agrid;

  local_rho  = rho;

  local_j[0] = j[0];
  local_j[1] = j[1];
  local_j[2] = j[2];

  for (i=0; i<6; i++) 
    local_pi[i] = pi[i];

  trace = local_pi[0] + local_pi[2] + local_pi[5];

#ifdef D3Q19

  double rho_times_coeff;
  double tmp1,tmp2;

  /* update the q=0 sublattice */
  lbfluid[0][0][index] = 1./3. * (local_rho-avg_rho) - 1./2.*trace;

  /* update the q=1 sublattice */
  rho_times_coeff = 1./18. * (local_rho-avg_rho);

  lbfluid[0][1][index] = rho_times_coeff + 1./6.*local_j[0] + 1./4.*local_pi[0] - 1./12.*trace;
  lbfluid[0][2][index] = rho_times_coeff - 1./6.*local_j[0] + 1./4.*local_pi[0] - 1./12.*trace;
  lbfluid[0][3][index] = rho_times_coeff + 1./6.*local_j[1] + 1./4.*local_pi[2] - 1./12.*trace;
  lbfluid[0][4][index] = rho_times_coeff - 1./6.*local_j[1] + 1./4.*local_pi[2] - 1./12.*trace;
  lbfluid[0][5][index] = rho_times_coeff + 1./6.*local_j[2] + 1./4.*local_pi[5] - 1./12.*trace;
  lbfluid[0][6][index] = rho_times_coeff - 1./6.*local_j[2] + 1./4.*local_pi[5] - 1./12.*trace;

  /* update the q=2 sublattice */
  rho_times_coeff = 1./36. * (local_rho-avg_rho);

  tmp1 = local_pi[0] + local_pi[2];
  tmp2 = 2.0*local_pi[1];

  lbfluid[0][7][index]  = rho_times_coeff + 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][8][index]  = rho_times_coeff - 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][9][index]  = rho_times_coeff + 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
  lbfluid[0][10][index] = rho_times_coeff - 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

  tmp1 = local_pi[0] + local_pi[5];
  tmp2 = 2.0*local_pi[3];

  lbfluid[0][11][index] = rho_times_coeff + 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][12][index] = rho_times_coeff - 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][13][index] = rho_times_coeff + 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
  lbfluid[0][14][index] = rho_times_coeff - 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

  tmp1 = local_pi[2] + local_pi[5];
  tmp2 = 2.0*local_pi[4];

  lbfluid[0][15][index] = rho_times_coeff + 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][16][index] = rho_times_coeff - 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
  lbfluid[0][17][index] = rho_times_coeff + 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
  lbfluid[0][18][index] = rho_times_coeff - 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

#else

  int i;
  double tmp=0.0;
  double (*c)[3] = lbmodel.c;
  double (*coeff)[4] = lbmodel.coeff;

  for (i=0;i<n_veloc;i++) {

    tmp = local_pi[0]*c[i][0]*c[i][0]
      + (2.0*local_pi[1]*c[i][0]+local_pi[2]*c[i][1])*c[i][1]
      + (2.0*(local_pi[3]*c[i][0]+local_pi[4]*c[i][1])+local_pi[5]*c[i][2])*c[i][2];

    lbfluid[0][i][index] =  coeff[i][0] * (local_rho-avg_rho);
    lbfluid[0][i][index] += coeff[i][1] * scalar(local_j,c[i]);
    lbfluid[0][i][index] += coeff[i][2] * tmp;
    lbfluid[0][i][index] += coeff[i][3] * trace;

  }

#endif
}
  
/*@}*/


/** Calculation of hydrodynamic modes */
void lb_calc_modes(index_t index, double *mode) {

#ifdef D3Q19
  double n0, n1p, n1m, n2p, n2m, n3p, n3m, n4p, n4m, n5p, n5m, n6p, n6m, n7p, n7m, n8p, n8m, n9p, n9m;

  n0  = lbfluid[0][0][index];
  n1p = lbfluid[0][1][index] + lbfluid[0][2][index];
  n1m = lbfluid[0][1][index] - lbfluid[0][2][index];
  n2p = lbfluid[0][3][index] + lbfluid[0][4][index];
  n2m = lbfluid[0][3][index] - lbfluid[0][4][index];
  n3p = lbfluid[0][5][index] + lbfluid[0][6][index];
  n3m = lbfluid[0][5][index] - lbfluid[0][6][index];
  n4p = lbfluid[0][7][index] + lbfluid[0][8][index];
  n4m = lbfluid[0][7][index] - lbfluid[0][8][index];
  n5p = lbfluid[0][9][index] + lbfluid[0][10][index];
  n5m = lbfluid[0][9][index] - lbfluid[0][10][index];
  n6p = lbfluid[0][11][index] + lbfluid[0][12][index];
  n6m = lbfluid[0][11][index] - lbfluid[0][12][index];
  n7p = lbfluid[0][13][index] + lbfluid[0][14][index];
  n7m = lbfluid[0][13][index] - lbfluid[0][14][index];
  n8p = lbfluid[0][15][index] + lbfluid[0][16][index];
  n8m = lbfluid[0][15][index] - lbfluid[0][16][index];
  n9p = lbfluid[0][17][index] + lbfluid[0][18][index];
  n9m = lbfluid[0][17][index] - lbfluid[0][18][index];
//  printf("n: ");
//  for (i=0; i<19; i++)
//    printf("%f ", lbfluid[1][i][index]);
//  printf("\n");
  
  /* mass mode */
  mode[0] = n0 + n1p + n2p + n3p + n4p + n5p + n6p + n7p + n8p + n9p;
  
  /* momentum modes */
  mode[1] = n1m + n4m + n5m + n6m + n7m;
  mode[2] = n2m + n4m - n5m + n8m + n9m;
  mode[3] = n3m + n6m - n7m + n8m - n9m;

  /* stress modes */
  mode[4] = -n0 + n4p + n5p + n6p + n7p + n8p + n9p;
  mode[5] = n1p - n2p + n6p + n7p - n8p - n9p;
  mode[6] = n1p + n2p - n6p - n7p - n8p - n9p - 2.*(n3p - n4p - n5p);
  mode[7] = n4p - n5p;
  mode[8] = n6p - n7p;
  mode[9] = n8p - n9p;

#ifndef OLD_FLUCT
  /* kinetic modes */
  mode[10] = -2.*n1m + n4m + n5m + n6m + n7m;
  mode[11] = -2.*n2m + n4m - n5m + n8m + n9m;
  mode[12] = -2.*n3m + n6m - n7m + n8m - n9m;
  mode[13] = n4m + n5m - n6m - n7m;
  mode[14] = n4m - n5m - n8m - n9m;
  mode[15] = n6m - n7m - n8m + n9m;
  mode[16] = n0 + n4p + n5p + n6p + n7p + n8p + n9p 
             - 2.*(n1p + n2p + n3p);
  mode[17] = - n1p + n2p + n6p + n7p - n8p - n9p;
  mode[18] = - n1p - n2p -n6p - n7p - n8p - n9p
             + 2.*(n3p + n4p + n5p);
#endif

#else
  int i, j;
  for (i=0; i<n_veloc; i++) {
    mode[i] = 0.0;
    for (j=0; j<n_veloc; j++) {
      mode[i] += lbmodel.e[i][j]*lbfluid[0][i][index];
    }
  }
#endif

}

/** Streaming and calculation of modes (pull scheme) */
inline void lb_pull_calc_modes(index_t index, double *mode) {

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];

  double n[19];
  n[0]  = lbfluid[0][0][index];
  n[1]  = lbfluid[0][1][index-1];
  n[2]  = lbfluid[0][2][index+1];
  n[3]  = lbfluid[0][3][index-yperiod];
  n[4]  = lbfluid[0][4][index+yperiod];
  n[5]  = lbfluid[0][5][index-zperiod];
  n[6]  = lbfluid[0][6][index+zperiod];
  n[7]  = lbfluid[0][7][index-(1+yperiod)];
  n[8]  = lbfluid[0][8][index+(1+yperiod)];
  n[9]  = lbfluid[0][9][index-(1-yperiod)];
  n[10] = lbfluid[0][10][index+(1-yperiod)];
  n[11] = lbfluid[0][11][index-(1+zperiod)];
  n[12] = lbfluid[0][12][index+(1+zperiod)];
  n[13] = lbfluid[0][13][index-(1-zperiod)];
  n[14] = lbfluid[0][14][index+(1-zperiod)];
  n[15] = lbfluid[0][15][index-(yperiod+zperiod)];
  n[16] = lbfluid[0][16][index+(yperiod+zperiod)];
  n[17] = lbfluid[0][17][index-(yperiod-zperiod)];
  n[18] = lbfluid[0][18][index+(yperiod-zperiod)];

#ifdef D3Q19
  /* mass mode */
  mode[ 0] =   n[ 0] + n[ 1] + n[ 2] + n[ 3] + n[4] + n[5] + n[6]
             + n[ 7] + n[ 8] + n[ 9] + n[10]
             + n[11] + n[12] + n[13] + n[14]
             + n[15] + n[16] + n[17] + n[18];

  /* momentum modes */
  mode[ 1] =   n[ 1] - n[ 2] 
             + n[ 7] - n[ 8] + n[ 9] - n[10] + n[11] - n[12] + n[13] - n[14];
  mode[ 2] =   n[ 3] - n[ 4]
             + n[ 7] - n[ 8] - n[ 9] + n[10] + n[15] - n[16] + n[17] - n[18];
  mode[ 3] =   n[ 5] - n[ 6]
             + n[11] - n[12] - n[13] + n[14] + n[15] - n[16] - n[17] + n[18];

  /* stress modes */
  mode[ 4] = - n[ 0] 
             + n[ 7] + n[ 8] + n[ 9] + n[10] 
             + n[11] + n[12] + n[13] + n[14] 
             + n[15] + n[16] + n[17] + n[18];
  mode[ 5] =   n[ 1] + n[ 2] - n[ 3] - n[4]
             + n[11] + n[12] + n[13] + n[14] - n[15] - n[16] - n[17] - n[18];
  mode[ 6] =   n[ 1] + n[ 2] + n[ 3] + n[ 4] 
             - n[11] - n[12] - n[13] - n[14] - n[15] - n[16] - n[17] - n[18]
             - 2.*(n[5] + n[6] - n[7] - n[8] - n[9] - n[10]);
  mode[ 7] =   n[ 7] + n[ 8] - n[ 9] - n[10];
  mode[ 8] =   n[11] + n[12] - n[13] - n[14];
  mode[ 9] =   n[15] + n[16] - n[17] - n[18];

  /* kinetic modes */
  mode[10] = 2.*(n[2] - n[1]) 
             + n[7] - n[8] + n[9] - n[10] + n[11] - n[12] + n[13] - n[14];
  mode[11] = 2.*(n[4] - n[3])
             + n[7] - n[8] - n[9] + n[10] + n[15] - n[16] + n[17] - n[18];
  mode[12] = 2.*(n[6] - n[5])
             + n[11] - n[12] - n[13] + n[14] + n[15] - n[16] - n[17] + n[18];
  mode[13] =   n[ 7] - n[ 8] + n[ 9] - n[10] - n[11] + n[12] - n[13] + n[14];
  mode[14] =   n[ 7] - n[ 8] - n[ 9] + n[10] - n[15] + n[16] - n[17] + n[18];
  mode[15] =   n[11] - n[12] - n[13] + n[14] - n[15] + n[16] + n[17] - n[18];
  mode[16] =   n[ 0]
             + n[ 7] + n[ 8] + n[ 9] + n[10] 
             + n[11] + n[12] + n[13] + n[14] 
             + n[15] + n[16] + n[17] + n[18]
             - 2.*(n[1] + n[2] + n[3] + n[4] + n[5] + n[6]);
  mode[17] =   n[ 3] + n[ 4] - n[ 1] - n[ 2] 
             + n[11] + n[12] + n[13] + n[14] 
             - n[15] - n[16] - n[17] - n[18];
  mode[18] = - n[ 1] - n[ 2] - n[ 3] - n[ 4] 
             - n[11] - n[12] - n[13] - n[14] - n[15] - n[16] - n[17] - n[18]
             + 2.*(n[5] + n[6] + n[7] + n[8] + n[9] + n[10]);
#else
  int i, j;
  double **e = lbmodel.e;
  for (i=0; i<n_veloc; i++) {
    mode[i] = 0.0;
    for (j=0; j<n_veloc; j++) {
      mode[i] += e[i][j]*n[j];
    }
  }
#endif
}

inline void lb_relax_modes(index_t index, double *mode) {

  double rho, j[3], pi_eq[6];

  /* re-construct the real density 
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = mode[0] + lbpar.rho[0]*agrid*agrid*agrid;

  j[0] = mode[1];
  j[1] = mode[2];
  j[2] = mode[3];

  /* if forces are present, the momentum density is redefined to
   * inlcude one half-step of the force action.  See the
   * Chapman-Enskog expansion in [Ladd & Verberg]. */
#ifndef EXTERNAL_FORCES
  if (lbfields[index].has_force) 
#endif
  {
    j[0] += 0.5*lbfields[index].force[0];
    j[1] += 0.5*lbfields[index].force[1];
    j[2] += 0.5*lbfields[index].force[2];
  }

  /* equilibrium part of the stress modes */
  pi_eq[0] = scalar(j,j)/rho;
  pi_eq[1] = (SQR(j[0])-SQR(j[1]))/rho;
  pi_eq[2] = (scalar(j,j) - 3.0*SQR(j[2]))/rho;
  pi_eq[3] = j[0]*j[1]/rho;
  pi_eq[4] = j[0]*j[2]/rho;
  pi_eq[5] = j[1]*j[2]/rho;

  /* relax the stress modes */  
  mode[4] = pi_eq[0] + gamma_bulk*(mode[4] - pi_eq[0]);
  mode[5] = pi_eq[1] + gamma_shear*(mode[5] - pi_eq[1]);
  mode[6] = pi_eq[2] + gamma_shear*(mode[6] - pi_eq[2]);
  mode[7] = pi_eq[3] + gamma_shear*(mode[7] - pi_eq[3]);
  mode[8] = pi_eq[4] + gamma_shear*(mode[8] - pi_eq[4]);
  mode[9] = pi_eq[5] + gamma_shear*(mode[9] - pi_eq[5]);

#ifndef OLD_FLUCT
  /* relax the ghost modes (project them out) */
  /* ghost modes have no equilibrium part due to orthogonality */
  mode[10] = gamma_odd*mode[10];
  mode[11] = gamma_odd*mode[11];
  mode[12] = gamma_odd*mode[12];
  mode[13] = gamma_odd*mode[13];
  mode[14] = gamma_odd*mode[14];
  mode[15] = gamma_odd*mode[15];
  mode[16] = gamma_even*mode[16];
  mode[17] = gamma_even*mode[17];
  mode[18] = gamma_even*mode[18];
#endif

}

inline void lb_thermalize_modes(index_t index, double *mode) {
    double fluct[6];
#ifdef GAUSSRANDOM
    double rootrho_gauss = sqrt(fabs(mode[0]+lbpar.rho[0]*agrid*agrid*agrid));

    /* stress modes */
    mode[4] += (fluct[0] = rootrho_gauss*lb_phi[4]*gaussian_random());
    mode[5] += (fluct[1] = rootrho_gauss*lb_phi[5]*gaussian_random());
    mode[6] += (fluct[2] = rootrho_gauss*lb_phi[6]*gaussian_random());
    mode[7] += (fluct[3] = rootrho_gauss*lb_phi[7]*gaussian_random());
    mode[8] += (fluct[4] = rootrho_gauss*lb_phi[8]*gaussian_random());
    mode[9] += (fluct[5] = rootrho_gauss*lb_phi[9]*gaussian_random());
    
#ifndef OLD_FLUCT
    /* ghost modes */
    mode[10] += rootrho_gauss*lb_phi[10]*gaussian_random();
    mode[11] += rootrho_gauss*lb_phi[11]*gaussian_random();
    mode[12] += rootrho_gauss*lb_phi[12]*gaussian_random();
    mode[13] += rootrho_gauss*lb_phi[13]*gaussian_random();
    mode[14] += rootrho_gauss*lb_phi[14]*gaussian_random();
    mode[15] += rootrho_gauss*lb_phi[15]*gaussian_random();
    mode[16] += rootrho_gauss*lb_phi[16]*gaussian_random();
    mode[17] += rootrho_gauss*lb_phi[17]*gaussian_random();
    mode[18] += rootrho_gauss*lb_phi[18]*gaussian_random();
#endif

#elif defined (GAUSSRANDOMCUT)
    double rootrho_gauss = sqrt(fabs(mode[0]+lbpar.rho[0]*agrid*agrid*agrid));

    /* stress modes */
    mode[4] += (fluct[0] = rootrho_gauss*lb_phi[4]*gaussian_random_cut());
    mode[5] += (fluct[1] = rootrho_gauss*lb_phi[5]*gaussian_random_cut());
    mode[6] += (fluct[2] = rootrho_gauss*lb_phi[6]*gaussian_random_cut());
    mode[7] += (fluct[3] = rootrho_gauss*lb_phi[7]*gaussian_random_cut());
    mode[8] += (fluct[4] = rootrho_gauss*lb_phi[8]*gaussian_random_cut());
    mode[9] += (fluct[5] = rootrho_gauss*lb_phi[9]*gaussian_random_cut());
    
#ifndef OLD_FLUCT
    /* ghost modes */
    mode[10] += rootrho_gauss*lb_phi[10]*gaussian_random_cut();
    mode[11] += rootrho_gauss*lb_phi[11]*gaussian_random_cut();
    mode[12] += rootrho_gauss*lb_phi[12]*gaussian_random_cut();
    mode[13] += rootrho_gauss*lb_phi[13]*gaussian_random_cut();
    mode[14] += rootrho_gauss*lb_phi[14]*gaussian_random_cut();
    mode[15] += rootrho_gauss*lb_phi[15]*gaussian_random_cut();
    mode[16] += rootrho_gauss*lb_phi[16]*gaussian_random_cut();
    mode[17] += rootrho_gauss*lb_phi[17]*gaussian_random_cut();
    mode[18] += rootrho_gauss*lb_phi[18]*gaussian_random_cut();
#endif
    
#elif defined (FLATNOISE)
    double rootrho = sqrt(fabs(12.0*(mode[0]+lbpar.rho[0]*agrid*agrid*agrid)));

    /* stress modes */
    mode[4] += (fluct[0] = rootrho*lb_phi[4]*(d_random()-0.5));
    mode[5] += (fluct[1] = rootrho*lb_phi[5]*(d_random()-0.5));
    mode[6] += (fluct[2] = rootrho*lb_phi[6]*(d_random()-0.5));
    mode[7] += (fluct[3] = rootrho*lb_phi[7]*(d_random()-0.5));
    mode[8] += (fluct[4] = rootrho*lb_phi[8]*(d_random()-0.5));
    mode[9] += (fluct[5] = rootrho*lb_phi[9]*(d_random()-0.5));
    
#ifndef OLD_FLUCT
    /* ghost modes */
    mode[10] += rootrho*lb_phi[10]*(d_random()-0.5);
    mode[11] += rootrho*lb_phi[11]*(d_random()-0.5);
    mode[12] += rootrho*lb_phi[12]*(d_random()-0.5);
    mode[13] += rootrho*lb_phi[13]*(d_random()-0.5);
    mode[14] += rootrho*lb_phi[14]*(d_random()-0.5);
    mode[15] += rootrho*lb_phi[15]*(d_random()-0.5);
    mode[16] += rootrho*lb_phi[16]*(d_random()-0.5);
    mode[17] += rootrho*lb_phi[17]*(d_random()-0.5);
    mode[18] += rootrho*lb_phi[18]*(d_random()-0.5);
#endif
#else
#error No noise type defined for the CPU LB
#endif //GAUSSRANDOM
    
#ifdef ADDITIONAL_CHECKS
    rancounter += 15;
#endif
}

inline void lb_apply_forces(index_t index, double* mode) {

  double rho, *f, u[3], C[6];
  
  f = lbfields[index].force;

  rho = mode[0] + lbpar.rho[0]*agrid*agrid*agrid;

  /* hydrodynamic momentum density is redefined when external forces present */
  u[0] = (mode[1] + 0.5*f[0])/rho;
  u[1] = (mode[2] + 0.5*f[1])/rho;
  u[2] = (mode[3] + 0.5*f[2])/rho;

  C[0] = (1.+gamma_bulk)*u[0]*f[0] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[2] = (1.+gamma_bulk)*u[1]*f[1] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[5] = (1.+gamma_bulk)*u[2]*f[2] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[1] = 1./2.*(1.+gamma_shear)*(u[0]*f[1]+u[1]*f[0]);
  C[3] = 1./2.*(1.+gamma_shear)*(u[0]*f[2]+u[2]*f[0]);
  C[4] = 1./2.*(1.+gamma_shear)*(u[1]*f[2]+u[2]*f[1]);

  /* update momentum modes */
  mode[1] += f[0];
  mode[2] += f[1];
  mode[3] += f[2];

  /* update stress modes */
  mode[4] += C[0] + C[2] + C[5];
//  mode[5] += 2.*C[0] - C[2] - C[5];
//  mode[6] += C[2] - C[5];
//ulf vorschlag: mode[6] += C[2] + C[5] - 2.*C[0];

  mode[5] += C[0] - C[2];
  mode[6] += C[0] + C[2] - 2.*C[5];
  mode[7] += C[1];
  mode[8] += C[3];
  mode[9] += C[4];

  /* reset force */
#ifdef EXTERNAL_FORCES
  // unit conversion: force density
  lbfields[index].force[0] = lbpar.ext_force[0]*pow(lbpar.agrid,4)*tau*tau;
  lbfields[index].force[1] = lbpar.ext_force[1]*pow(lbpar.agrid,4)*tau*tau;
  lbfields[index].force[2] = lbpar.ext_force[2]*pow(lbpar.agrid,4)*tau*tau;
#else
  lbfields[index].force[0] = 0.0;
  lbfields[index].force[1] = 0.0;
  lbfields[index].force[2] = 0.0;
  lbfields[index].has_force = 0;
#endif

}

inline void lb_calc_n_from_modes(index_t index, double *mode) {

  int i;
  double *w = lbmodel.w;

#ifdef D3Q19
  double (*e)[19] = d3q19_modebase;
  double m[19];

  /* normalization factors enter in the back transformation */
  for (i=0;i<n_veloc;i++) {
    m[i] = (1./e[19][i])*mode[i];
  }

  lbfluid[0][ 0][index] = m[0] - m[4] + m[16];
  lbfluid[0][ 1][index] = m[0] + m[1] + m[5] + m[6] - m[17] - m[18] - 2.*(m[10] + m[16]);
  lbfluid[0][ 2][index] = m[0] - m[1] + m[5] + m[6] - m[17] - m[18] + 2.*(m[10] - m[16]);
  lbfluid[0][ 3][index] = m[0] + m[2] - m[5] + m[6] + m[17] - m[18] - 2.*(m[11] + m[16]);
  lbfluid[0][ 4][index] = m[0] - m[2] - m[5] + m[6] + m[17] - m[18] + 2.*(m[11] - m[16]);
  lbfluid[0][ 5][index] = m[0] + m[3] - 2.*(m[6] + m[12] + m[16] - m[18]);
  lbfluid[0][ 6][index] = m[0] - m[3] - 2.*(m[6] - m[12] + m[16] - m[18]);
  lbfluid[0][ 7][index] = m[0] + m[ 1] + m[ 2] + m[ 4] + 2.*m[6]
        + m[7] + m[10] + m[11] + m[13] + m[14] + m[16] + 2.*m[18];
  lbfluid[0][ 8][index] = m[0] - m[ 1] - m[ 2] + m[ 4] + 2.*m[6]
        + m[7] - m[10] - m[11] - m[13] - m[14] + m[16] + 2.*m[18];
  lbfluid[0][ 9][index] = m[0] + m[ 1] - m[ 2] + m[ 4] + 2.*m[6]
        - m[7] + m[10] - m[11] + m[13] - m[14] + m[16] + 2.*m[18];
  lbfluid[0][10][index] = m[0] - m[ 1] + m[ 2] + m[ 4] + 2.*m[6]
        - m[7] - m[10] + m[11] - m[13] + m[14] + m[16] + 2.*m[18];
  lbfluid[0][11][index] = m[0] + m[ 1] + m[ 3] + m[ 4] + m[ 5] - m[ 6]
        + m[8] + m[10] + m[12] - m[13] + m[15] + m[16] + m[17] - m[18];
  lbfluid[0][12][index] = m[0] - m[ 1] - m[ 3] + m[ 4] + m[ 5] - m[ 6]
        + m[8] - m[10] - m[12] + m[13] - m[15] + m[16] + m[17] - m[18];
  lbfluid[0][13][index] = m[0] + m[ 1] - m[ 3] + m[ 4] + m[ 5] - m[ 6]
        - m[8] + m[10] - m[12] - m[13] - m[15] + m[16] + m[17] - m[18];
  lbfluid[0][14][index] = m[0] - m[ 1] + m[ 3] + m[ 4] + m[ 5] - m[ 6]
        - m[8] - m[10] + m[12] + m[13] + m[15] + m[16] + m[17] - m[18];
  lbfluid[0][15][index] = m[0] + m[ 2] + m[ 3] + m[ 4] - m[ 5] - m[ 6]
        + m[9] + m[11] + m[12] - m[14] - m[15] + m[16] - m[17] - m[18];
  lbfluid[0][16][index] = m[0] - m[ 2] - m[ 3] + m[ 4] - m[ 5] - m[ 6]
        + m[9] - m[11] - m[12] + m[14] + m[15] + m[16] - m[17] - m[18];
  lbfluid[0][17][index] = m[0] + m[ 2] - m[ 3] + m[ 4] - m[ 5] - m[ 6]
        - m[9] + m[11] - m[12] - m[14] + m[15] + m[16] - m[17] - m[18];
  lbfluid[0][18][index] = m[0] - m[ 2] + m[ 3] + m[ 4] - m[ 5] - m[ 6]
        - m[9] - m[11] + m[12] + m[14] - m[15] + m[16] - m[17] - m[18];

  /* weights enter in the back transformation */
  for (i=0;i<n_veloc;i++) {
    lbfluid[0][i][index] *= w[i];
  }

#else
  int j;
  double **e = lbmodel.e;
  for (i=0; i<n_veloc;i++) {
    lbfluid[0][i][index] = 0.0;
    for (j=0;j<n_veloc;j++) {
      lbfluid[0][i][index] += mode[j]*e[j][i]/e[19][j];
    }
    lbfluid[0][i][index] *= w[i];
  }
#endif

}

inline void lb_calc_n_from_modes_push(index_t index, double *m) {
    int i;

#ifdef D3Q19
    int yperiod = lblattice.halo_grid[0];
    int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
    index_t next[19];
    next[0]  = index;
    next[1]  = index + 1;
    next[2]  = index - 1;
    next[3]  = index + yperiod;
    next[4]  = index - yperiod;
    next[5]  = index + zperiod;
    next[6]  = index - zperiod;
    next[7]  = index + (1 + yperiod);
    next[8]  = index - (1 + yperiod);
    next[9]  = index + (1 - yperiod);
    next[10] = index - (1 - yperiod);
    next[11] = index + (1 + zperiod);
    next[12] = index - (1 + zperiod);
    next[13] = index + (1 - zperiod);
    next[14] = index - (1 - zperiod);
    next[15] = index + (yperiod + zperiod);
    next[16] = index - (yperiod + zperiod);
    next[17] = index + (yperiod - zperiod);
    next[18] = index - (yperiod - zperiod);

    /* normalization factors enter in the back transformation */
    for (i=0;i<n_veloc;i++) {
      m[i] = (1./d3q19_modebase[19][i])*m[i];
    }

#ifndef OLD_FLUCT
    lbfluid[1][ 0][next[0]] = m[0] - m[4] + m[16];
    lbfluid[1][ 1][next[1]] = m[0] + m[1] + m[5] + m[6] - m[17] - m[18] - 2.*(m[10] + m[16]);
    lbfluid[1][ 2][next[2]] = m[0] - m[1] + m[5] + m[6] - m[17] - m[18] + 2.*(m[10] - m[16]);
    lbfluid[1][ 3][next[3]] = m[0] + m[2] - m[5] + m[6] + m[17] - m[18] - 2.*(m[11] + m[16]);
    lbfluid[1][ 4][next[4]] = m[0] - m[2] - m[5] + m[6] + m[17] - m[18] + 2.*(m[11] - m[16]);
    lbfluid[1][ 5][next[5]] = m[0] + m[3] - 2.*(m[6] + m[12] + m[16] - m[18]);
    lbfluid[1][ 6][next[6]] = m[0] - m[3] - 2.*(m[6] - m[12] + m[16] - m[18]);
    lbfluid[1][ 7][next[7]] = m[0] + m[ 1] + m[ 2] + m[ 4] + 2.*m[6] + m[7] + m[10] + m[11] + m[13] + m[14] + m[16] + 2.*m[18];
    lbfluid[1][ 8][next[8]] = m[0] - m[ 1] - m[ 2] + m[ 4] + 2.*m[6] + m[7] - m[10] - m[11] - m[13] - m[14] + m[16] + 2.*m[18];
    lbfluid[1][ 9][next[9]] = m[0] + m[ 1] - m[ 2] + m[ 4] + 2.*m[6] - m[7] + m[10] - m[11] + m[13] - m[14] + m[16] + 2.*m[18];
    lbfluid[1][10][next[10]] = m[0] - m[ 1] + m[ 2] + m[ 4] + 2.*m[6] - m[7] - m[10] + m[11] - m[13] + m[14] + m[16] + 2.*m[18];
    lbfluid[1][11][next[11]] = m[0] + m[ 1] + m[ 3] + m[ 4] + m[ 5] - m[ 6] + m[8] + m[10] + m[12] - m[13] + m[15] + m[16] + m[17] - m[18];
    lbfluid[1][12][next[12]] = m[0] - m[ 1] - m[ 3] + m[ 4] + m[ 5] - m[ 6] + m[8] - m[10] - m[12] + m[13] - m[15] + m[16] + m[17] - m[18];
    lbfluid[1][13][next[13]] = m[0] + m[ 1] - m[ 3] + m[ 4] + m[ 5] - m[ 6] - m[8] + m[10] - m[12] - m[13] - m[15] + m[16] + m[17] - m[18];
    lbfluid[1][14][next[14]] = m[0] - m[ 1] + m[ 3] + m[ 4] + m[ 5] - m[ 6] - m[8] - m[10] + m[12] + m[13] + m[15] + m[16] + m[17] - m[18];
    lbfluid[1][15][next[15]] = m[0] + m[ 2] + m[ 3] + m[ 4] - m[ 5] - m[ 6] + m[9] + m[11] + m[12] - m[14] - m[15] + m[16] - m[17] - m[18];
    lbfluid[1][16][next[16]] = m[0] - m[ 2] - m[ 3] + m[ 4] - m[ 5] - m[ 6] + m[9] - m[11] - m[12] + m[14] + m[15] + m[16] - m[17] - m[18];
    lbfluid[1][17][next[17]] = m[0] + m[ 2] - m[ 3] + m[ 4] - m[ 5] - m[ 6] - m[9] + m[11] - m[12] - m[14] + m[15] + m[16] - m[17] - m[18];
    lbfluid[1][18][next[18]] = m[0] - m[ 2] + m[ 3] + m[ 4] - m[ 5] - m[ 6] - m[9] - m[11] + m[12] + m[14] - m[15] + m[16] - m[17] - m[18];
#else
    lbfluid[1][ 0][next[0]] = m[0] - m[4];
    lbfluid[1][ 1][next[1]] = m[0] + m[1] + m[5] + m[6];
    lbfluid[1][ 2][next[2]] = m[0] - m[1] + m[5] + m[6];
    lbfluid[1][ 3][next[3]] = m[0] + m[2] - m[5] + m[6];
    lbfluid[1][ 4][next[4]] = m[0] - m[2] - m[5] + m[6];
    lbfluid[1][ 5][next[5]] = m[0] + m[3] - 2.*m[6];
    lbfluid[1][ 6][next[6]] = m[0] - m[3] - 2.*m[6];
    lbfluid[1][ 7][next[7]] = m[0] + m[1] + m[2] + m[4] + 2.*m[6] + m[7];
    lbfluid[1][ 8][next[8]] = m[0] - m[1] - m[2] + m[4] + 2.*m[6] + m[7];
    lbfluid[1][ 9][next[9]] = m[0] + m[1] - m[2] + m[4] + 2.*m[6] - m[7];
    lbfluid[1][10][next[10]] = m[0] - m[1] + m[2] + m[4] + 2.*m[6] - m[7];
    lbfluid[1][11][next[11]] = m[0] + m[1] + m[3] + m[4] + m[5] - m[6] + m[8];
    lbfluid[1][12][next[12]] = m[0] - m[1] - m[3] + m[4] + m[5] - m[6] + m[8];
    lbfluid[1][13][next[13]] = m[0] + m[1] - m[3] + m[4] + m[5] - m[6] - m[8];
    lbfluid[1][14][next[14]] = m[0] - m[1] + m[3] + m[4] + m[5] - m[6] - m[8];
    lbfluid[1][15][next[15]] = m[0] + m[2] + m[3] + m[4] - m[5] - m[6] + m[9];
    lbfluid[1][16][next[16]] = m[0] - m[2] - m[3] + m[4] - m[5] - m[6] + m[9];
    lbfluid[1][17][next[17]] = m[0] + m[2] - m[3] + m[4] - m[5] - m[6] - m[9];
    lbfluid[1][18][next[18]] = m[0] - m[2] + m[3] + m[4] - m[5] - m[6] - m[9];
#endif

    /* weights enter in the back transformation */
    for (i=0;i<n_veloc;i++) {
      lbfluid[1][i][next[i]] *= lbmodel.w[i];
    }
#else
  int j;
  double **e = lbmodel.e;
  index_t next[n_veloc];
  for (i=0; i<n_veloc;i++) {
    next[i] = get_linear_index(c[i][0],c[i][1],c[i][2],lblattic.halo_grid);
    lbfluid[1][i][next[i]] = 0.0;
    for (j=0;j<n_veloc;j++) {
      lbfluid[1][i][next[i]] += mode[j]*e[j][i]/e[19][j];
    }
    lbfluid[1][i][index] *= w[i];
  }
#endif

}

/* Collisions and streaming (push scheme) */
inline void lb_collide_stream() {
    index_t index;
    int x, y, z;
    double modes[19];

    /* loop over all lattice cells (halo excluded) */
#ifdef LB_BOUNDARIES
    for (int i =0; i < n_lb_boundaries; i++) {
      lb_boundaries[i].force[0]=0.;
      lb_boundaries[i].force[1]=0.;
      lb_boundaries[i].force[2]=0.;
    }
#endif

    index = lblattice.halo_offset;
    for (z=1; z<=lblattice.grid[2]; z++) {
      for (y=1; y<=lblattice.grid[1]; y++) {
	for (x=1; x<=lblattice.grid[0]; x++) {
	  
#ifdef LB_BOUNDARIES
	  if (!lbfields[index].boundary)
#endif
	  {
	
	    /* calculate modes locally */
	    lb_calc_modes(index, modes);

	    /* deterministic collisions */
	    lb_relax_modes(index, modes);

	    /* fluctuating hydrodynamics */
	    if (fluct) lb_thermalize_modes(index, modes);

	    /* apply forces */
#ifdef EXTERNAL_FORCES
	    lb_apply_forces(index, modes);
#else
	    if (lbfields[index].has_force) lb_apply_forces(index, modes);
#endif

	    /* transform back to populations and streaming */
	    lb_calc_n_from_modes_push(index, modes);

	  }
#ifdef LB_BOUNDARIES
	  else {

/*      Here collision in the boundary nodes
 *      can be included, if this is necessary */
//	    lb_boundary_collisions(index, modes);

	  }
#endif

	  ++index; /* next node */
	}

	index += 2; /* skip halo region */
      }
      
      index += 2*lblattice.halo_grid[0]; /* skip halo region */
    }

    /* exchange halo regions */
    halo_push_communication();

#ifdef LB_BOUNDARIES
    /* boundary conditions for links */
    lb_bounce_back();
#endif

   /* swap the pointers for old and new population fields */
    double **tmp;
    tmp = lbfluid[0];
    lbfluid[0] = lbfluid[1];
    lbfluid[1] = tmp;

    /* halo region is invalid after update */
    lbpar.resend_halo = 1;
}

/** Streaming and collisions (pull scheme) */
inline void lb_stream_collide() {
    index_t index;
    int x, y, z;
    double modes[19];

    /* exchange halo regions */
    halo_communication(&update_halo_comm,(char*)**lbfluid);
#ifdef ADDITIONAL_CHECKS
    lb_check_halo_regions();
#endif

    /* loop over all lattice cells (halo excluded) */
    index = lblattice.halo_offset;
    for (z=1; z<=lblattice.grid[2]; z++) {
      for (y=1; y<=lblattice.grid[1]; y++) {
	for (x=1; x<=lblattice.grid[0]; x++) {
	  
	  {

	    /* stream (pull) and calculate modes */
	    lb_pull_calc_modes(index, modes);
  
	    /* deterministic collisions */
	    lb_relax_modes(index, modes);
    
	    /* fluctuating hydrodynamics */
	    if (fluct) lb_thermalize_modes(index, modes);
  
	    /* apply forces */
	    if (lbfields[index].has_force) lb_apply_forces(index, modes);
    
	    /* calculate new particle populations */
	    lb_calc_n_from_modes(index, modes);

	  }

	  ++index; /* next node */
	}

	index += 2; /* skip halo region */
      }
      
      index += 2*lblattice.halo_grid[0]; /* skip halo region */
    }

    /* swap the pointers for old and new population fields */
    //fprintf(stderr,"swapping pointers\n");
    double **tmp = lbfluid[0];
    lbfluid[0] = lbfluid[1];
    lbfluid[1] = tmp;

    /* halo region is invalid after update */
    lbpar.resend_halo = 1;
      
}

/***********************************************************************/
/** \name Update step for the lattice Boltzmann fluid                  */
/***********************************************************************/
/*@{*/
/*@}*/


/** Update the lattice Boltzmann fluid.  
 *
 * This function is called from the integrator. Since the time step
 * for the lattice dynamics can be coarser than the MD time step, we
 * monitor the time since the last lattice update.
 */
void lattice_boltzmann_update() {

  int factor = (int)round(tau/time_step);

  fluidstep += 1;

  if (fluidstep>=factor) {
    fluidstep=0;

#ifdef PULL
    lb_stream_collide();
#else 
    lb_collide_stream();
#endif
  }
  
}

/***********************************************************************/
/** \name Coupling part */
/***********************************************************************/
/*@{*/


/** Coupling of a single particle to viscous fluid with Stokesian friction.
 * 
 * Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
 *
 * @param p          The coupled particle (Input).
 * @param force      Coupling force between particle and fluid (Output).
 */
inline void lb_viscous_coupling(Particle *p, double force[3]) {
  int x,y,z;
  index_t node_index[8];
  double delta[6];
  double *local_f, interpolated_u[3],delta_j[3];

#if 0 // I have no idea what this should be for!
  if(!(p->l.ext_flag & COORD_FIXED(0)) && !(p->l.ext_flag & COORD_FIXED(1)) && !(p->l.ext_flag & COORD_FIXED(2)))
#endif
  
  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));

  /* determine elementary lattice cell surrounding the particle 
     and the relative position of the particle in this cell */ 
  map_position_to_lattice(&lblattice,p->r.p,node_index,delta);
  
  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB delta=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f) pos=(%.3f,%.3f,%.3f)\n",this_node,delta[0],delta[1],delta[2],delta[3],delta[4],delta[5],p->r.p[0],p->r.p[1],p->r.p[2]));

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  lb_lbfluid_get_interpolated_velocity(p->r.p, interpolated_u);

  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB u = (%.16e,%.3e,%.3e) v = (%.16e,%.3e,%.3e)\n",this_node,interpolated_u[0],interpolated_u[1],interpolated_u[2],p->m.v[0],p->m.v[1],p->m.v[2]));

  /* calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units 
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
#ifdef LB_ELECTROHYDRODYNAMICS
  force[0] = - lbpar.friction[0] * (p->m.v[0]/time_step - interpolated_u[0] - p->p.mu_E[0]);
  force[1] = - lbpar.friction[0] * (p->m.v[1]/time_step - interpolated_u[1] - p->p.mu_E[1]);
  force[2] = - lbpar.friction[0] * (p->m.v[2]/time_step - interpolated_u[2] - p->p.mu_E[2]);
#else
  force[0] = - lbpar.friction[0] * (p->m.v[0]/time_step - interpolated_u[0]);
  force[1] = - lbpar.friction[0] * (p->m.v[1]/time_step - interpolated_u[1]);
  force[2] = - lbpar.friction[0] * (p->m.v[2]/time_step - interpolated_u[2]);
#endif


  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_drag = (%.6e,%.3e,%.3e)\n",this_node,force[0],force[1],force[2]));

  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_random = (%.6e,%.3e,%.3e)\n",this_node,p->lc.f_random[0],p->lc.f_random[1],p->lc.f_random[2]));

  force[0] = force[0] + p->lc.f_random[0];
  force[1] = force[1] + p->lc.f_random[1];
  force[2] = force[2] + p->lc.f_random[2];

  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_tot = (%.6e,%.3e,%.3e)\n",this_node,force[0],force[1],force[2]));
      
  /* transform momentum transfer to lattice units
     (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  delta_j[0] = - force[0]*time_step*tau/agrid;
  delta_j[1] = - force[1]*time_step*tau/agrid;
  delta_j[2] = - force[2]*time_step*tau/agrid;
  
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

}

int lb_lbfluid_get_interpolated_velocity(double* p, double* v) {
  index_t node_index[8], index;
  double delta[6];
  double local_rho, local_j[3], interpolated_u[3];
  double modes[19];
  int x,y,z;
  double pos[3];

#ifdef LB_BOUNDARIES
  double lbboundary_mindist, distvec[3];
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
  map_position_to_lattice(&lblattice,pos,node_index,delta);

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  interpolated_u[0] = interpolated_u[1] = interpolated_u[2] = 0.0 ;

  for (z=0;z<2;z++) {
    for (y=0;y<2;y++) {
      for (x=0;x<2;x++) {
        	
        index = node_index[(z*2+y)*2+x];
        
#ifdef LB_BOUNDARIES
        if (lbfields[index].boundary) {
          local_rho=lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;
          local_j[0] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
          local_j[1] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
          local_j[2] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
        } else {
          lb_calc_modes(index, modes);
          local_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid + modes[0];
          local_j[0] = modes[1];
          local_j[1] = modes[2];
          local_j[2] = modes[3];
        }
#else 
        lb_calc_modes(index, modes);
        local_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid + modes[0];
        local_j[0] = modes[1];
        local_j[1] = modes[2];
        local_j[2] = modes[3];
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


/** Calculate particle lattice interactions.
 * So far, only viscous coupling with Stokesian friction is
 * implemented.
 * Include all particle-lattice forces in this function.
 * The function is called from \ref force_calc.
 *
 * Parallelizing the fluid particle coupling is not straightforward
 * because drawing of random numbers makes the whole thing nonlocal.
 * One way to do it is to treat every particle only on one node, i.e.
 * the random numbers need not be communicated. The particles that are 
 * not fully inside the local lattice are taken into account via their
 * ghost images on the neighbouring nodes. But this requires that the 
 * correct values of the surrounding lattice nodes are available on 
 * the respective node, which means that we have to communicate the 
 * halo regions before treating the ghost particles. Moreover, after 
 * determining the ghost couplings, we have to communicate back the 
 * halo region such that all local lattice nodes have the correct values.
 * Thus two communication phases are involved which will most likely be 
 * the bottleneck of the computation.
 *
 * Another way of dealing with the particle lattice coupling is to 
 * treat a particle and all of it's images explicitly. This requires the
 * communication of the random numbers used in the calculation of the 
 * coupling force. The problem is now that, if random numbers have to 
 * be redrawn, we cannot efficiently determine which particles and which 
 * images have to be re-calculated. We therefore go back to the outset
 * and go through the whole system again until no failure occurs during
 * such a sweep. In the worst case, this is very inefficient because
 * many things are recalculated although they actually don't need.
 * But we can assume that this happens extremely rarely and then we have
 * on average only one communication phase for the random numbers, which
 * probably makes this method preferable compared to the above one.
 */
void calc_particle_lattice_ia() {
  int i, c, np;
  Cell *cell ;
  Particle *p ;
  double force[3];


  if (transfer_momentum) {

    if (lbpar.resend_halo) { /* first MD step after last LB update */
      
      /* exchange halo regions (for fluid-particle coupling) */
      halo_communication(&update_halo_comm, (char*)**lbfluid);
#ifdef ADDITIONAL_CHECKS
      lb_check_halo_regions();
#endif
      
      /* halo is valid now */
      lbpar.resend_halo = 0;

      /* all fields have to be recalculated */
      for (i=0; i<lblattice.halo_grid_volume; ++i) {
        lbfields[i].recalc_fields = 1;
      }

    }
      
    /* draw random numbers for local particles */
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;
      for (i=0;i<np;i++) {
#ifdef GAUSSRANDOM
        p[i].lc.f_random[0] = lb_coupl_pref2*gaussian_random();
        p[i].lc.f_random[1] = lb_coupl_pref2*gaussian_random();
        p[i].lc.f_random[2] = lb_coupl_pref2*gaussian_random();
#elif defined (GAUSSRANDOMCUT)
        p[i].lc.f_random[0] = lb_coupl_pref2*gaussian_random_cut();
        p[i].lc.f_random[1] = lb_coupl_pref2*gaussian_random_cut();
        p[i].lc.f_random[2] = lb_coupl_pref2*gaussian_random_cut();
#elif defined (FLATNOISE)
        p[i].lc.f_random[0] = lb_coupl_pref*(d_random()-0.5);
        p[i].lc.f_random[1] = lb_coupl_pref*(d_random()-0.5);
        p[i].lc.f_random[2] = lb_coupl_pref*(d_random()-0.5);
#else
#error No noise type defined for the CPU LB
#endif

#ifdef ADDITIONAL_CHECKS
        rancounter += 3;
#endif
      }
    }
    
    /* communicate the random numbers */
    ghost_communicator(&cell_structure.ghost_lbcoupling_comm) ;
    
    /* local cells */
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {

        lb_viscous_coupling(&p[i],force);

        /* add force to the particle */
        p[i].f.f[0] += force[0];
        p[i].f.f[1] += force[1];
        p[i].f.f[2] += force[2];

        ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.6e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
  
      }

    }

    /* ghost cells */
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

          ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: LB coupling of ghost particle:\n",this_node));

          lb_viscous_coupling(&p[i],force);

          /* ghosts must not have the force added! */

          ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.6e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));

        }
      }
    }

  }
}

/***********************************************************************/

/** Calculate the average density of the fluid in the system.
 * This function has to be called after changing the density of
 * a local lattice site in order to set lbpar.rho consistently. */
void lb_calc_average_rho() {

  index_t index;
  int x, y, z;
  double rho, local_rho, sum_rho;

  rho = 0.0;
  local_rho = 0.0;
  index = 0;
  for (z=1; z<=lblattice.grid[2]; z++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (x=1; x<=lblattice.grid[0]; x++) {

        lb_calc_local_rho(index, &rho);
        local_rho += rho;

        index++;
      }
      index += 2;
    }
    index += 2*lblattice.halo_grid[0];
  }

  MPI_Allreduce(&rho, &sum_rho, 1, MPI_DOUBLE, MPI_SUM, comm_cart);

  /* calculate average density in MD units */
  // TODO!!!
  lbpar.rho[0] = sum_rho / (box_l[0]*box_l[1]*box_l[2]);
  
}

/*@}*/

#ifdef ADDITIONAL_CHECKS
static int compare_buffers(double *buf1, double *buf2, int size) {
  int ret;
  if (memcmp(buf1,buf2,size)) {
    char *errtxt;
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{102 Halo buffers are not identical} ");
    ret = 1;
  } else {
    ret = 0;
  }
  return ret;
}

/** Checks consistency of the halo regions (ADDITIONAL_CHECKS)
    This function can be used as an additional check. It test whether the 
    halo regions have been exchanged correctly.
*/
static void lb_check_halo_regions()
{

  index_t index;
  int i,x,y,z, s_node, r_node, count=n_veloc;
  double *s_buffer, *r_buffer;
  MPI_Status status[2];

  r_buffer = (double*) malloc(count*sizeof(double));
  s_buffer = (double*) malloc(count*sizeof(double));

  if (PERIODIC(0)) {
    for (z=0;z<lblattice.halo_grid[2];++z) {
      for (y=0;y<lblattice.halo_grid[1];++y) {

	index  = get_linear_index(0,y,z,lblattice.halo_grid);
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[1];
	r_node = node_neighbors[0];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(lblattice.grid[0],y,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(lblattice.grid[0],y,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  if (compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld y=%d z=%d\n",0,index,y,z);
	  }
	}

	index = get_linear_index(lblattice.grid[0]+1,y,z,lblattice.halo_grid); 
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[0];
	r_node = node_neighbors[1];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(1,y,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(1,y,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) r_buffer[i] = lbfluid[0][i][index];
	  if (compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld y=%d z=%d\n",0,index,y,z);	  
	  }
	}
	
      }      
    }
  }

  if (PERIODIC(1)) {
    for (z=0;z<lblattice.halo_grid[2];++z) {
      for (x=0;x<lblattice.halo_grid[0];++x) {

	index = get_linear_index(x,0,z,lblattice.halo_grid);
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[3];
	r_node = node_neighbors[2];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(x,lblattice.grid[1],z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(x,lblattice.grid[1],z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) r_buffer[i] = lbfluid[0][i][index];
	  if (compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld x=%d z=%d\n",1,index,x,z);
	  }
	}

      }
      for (x=0;x<lblattice.halo_grid[0];++x) {

	index = get_linear_index(x,lblattice.grid[1]+1,z,lblattice.halo_grid);
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[2];
	r_node = node_neighbors[3];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(x,1,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(x,1,z,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) r_buffer[i] = lbfluid[0][i][index];
	  if (compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld x=%d z=%d\n",1,index,x,z);
	  }
	}

      }
    }
  }

  if (PERIODIC(2)) {
    for (y=0;y<lblattice.halo_grid[1];++y) {
      for (x=0;x<lblattice.halo_grid[0];++x) {

	index = get_linear_index(x,y,0,lblattice.halo_grid);
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[5];
	r_node = node_neighbors[4];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(x,y,lblattice.grid[2],lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(x,y,lblattice.grid[2],lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) r_buffer[i] = lbfluid[0][i][index];
	  if (compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld x=%d y=%d z=%d\n",2,index,x,y,lblattice.grid[2]);  
	  }
	}

      }
    }
    for (y=0;y<lblattice.halo_grid[1];++y) {
      for (x=0;x<lblattice.halo_grid[0];++x) {

	index = get_linear_index(x,y,lblattice.grid[2]+1,lblattice.halo_grid);
	for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];

	s_node = node_neighbors[4];
	r_node = node_neighbors[5];
	if (n_nodes > 1) {
	  MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		       comm_cart, status);
	  index = get_linear_index(x,y,1,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) s_buffer[i] = lbfluid[0][i][index];
	  compare_buffers(s_buffer,r_buffer,count*sizeof(double));
	} else {
	  index = get_linear_index(x,y,1,lblattice.halo_grid);
	  for (i=0;i<n_veloc;i++) r_buffer[i] = lbfluid[0][i][index];
	  if(compare_buffers(s_buffer,r_buffer,count*sizeof(double))) {
	    fprintf(stderr,"buffers differ in dir=%d at index=%ld x=%d y=%d\n",2,index,x,y);
	  }
	}
      
      }
    }
  }

  free(r_buffer);
  free(s_buffer);

  //if (check_runtime_errors());
  //else fprintf(stderr,"halo check successful\n");

}
#endif /* ADDITIONAL_CHECKS */

#if 0 /* These debug functions are used nowhere. If you need it, here they are.
        Remove this comment line and the matching #endif.
        The functions in question are:
            lb_lattice_sum
            lb_check_mode_transformation
            lb_init_mode_transformation
            lb_check_negative_n
        */
#ifdef ADDITIONAL_CHECKS
static void lb_lattice_sum() {

    int n_veloc = lbmodel.n_veloc;
    double *w   = lbmodel.w;
    double (*v)[3]  = lbmodel.c;
    
    //int n_veloc = 14;
    //double w[14]    = { 7./18., 
    //                    1./12., 1./12., 1./12., 1./12., 1./18.,
    //                    1./36., 1./36., 1./36., 1./36., 
    //                    1./36., 1./36., 1./36., 1./36. };
    //double v[14][3] = { { 0., 0., 0. },
    //                    { 1., 0., 0. },
    //                    {-1., 0., 0. },
    //                    { 0., 1., 0. },
    //		        { 0.,-1., 0. },
    //                    { 0., 0., 1. },
    //                    { 1., 1., 0. },
    //                    {-1.,-1., 0. },
    //                    { 1.,-1., 0. },
    //                    {-1., 1., 0. },
    //                    { 1., 0., 1. },
    //                    {-1., 0., 1. },
    //                    { 0., 1., 1. },
    //                    { 0.,-1., 1. } };

    int i,a,b,c,d,e;
    double sum1,sum2,sum3,sum4,sum5;
    int count=0;

    for (a=0; a<3; a++) 
      {
	sum1 = 0.0;
	for (i=0; i<n_veloc; ++i) {
	  if (v[i][2] < 0) sum1 += w[i]*v[i][a];      
	}
	if (fabs(sum1) > ROUND_ERROR_PREC) {
	  count++; fprintf(stderr,"(%d) %f\n",a,sum1);
	}
      }

    for (a=0; a<3; a++)
      for (b=0; b<3; b++) 
	{
	  sum2 = 0.0;
	  for (i=0; i<n_veloc; ++i) {
	    if (v[i][2] < 0) sum2 += w[i]*v[i][a]*v[i][b];      
	  }
	  if (sum2!=0.0) {
	    count++; fprintf(stderr,"(%d,%d) %f\n",a,b,sum2);
	  }
	}

    for (a=0; a<3; a++)
      for (b=0; b<3; b++) 
	for (c=0; c<3; c++) 
	  {
	    sum3 = 0.0;
	    for (i=0; i<n_veloc; ++i) {
	      if (v[i][2] < 0) sum3 += w[i]*v[i][a]*v[i][b]*v[i][c];      
	    }
	    if (sum3!=0.0) {
	      count++; fprintf(stderr,"(%d,%d,%d) %f\n",a,b,c,sum3);
	    }
	  }

    for (a=0; a<3; a++)
      for (b=0; b<3; b++)
	for (c=0; c<3; c++)
	  for (d=0; d<3; d++)
	    {
	      sum4 = 0.0;
	      for (i=0; i<n_veloc; ++i) {
		if (v[i][2] < 0) sum4 += w[i]*v[i][a]*v[i][b]*v[i][c]*v[i][d];      
	      }
	      if (fabs(sum4) > ROUND_ERROR_PREC) { 
		  count++; fprintf(stderr,"(%d,%d,%d,%d) %f\n",a,b,c,d,sum4); 
	      }
	    }

    for (a=0; a<3; a++)
      for (b=0; b<3; b++)
	for (c=0; c<3; c++)
	  for (d=0; d<3; d++)
	    for (e=0; e<3; e++) 
	      {
		sum5 = 0.0;
		for (i=0; i<n_veloc; ++i) {
		  if (v[i][2] < 0) sum5 += w[i]*v[i][a]*v[i][b]*v[i][c]*v[i][d]*v[i][e];      
		}
		if (fabs(sum5) > ROUND_ERROR_PREC) { 
		  count++; fprintf(stderr,"(%d,%d,%d,%d,%d) %f\n",a,b,c,d,e,sum5);
		}
	      }

    fprintf(stderr,"%d non-null entries\n",count);

}
#endif /* #ifdef ADDITIONAL_CHECKS */

#ifdef ADDITIONAL_CHECKS
static void lb_check_mode_transformation(index_t index, double *mode) {

  /* check if what I think is right */

  int i;
  double *w = lbmodel.w;
  double (*e)[19] = d3q19_modebase;
  double sum_n=0.0, sum_m=0.0;
  double n_eq[19];
  double m_eq[19];
  // unit conversion: mass density
  double avg_rho = lbpar.rho*lbpar.agrid*lbpar.agrid*lbpar.agrid;
  double (*c)[3] = lbmodel.c;

  m_eq[0] = mode[0];
  m_eq[1] = mode[1];
  m_eq[2] = mode[2];
  m_eq[3] = mode[3];

  double rho = mode[0] + avg_rho;
  double *j  = mode+1;

  /* equilibrium part of the stress modes */
  /* remember that the modes have (\todo not?) been normalized! */
  m_eq[4] = /*1./6.*/scalar(j,j)/rho;
  m_eq[5] = /*1./4.*/(SQR(j[0])-SQR(j[1]))/rho;
  m_eq[6] = /*1./12.*/(scalar(j,j) - 3.0*SQR(j[2]))/rho;
  m_eq[7] = j[0]*j[1]/rho;
  m_eq[8] = j[0]*j[2]/rho;
  m_eq[9] = j[1]*j[2]/rho;

  for (i=10;i<n_veloc;i++) {
    m_eq[i] = 0.0;
  }

  for (i=0;i<n_veloc;i++) {
    n_eq[i] = w[i]*((rho-avg_rho) + 3.*scalar(j,c[i]) + 9./2.*SQR(scalar(j,c[i]))/rho - 3./2.*scalar(j,j)/rho);
  } 

  for (i=0;i<n_veloc;i++) {
    sum_n += SQR(lbfluid[0][i][index]-n_eq[i])/w[i];
    sum_m += SQR(mode[i]-m_eq[i])/e[19][i];
  }

  if (fabs(sum_n-sum_m)>ROUND_ERROR_PREC) {    
    fprintf(stderr,"Attention: sum_n=%f sum_m=%f %e\n",sum_n,sum_m,fabs(sum_n-sum_m));
  }

}

static void lb_init_mode_transformation() {

#ifdef D3Q19
  int i, j, k, l;
  int n_veloc = 14;
  double w[14]    = { 7./18., 
                      1./12., 1./12., 1./12., 1./12., 1./18.,
                      1./36., 1./36., 1./36., 1./36., 
                      1./36., 1./36., 1./36., 1./36. };
  double c[14][3] = { { 0., 0., 0. },
                      { 1., 0., 0. },
                      {-1., 0., 0. },
                      { 0., 1., 0. },
            		      { 0.,-1., 0. },
                      { 0., 0., 1. },
                      { 1., 1., 0. },
                      {-1.,-1., 0. },
                      { 1.,-1., 0. },
                      {-1., 1., 0. },
                      { 1., 0., 1. },
                      {-1., 0., 1. },
                      { 0., 1., 1. },
                      { 0.,-1., 1. } };

  double b[19][14];
  double e[14][14];
  double proj, norm[14];

  /* construct polynomials from the discrete velocity vectors */
  for (i=0;i<n_veloc;i++) {
    b[0][i]  = 1;
    b[1][i]  = c[i][0];
    b[2][i]  = c[i][1];
    b[3][i]  = c[i][2];
    b[4][i]  = scalar(c[i],c[i]);
    b[5][i]  = c[i][0]*c[i][0]-c[i][1]*c[i][1];
    b[6][i]  = scalar(c[i],c[i])-3*c[i][2]*c[i][2];
    //b[5][i]  = 3*c[i][0]*c[i][0]-scalar(c[i],c[i]);
    //b[6][i]  = c[i][1]*c[i][1]-c[i][2]*c[i][2];
    b[7][i]  = c[i][0]*c[i][1];
    b[8][i]  = c[i][0]*c[i][2];
    b[9][i]  = c[i][1]*c[i][2];
    b[10][i] = 3*scalar(c[i],c[i])*c[i][0];
    b[11][i] = 3*scalar(c[i],c[i])*c[i][1];
    b[12][i] = 3*scalar(c[i],c[i])*c[i][2];
    b[13][i] = (c[i][1]*c[i][1]-c[i][2]*c[i][2])*c[i][0];
    b[14][i] = (c[i][0]*c[i][0]-c[i][2]*c[i][2])*c[i][1];
    b[15][i] = (c[i][0]*c[i][0]-c[i][1]*c[i][1])*c[i][2];
    b[16][i] = 3*scalar(c[i],c[i])*scalar(c[i],c[i]);
    b[17][i] = 2*scalar(c[i],c[i])*b[5][i];
    b[18][i] = 2*scalar(c[i],c[i])*b[6][i];
  }

  for (i=0;i<n_veloc;i++) {
    b[0][i]  = 1;
    b[1][i]  = c[i][0];
    b[2][i]  = c[i][1];
    b[3][i]  = c[i][2];
    b[4][i]  = scalar(c[i],c[i]);
    b[5][i]  = SQR(c[i][0])-SQR(c[i][1]);
    b[6][i]  = c[i][0]*c[i][1];
    b[7][i]  = c[i][0]*c[i][2];
    b[8][i]  = c[i][1]*c[i][2];
    b[9][i]  = scalar(c[i],c[i])*c[i][0];
    b[10][i] = scalar(c[i],c[i])*c[i][1];
    b[11][i] = scalar(c[i],c[i])*c[i][2];
    b[12][i] = (c[i][0]*c[i][0]-c[i][1]*c[i][1])*c[i][2];
    b[13][i] = scalar(c[i],c[i])*scalar(c[i],c[i]);
  }

  /* Gram-Schmidt orthogonalization procedure */
  for (j=0;j<n_veloc;j++) {
    for (i=0;i<n_veloc;i++) e[j][i] = b[j][i];
    for (k=0;k<j;k++) {
      proj = 0.0;
      for (l=0;l<n_veloc;l++) {
	proj += w[l]*e[k][l]*b[j][l];
      }
      if (j==13) fprintf(stderr,"%d %f\n",k,proj/norm[k]);
      for (i=0;i<n_veloc;i++) e[j][i] -= proj/norm[k]*e[k][i];
    }
    norm[j] = 0.0;
    for (i=0;i<n_veloc;i++) norm[j] += w[i]*SQR(e[j][i]);
  }
  
  fprintf(stderr,"e[%d][%d] = {\n",n_veloc,n_veloc);
  for (i=0;i<n_veloc;i++) {
    fprintf(stderr,"{ % .3f",e[i][0]);
    for (j=1;j<n_veloc;j++) {
      fprintf(stderr,", % .3f",e[i][j]);
    }
    fprintf(stderr," } %.9f\n",norm[i]);
  }
  fprintf(stderr,"};\n");

  /* projections on lattice tensors */
  for (i=0;i<n_veloc;i++) {
    proj = 0.0;
    for (k=0;k<n_veloc;k++) {
      proj += e[i][k] * w[k] * 1;
    }
    fprintf(stderr, "%.6f",proj);
    
    for (j=0;j<3;j++) {
      proj = 0.0;
      for (k=0;k<n_veloc;k++) {
	proj += e[i][k] * w[k] * c[k][j];
      }
      fprintf(stderr, " %.6f",proj);
    }

    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	proj=0.0;
	for (l=0;l<n_veloc;l++) {
	  proj += e[i][l] * w[l] * c[l][j] * c[l][k];
	}
	fprintf(stderr, " %.6f",proj);
      }
    }

    fprintf(stderr,"\n");

  }

  //proj = 0.0;
  //for (k=0;k<n_veloc;k++) {
  //  proj += c[k][2] * w[k] * 1;
  //}
  //fprintf(stderr,"%.6f",proj);
  //
  //proj = 0.0;
  //for (k=0;k<n_veloc;k++) {
  //  proj += c[k][2] * w[k] * c[k][2];
  //}
  //fprintf(stderr," %.6f",proj);
  //
  //proj = 0.0;
  //for (k=0;k<n_veloc;k++) {
  //  proj += c[k][2] * w[k] * c[k][2] * c[k][2];
  //}
  //fprintf(stderr," %.6f",proj);
  //
  //fprintf(stderr,"\n");

#else /* not D3Q19 */
  int i, j, k, l;
  double b[9][9];
  double e[9][9];
  double proj, norm[9];

  double c[9][2] = { { 0, 0 },
		     { 1, 0 },
		     {-1, 0 },
                     { 0, 1 },
                     { 0,-1 },
		     { 1, 1 },
		     {-1,-1 },
		     { 1,-1 },
		     {-1, 1 } };

  double w[9] = { 4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36 };

  n_veloc = 9;

  /* construct polynomials from the discrete velocity vectors */
  for (i=0;i<n_veloc;i++) {
    b[0][i] = 1;
    b[1][i] = c[i][0];
    b[2][i] = c[i][1];
    b[3][i] = 3*(SQR(c[i][0]) + SQR(c[i][1]));
    b[4][i] = c[i][0]*c[i][0]-c[i][1]*c[i][1];
    b[5][i] = c[i][0]*c[i][1];
    b[6][i] = 3*(SQR(c[i][0])+SQR(c[i][1]))*c[i][0];
    b[7][i] = 3*(SQR(c[i][0])+SQR(c[i][1]))*c[i][1];
    b[8][i] = (b[3][i]-5)*b[3][i]/2;
  }

  /* Gram-Schmidt orthogonalization procedure */
  for (j=0;j<n_veloc;j++) {
    for (i=0;i<n_veloc;i++) e[j][i] = b[j][i];
    for (k=0;k<j;k++) {
      proj = 0.0;
      for (l=0;l<n_veloc;l++) {
	proj += w[l]*e[k][l]*b[j][l];
      }
      for (i=0;i<n_veloc;i++) e[j][i] -= proj/norm[k]*e[k][i];
    }
    norm[j] = 0.0;
    for (i=0;i<n_veloc;i++) norm[j] += w[i]*SQR(e[j][i]);
  }
  
  fprintf(stderr,"e[%d][%d] = {\n",n_veloc,n_veloc);
  for (i=0;i<n_veloc;i++) {
    fprintf(stderr,"{ % .1f",e[i][0]);
    for (j=1;j<n_veloc;j++) {
      fprintf(stderr,", % .1f",e[i][j]);
    }
    fprintf(stderr," } %.2f\n",norm[i]);
  }
  fprintf(stderr,"};\n");

#endif

}
#endif /* ADDITIONAL_CHECKS */

#ifdef ADDITIONAL_CHECKS

/** Check for negative populations.  
 
    Checks for negative populations and increases failcounter for each
    occurence.
 
    @param  index Index of the local lattice site (Input).
    @return Number of negative populations on the local lattice site. */
static int lb_check_negative_n(index_t index)
{
  int i, localfails=0;

  for (i=0; i<n_veloc; i++) {
    if (lbfluid[0][i][index]+lbmodel.coeff[i][0]*lbpar.rho < 0.0) {
      ++localfails;
      ++failcounter;
      fprintf(stderr,"%d: Negative population n[%d]=%le (failcounter=%d, rancounter=%d).\n   Check your parameters if this occurs too often!\n",this_node,i,lbmodel.coeff[i][0]*lbpar.rho+lbfluid[0][i][index],failcounter,rancounter);
      break;
   }
  }

  return localfails;
}
#endif /* ADDITIONAL_CHECKS */
#endif /* #if 0 */
/* Here, the unused "ADDITIONAL_CHECKS functions end. */

#endif
