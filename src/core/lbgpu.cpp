/* 
   Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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

/** \file lbgpu.cpp
 *
 * C file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.hpp.
 */

#include "config.hpp"

#ifdef LB_GPU

#include "lbgpu.hpp"
#include "utils.hpp"
#include "communication.hpp"
#include "thermostat.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "global.hpp"
#include "lbboundaries.hpp"
#include "cuda_interface.hpp"
#include "statistics.hpp"
#include "partCfg_global.hpp"
#include "lb.hpp"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <random>

#ifndef D3Q19
#error The implementation only works for D3Q19 so far!
#endif

#if (LB_COMPONENTS==1) 
#    define  SC0 {0.0}
#    define  SC20 {0.0, 0.0}
#    define  SC1 {1.0}
#    define  SCM1 {-1.0}
#endif 
#if (LB_COMPONENTS==2) 
#    define  SC0 { 0.0 , 0.0 } 
#    define  SC20 {0.0, 0.0, 0.0, 0.0}
#    define  SC1 { 1.0 , 1.0 } 
#    define  SCM1 { -1.0, -1.0 } 
#endif 


/** Struct holding the Lattice Boltzmann parameters */
// LB_parameters_gpu lbpar_gpu = { .rho=SC0, .mu=SC0, .viscosity=SC0, .gamma_shear=SC0, .gamma_bulk=SC0,
//                                 .gamma_odd=SC0,.gamma_even=SC0, .agrid=0.0, .tau=-1.0, .friction=SC0, .time_step=0.0, .lb_coupl_pref=SC1 ,
//                                 .lb_coupl_pref2=SC0, .bulk_viscosity=SCM1, .dim_x=0, .dim_y=0, .dim_z=0, .number_of_nodes=0, 
//                                 .number_of_particles=0, .fluct=0, .calc_val=1, .external_force=0, .ext_force={0.0, 0.0, 0.0}, 
//                                  .your_seed=12345, .reinit=0,
// #ifdef SHANCHEN
//                                 .coupling=SC20, .gamma_mobility=SC1
// #endif
// };
LB_parameters_gpu lbpar_gpu = {
  // rho
  SC0, 
  // mu
  SC0,
  // viscosity
  SC0,
  // gamma_shear
  SC0,
  // gamma_bulk
  SC0, 
  // gamma_odd
  SC0,
  // gamma_even
  SC0,
  //is_TRT
  false,
  // friction
  SC0,
  // lb_couple_switch
  LB_COUPLE_TWO_POINT,
  // lb_coupl_pref
  SC0,
  // lb_coupl_pref2
  SC0,
  // bulk_viscosity
  SCM1,
  // agrid
  0.0,
  // tau
  -1.0,
  // time_step
  0.0,
  // dim_x;
  0,
  // dim_y;
  0,
  // dim_z;
  0,
  // number_of_nodes
  0,
  // number_of_particles
  0,
#ifdef LB_BOUNDARIES_GPU
  // number_of_boundnodes
  0,
#endif
  // fluct
  0,
  // calc_val
  1,
  // external_force
  0,
  // ext_force
  {0.0, 0.0, 0.0},
  // your_seed
  12345,
  // reinit
  0,
#ifdef SHANCHEN
  // gamma_mobility
  SC1,
  // mobility
  SC1,
  // coupling
  SC20,
  // remove_momentum
  0
#endif // SHANCHEN  
};

/** this is the array that stores the hydrodynamic fields for the output */
LB_rho_v_pi_gpu *host_values = nullptr;

LB_nodes_gpu *host_nodes = nullptr;

/** Flag indicating momentum exchange between particles and fluid */
int transfer_momentum_gpu = 0;

static int max_ran = 1000000;
/*@}*/
//static double tau;

/** measures the MD time since the last fluid update */
static int fluidstep = 0;

/** c_sound_square in LB units*/
static float c_sound_sq = 1.0f/3.0f;

//clock_t start, end;
int i;


int n_extern_node_force_densities = 0;
LB_extern_nodeforcedensity_gpu *host_extern_node_force_densities = nullptr;
int ek_initialized = 0;

/*-----------------------------------------------------------*/
/** main of lb_gpu_programm */
/*-----------------------------------------------------------*/
#ifdef SHANCHEN
/* called from forces.cpp. This is at the beginning of the force
   calculation loop, so we increment the fluidstep counter here,
   and we reset it only when the last call to a LB function
   [lattice_boltzmann_update_gpu()] is performed within integrate_vv()
 */
void lattice_boltzmann_calc_shanchen_gpu(void){

  int factor = (int)round(lbpar_gpu.tau/time_step);

  if (fluidstep+1 >= factor) 
     lb_calc_shanchen_GPU();
}
#endif //SHANCHEN

/** lattice boltzmann update gpu called from integrate.cpp
*/

void lattice_boltzmann_update_gpu() {

  int factor = (int)round(lbpar_gpu.tau/time_step);

  fluidstep += 1;

  if (fluidstep>=factor) {

    fluidstep=0; 
    lb_integrate_GPU();
#ifdef SHANCHEN
    if(lbpar_gpu.remove_momentum) lb_remove_fluid_momentum_GPU();
#endif
    LB_TRACE (fprintf(stderr,"lb_integrate_GPU \n"));

  }
}

/** (re-) allocation of the memory needed for the particles (cpu part)
*/
void lb_realloc_particles_gpu(){

  lbpar_gpu.number_of_particles = n_part;
  LB_TRACE (printf("#particles realloc\t %u \n", lbpar_gpu.number_of_particles));
  lbpar_gpu.your_seed = (unsigned int)std::random_device{}();

  LB_TRACE (fprintf(stderr,"test your_seed %u \n", lbpar_gpu.your_seed));

  lb_realloc_particles_GPU_leftovers(&lbpar_gpu);
}

/** (Re-)initializes the fluid according to the given value of rho. */
void lb_reinit_fluid_gpu() {

  lb_reinit_parameters_gpu();
  if(lbpar_gpu.number_of_nodes != 0){
    lb_reinit_GPU(&lbpar_gpu);
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    lbpar_gpu.reinit = 1;
  }

  LB_TRACE (fprintf(stderr,"lb_reinit_fluid_gpu \n"));
}

/** Release the fluid. */
/*not needed in Espresso but still not deleted.
  Despite the name (TODO: change it), it releases 
  only the fluid-related memory on the gpu.*/
void lb_release_gpu(){

  if(host_nodes !=nullptr) { free(host_nodes); host_nodes=nullptr ;} 
  if(host_values!=nullptr) { free(host_values); host_values=nullptr;}
//  if(host_forces!=nullptr) free(host_forces);
//  if(host_data  !=nullptr) free(host_data);
}
/** (Re-)initializes the fluid. */
void lb_reinit_parameters_gpu() {
  int ii;

  lbpar_gpu.time_step = (float)time_step;
  for(ii=0;ii<LB_COMPONENTS;++ii){
    lbpar_gpu.mu[ii] = 0.0;
  
    if (lbpar_gpu.viscosity[ii] > 0.0) {
      /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
      lbpar_gpu.gamma_shear[ii] = 1. - 2./(6.*lbpar_gpu.viscosity[ii]*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid) + 1.);   
    }
  
    if (lbpar_gpu.bulk_viscosity[ii] > 0.0) {
      /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
      lbpar_gpu.gamma_bulk[ii] = 1. - 2./(9.*lbpar_gpu.bulk_viscosity[ii]*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid) + 1.);
    }

    //By default, gamma_even and gamma_odd are chosen such that the MRT becomes
    //a TRT with ghost mode relaxation factors that minimize unphysical wall
    //slip at bounce-back boundaries. For the relation between the gammas
    //achieving this, consult
    //  D. d’Humières, I. Ginzburg, Comp. & Math. w. App. 58(5):823–840 (2009)
    //Note that the relaxation operator in Espresso is defined as
    //  m* = m_eq + gamma * (m - m_eq)
    //as opposed to this reference, where
    //  m* = m + lambda * (m - m_eq)
    
    if (lbpar_gpu.is_TRT) {
      lbpar_gpu.gamma_bulk[ii] = lbpar_gpu.gamma_shear[ii];
      lbpar_gpu.gamma_even[ii] = lbpar_gpu.gamma_shear[ii];
      lbpar_gpu.gamma_odd[ii] = -(7.0f*lbpar_gpu.gamma_even[ii]+1.0f)/(lbpar_gpu.gamma_even[ii]+7.0f);
      //lbpar_gpu.gamma_odd[ii] = lbpar_gpu.gamma_shear[ii]; //uncomment for BGK
    }
    
    //lbpar_gpu.gamma_even[ii] = 0.0; //uncomment for special case of BGK
    //lbpar_gpu.gamma_odd[ii] = 0.0;
    //lbpar_gpu.gamma_shear[ii] = 0.0;
    //lbpar_gpu.gamma_bulk[ii] = 0.0;
    
    //printf("gamma_shear=%e\n", lbpar_gpu.gamma_shear[ii]);
    //printf("gamma_bulk=%e\n", lbpar_gpu.gamma_bulk[ii]);
    //printf("TRT gamma_odd=%e\n", lbpar_gpu.gamma_odd[ii]);
    //printf("TRT gamma_even=%e\n", lbpar_gpu.gamma_even[ii]);
    //printf("\n");

#ifdef SHANCHEN
    if (lbpar_gpu.mobility[0] > 0.0) {
      lbpar_gpu.gamma_mobility[0] = 1. - 2./(6.*lbpar_gpu.mobility[0]*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid) + 1.);
    }
#endif
    if (temperature > 0.0) {  /* fluctuating hydrodynamics ? */
  
      lbpar_gpu.fluct = 1;
  	LB_TRACE (fprintf(stderr, "fluct on \n"));
      /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).*/
      /* Note that the modes are not normalized as in the paper here! */
      lbpar_gpu.mu[ii] = (float)temperature*lbpar_gpu.tau*lbpar_gpu.tau/c_sound_sq/(lbpar_gpu.agrid*lbpar_gpu.agrid); 
  
      /* lb_coupl_pref is stored in MD units (force)
       * Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
       * The factor 12 comes from the fact that we use random numbers
       * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
       * time_step comes from the discretization.
       */
  
      lbpar_gpu.lb_coupl_pref[ii] = sqrtf(12.f*2.f*lbpar_gpu.friction[ii]*(float)temperature/lbpar_gpu.time_step);
      lbpar_gpu.lb_coupl_pref2[ii] = sqrtf(2.f*lbpar_gpu.friction[ii]*(float)temperature/lbpar_gpu.time_step);
  
    } else {
      /* no fluctuations at zero temperature */
      lbpar_gpu.fluct = 0;
      lbpar_gpu.lb_coupl_pref[ii] = 0.0;
      lbpar_gpu.lb_coupl_pref2[ii] = 0.0;
    }
  	LB_TRACE (fprintf(stderr,"lb_reinit_prarameters_gpu \n"));
  }


#ifdef ELECTROKINETICS
  if (ek_initialized) {
    lbpar_gpu.dim_x = (unsigned int) round(box_l[0] / lbpar_gpu.agrid); //TODO code duplication with lb.c start
    lbpar_gpu.dim_y = (unsigned int) round(box_l[1] / lbpar_gpu.agrid);
    lbpar_gpu.dim_z = (unsigned int) round(box_l[2] / lbpar_gpu.agrid);
    
    unsigned int tmp[3];
    
    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;
    
    /* sanity checks */
    int dir;
    
    for (dir=0;dir<3;dir++) {
    /* check if box_l is compatible with lattice spacing */
      if (fabs(box_l[dir] - tmp[dir] * lbpar_gpu.agrid) > 1.0e-3) {
          runtimeErrorMsg() <<"Lattice spacing lbpar_gpu.agrid= "<< lbpar_gpu.agrid << " is incompatible with box_l[" << dir << "]=" << box_l[dir];
      }
    }
    
    lbpar_gpu.number_of_nodes = lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
    lbpar_gpu.tau = (float) time_step; //TODO code duplication with lb.c end
  }
#endif
  
	LB_TRACE (fprintf(stderr,"lb_reinit_prarameters_gpu \n"));

  reinit_parameters_GPU(&lbpar_gpu);
}

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init_gpu() {

  LB_TRACE(printf("Begin initialzing fluid on GPU\n"));
  /** set parameters for transfer to gpu */
  lb_reinit_parameters_gpu();

  lb_realloc_particles_gpu();

  lb_init_GPU(&lbpar_gpu);
  
  gpu_init_particle_comm();
  cuda_bcast_global_part_params();

  LB_TRACE(printf("Initialzing fluid on GPU successful\n"));
}

/*@}*/

int lb_lbnode_set_extforce_density_GPU(int ind[3], double f[3])
{
  if ( ind[0] < 0 || ind[0] >= int(lbpar_gpu.dim_x) ||
       ind[1] < 0 || ind[1] >= int(lbpar_gpu.dim_y) ||
       ind[2] < 0 || ind[2] >= int(lbpar_gpu.dim_z) )
    return ES_ERROR;

  unsigned int index =
    ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;

  size_t  size_of_extforces = (n_extern_node_force_densities+1)*sizeof(LB_extern_nodeforcedensity_gpu);
  host_extern_node_force_densities = (LB_extern_nodeforcedensity_gpu*) Utils::realloc(host_extern_node_force_densities, size_of_extforces);
  
  host_extern_node_force_densities[n_extern_node_force_densities].force_density[0] = (float)f[0];
  host_extern_node_force_densities[n_extern_node_force_densities].force_density[1] = (float)f[1];
  host_extern_node_force_densities[n_extern_node_force_densities].force_density[2] = (float)f[2];
  
  host_extern_node_force_densities[n_extern_node_force_densities].index = index;
  n_extern_node_force_densities++;
  
  if(lbpar_gpu.external_force_density == 0)lbpar_gpu.external_force_density = 1;

  lb_init_extern_nodeforcedensities_GPU(n_extern_node_force_densities, host_extern_node_force_densities, &lbpar_gpu);

  return ES_OK;
}

void lb_GPU_sanity_checks()
{
  if(this_node == 0){
    if (lbpar_gpu.agrid < 0.0) {
        runtimeErrorMsg() <<"Lattice Boltzmann agrid not set";
    }
    if (lbpar_gpu.tau < 0.0) {
        runtimeErrorMsg() <<"Lattice Boltzmann time step not set";
    }
    for(int i=0;i<LB_COMPONENTS;i++){
      if (lbpar_gpu.rho[0] < 0.0) {
          runtimeErrorMsg() <<"Lattice Boltzmann fluid density not set";
      }
      if (lbpar_gpu.viscosity[0] < 0.0) {
          runtimeErrorMsg() <<"Lattice Boltzmann fluid viscosity not set";
      }
    }
  }
}

int lb_lbfluid_save_checkpoint_wrapper(char* filename, int binary)
{
  return lb_lbfluid_save_checkpoint(filename, binary);
}

int lb_lbfluid_load_checkpoint_wrapper(char* filename, int binary)
{
  return lb_lbfluid_load_checkpoint(filename, binary);
}

void lb_lbfluid_particles_add_momentum(float momentum[3]) {
  auto & parts = partCfg();
  auto const n_part = parts.size();

  for (auto const &p : parts) {
    double new_velocity[3] = {
        p.m.v[0] +
            momentum[0] / p.p.mass / n_part,
        p.m.v[1] +
            momentum[1] / p.p.mass / n_part,
        p.m.v[2] +
            momentum[2] / p.p.mass / n_part};
    set_particle_v(p.p.identity, new_velocity);
  }
}

void lb_lbfluid_calc_linear_momentum ( float momentum[3], int include_particles, int include_lbfluid )
{
  auto linear_momentum = calc_linear_momentum(include_particles, include_lbfluid);
  momentum[0] = linear_momentum[0];
  momentum[1] = linear_momentum[1];
  momentum[2] = linear_momentum[2];
}

#endif /* LB_GPU */
