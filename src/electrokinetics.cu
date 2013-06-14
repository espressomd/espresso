/* 
   Copyright (C) 2010,2011,2012 The ESPResSo project

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

#include <cuda.h>
#include <cufft.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
//#include "lb-boundaries.h" //TODO: needed to get rid of the code duplication below
#include "config.h"
#include "electrokinetics.h"
#include "cuda_common.h"
#include "lbgpu.h"
#include "constraint.h"

#ifdef __cplusplus
}
#endif


#ifdef ELECTROKINETICS

#ifdef __cplusplus
extern "C" {
#endif

  /* TODO: get rid of this code duplication with lb-boundaries.h by solving the
           cuda-mpi incompatibility */
           
#define LATTICE_LB_GPU   2
extern int lattice_switch;
           
#ifdef EK_BOUNDARIES

  typedef struct {

    int type;
    double slip_pref;

    union {
      Constraint_wall wal;
      Constraint_sphere sph;
      Constraint_cylinder cyl;
      Constraint_rhomboid rhomboid;
      Constraint_pore pore;
    } c;
    
    double force[3];
    double velocity[3];
    
  #ifdef EK_BOUNDARIES
    float charge_density;
  #endif

  } LB_Boundary;

  extern int n_lb_boundaries;
  extern LB_Boundary *lb_boundaries;

  void lb_init_boundaries();
  
#endif
  /* end of code duplication */

  extern cudaStream_t stream[1];
  extern cudaError_t _err;

  #define PI_FLOAT 3.14159265358979323846f

  EK_parameters ek_parameters = { -1.0, -1.0,    0,
                                     0,    0,    0,
                                  -1.0, -1.0,  0.0,
                                   0.0,  0.0, -1.0,
                                  -1.0,    0
                                };
                                
  static __device__ __constant__ EK_parameters ek_parameters_gpu;
  EK_parameters *ek_parameters_gpu_pointer;
  LB_parameters_gpu *ek_lbparameters_gpu;
  CUDA_particle_data *particle_data_gpu;

  cufftHandle plan_fft;
  cufftHandle plan_ifft;
  
  bool initialized = false;
  
  extern LB_parameters_gpu lbpar_gpu;
  extern LB_node_force_gpu node_f;
  extern LB_nodes_gpu *current_nodes;

#ifdef __cplusplus
}
#endif


__device__ inline void atomicadd( float* address,
                                  float value
                                ) {

#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110

  #warning Using slower atomicAdd emulation
  
  //float-atomic-add from 
  //[url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]
  
  float old = value;
  while( ( old = atomicExch( address, atomicExch( address, 0.0f ) + old ) ) != 0.0f );
  
#else
  #error CUDA compute capability 1.1 or higher required
#endif
}


__device__ unsigned int ek_getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x +
         blockDim.x * blockIdx.x +
         threadIdx.x;
}


__device__ void rhoindex_linear2cartesian( unsigned int index,
                                           unsigned int * coord
                                         ) {

  coord[0]  = index % ek_parameters_gpu.dim_x;
  index    /= ek_parameters_gpu.dim_x;
  coord[1]  = index % ek_parameters_gpu.dim_y;
  coord[2]  = index / ek_parameters_gpu.dim_y;
}


__device__ unsigned int rhoindex_cartesian2linear( unsigned int x,
                                                   unsigned int y,
                                                   unsigned int z
                                                 ) {

  return z * ek_parameters_gpu.dim_y * ek_parameters_gpu.dim_x +
         y * ek_parameters_gpu.dim_x +
         x;
}


__device__ void jindex_linear2cartesian( unsigned int index,
                                         unsigned int * coord,
                                         unsigned int * c
                                       ) {

  coord[0]  = index % ek_parameters_gpu.dim_x;
  index    /= ek_parameters_gpu.dim_x;
  coord[1]  = index % ek_parameters_gpu.dim_y;
  index    /= ek_parameters_gpu.dim_y;
  coord[2]  = index % ek_parameters_gpu.dim_z;
  *c        = index / ek_parameters_gpu.dim_z;
}


__device__ unsigned int jindex_cartesian2linear( unsigned int x,
                                                 unsigned int y,
                                                 unsigned int z,
                                                 unsigned int c
                                               ) {
                                                 
  return c * ek_parameters_gpu.number_of_nodes + 
         z * ek_parameters_gpu.dim_y * ek_parameters_gpu.dim_x +
         y * ek_parameters_gpu.dim_x +
         x;
}


//TODO fluxindex fastest running might improve caching
__device__ unsigned int jindex_getByRhoLinear( unsigned int rho_index,
                                               unsigned int c
                                             ) {
                                               
  return c * ek_parameters_gpu.number_of_nodes +
         rho_index;
}


__device__ void ek_displacement( float * dx,
                                  LB_nodes_gpu n,
                                  unsigned int node_index,
                                  LB_parameters_gpu * ek_lbparameters_gpu
                               ) {
                                 
  float rho = ek_lbparameters_gpu->rho *
              ek_lbparameters_gpu->agrid *
              ek_lbparameters_gpu->agrid *
              ek_lbparameters_gpu->agrid;
  
  rho +=    n.vd[  0 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  1 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  2 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  3 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  4 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  5 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  6 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  7 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  8 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[  9 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 10 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 11 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 12 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 13 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 14 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 15 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 16 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 17 * ek_lbparameters_gpu->number_of_nodes + node_index ] +
            n.vd[ 18 * ek_lbparameters_gpu->number_of_nodes + node_index ];

  dx[0] = ( n.vd[  1 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[  2 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[  7 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[  8 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[  9 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 10 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 11 * ek_lbparameters_gpu->number_of_nodes + node_index ] - 
            n.vd[ 12 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 13 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 14 * ek_lbparameters_gpu->number_of_nodes + node_index ]   );
                 
  dx[1] = ( n.vd[  3 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[  4 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[  7 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[  8 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) -
          ( n.vd[  9 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 10 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 15 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 16 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 17 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 18 * ek_lbparameters_gpu->number_of_nodes + node_index ]   );
          
  dx[2] = ( n.vd[  5 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[  6 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 11 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 12 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) -
          ( n.vd[ 13 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 14 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) +
          ( n.vd[ 15 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 16 * ek_lbparameters_gpu->number_of_nodes + node_index ]   ) -
          ( n.vd[ 17 * ek_lbparameters_gpu->number_of_nodes + node_index ] -
            n.vd[ 18 * ek_lbparameters_gpu->number_of_nodes + node_index ]   );

  dx[0] *= 1.0f / rho;
  dx[1] *= 1.0f / rho;
  dx[2] *= 1.0f / rho;
}


__global__ void ek_calculate_quantities( unsigned int species_index,
                                         LB_nodes_gpu lb_node,
                                         LB_node_force_gpu node_f,
                                         LB_parameters_gpu *ek_lbparameters_gpu
                                       ) {
                                       
  unsigned int index = ek_getThreadIndex ();
  
  if(index < ek_parameters_gpu.number_of_nodes) {
  
    unsigned int coord[3];
    unsigned int neighborindex[9];
    float dx[3];
    int di[3];
    int node;
    float flux, force;
    float boltzmannfactor_local, boltzmannfactor_neighbor;
    
    rhoindex_linear2cartesian( index, coord );
    
    /* Calculate the diffusive fluxes between this node and its neighbors. Only 
       the 9 fluxes along the directions of the LB velocities c_i with i odd are
       stored with a node to avoid redundencies. */
       
    neighborindex[EK_LINK_U00] =
      rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu.dim_x,
         coord[1],
         coord[2]
      );
      
    neighborindex[EK_LINK_0U0] =
      rhoindex_cartesian2linear(
         coord[0],
        (coord[1] + 1) % ek_parameters_gpu.dim_y,
         coord[2]
      );
      
    neighborindex[EK_LINK_00U] =
      rhoindex_cartesian2linear(
         coord[0],
         coord[1],
        (coord[2] + 1) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_UU0] =
      rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu.dim_x,
        (coord[1] + 1) % ek_parameters_gpu.dim_y,
         coord[2]
      );
      
    neighborindex[EK_LINK_UD0] =
      rhoindex_cartesian2linear(
        (coord[0] + 1                          ) % ek_parameters_gpu.dim_x,
        (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
         coord[2]
      );
      
    neighborindex[EK_LINK_U0U] =
      rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu.dim_x,
         coord[1],
        (coord[2] + 1) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_U0D] =
      rhoindex_cartesian2linear(
        (coord[0] + 1                          ) % ek_parameters_gpu.dim_x,
         coord[1],
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_0UU] =
      rhoindex_cartesian2linear(
         coord[0],
        (coord[1] + 1) % ek_parameters_gpu.dim_y,
        (coord[2] + 1) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_0UD] =
      rhoindex_cartesian2linear(
         coord[0],
        (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
    
    
    /* diffusive contribution to flux and LB force*/
    
    /* TODO: take out all of the boltzmann factor based calculations and replace
             them with direct gradient evaluations. */
             
    boltzmannfactor_local = 
      exp( 1.0f / ek_parameters_gpu.T *
           ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
         );
    
    float tune1 = 1.0f; //needs scaling when other fore directions are back on
    float tune2 = 0.5f;
    float tune3 = 0.5f;
    float wallcorrection = 0.0f; //might have to be different
    
    //face in x
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
             ek_parameters_gpu.ext_force[0][species_index] * ek_parameters_gpu.agrid
           )
         );
         
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;
           
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_U00]]); //I think this is shouldn't be there (10.02.2013).

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U00)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]]) / ek_parameters_gpu.agrid * tune3;
//    force *= powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    force += force * wallcorrection * (lb_node.boundary[index] + lb_node.boundary[neighborindex[EK_LINK_U00]] != 0);
//    atomicadd(&node_f.force[index], force);
//    atomicadd(&node_f.force[neighborindex[EK_LINK_U00]], force);

    force  = ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*)ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
               ((cufftReal*)ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;
            
    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;
             
    atomicadd( &node_f.force[index],
                ek_parameters_gpu.rho[species_index][index] *
                ( force / 2 +
                  ek_parameters_gpu.ext_force[0][species_index] *
                  powf(ek_parameters_gpu.agrid, 1) *
                  ek_parameters_gpu.time_step *
                  ek_parameters_gpu.time_step
                )
              );
              
    atomicadd( &node_f.force[neighborindex[EK_LINK_U00]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
                force / 2 );
    
    //face in y
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
             ek_parameters_gpu.ext_force[1][species_index] * ek_parameters_gpu.agrid
           )
         );
         
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;
           
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_0U0]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
               flux * ek_parameters_gpu.time_step );
              
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]]) / ek_parameters_gpu.agrid * tune3;
//    force *= powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + index], force);
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]], force);

    force  = ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;
            
    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;
             
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               ( force / 2 +
                 ek_parameters_gpu.ext_force[1][species_index] *
                 powf(ek_parameters_gpu.agrid, 1) *
                 ek_parameters_gpu.time_step *
                 ek_parameters_gpu.time_step
               )
             );
              
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]],
               ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
               force / 2 );
               
    //face in z
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ( (cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
               ek_parameters_gpu.ext_force[2][species_index] *
               ek_parameters_gpu.agrid
             )
           );
           
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;
           
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_00U]]); //Still think shouldn't be there, but EOF fluctuates more without it. (13.02.2013)

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_00U)],
               flux * ek_parameters_gpu.time_step );
              
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 - (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]]) / ek_parameters_gpu.agrid * tune3;
//    force *= powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    force += force * wallcorrection * (lb_node.boundary[index] + lb_node.boundary[neighborindex[EK_LINK_00U]] != 0);
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + index], force);
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]], force);

    force  = ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;
            
    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;
             
    atomicadd( &node_f.force[2 * ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               ( force / 2 +
                 ek_parameters_gpu.ext_force[2][species_index] *
                 powf(ek_parameters_gpu.agrid, 1) *
                 ek_parameters_gpu.time_step *
                 ek_parameters_gpu.time_step
               )
             );
              
    atomicadd( &node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]],
               ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
               force / 2 );
    
    //edge in z
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UU0]] -
             ( ek_parameters_gpu.ext_force[0][species_index] +
               ek_parameters_gpu.ext_force[1][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
             
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f/boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
           
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_UU0]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
               flux * ek_parameters_gpu.time_step
             );
              
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1; //Pago says the 0.5 goes here. I doubt it.
//    atomicadd(&node_f.force[index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[neighborindex[EK_LINK_UU0]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UU0]], force / sqrt(2.0f));
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UD0]] -
             ( ek_parameters_gpu.ext_force[0][species_index] -
               ek_parameters_gpu.ext_force[1][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_UD0]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + index], -force / sqrt(2.0f));
//    atomicadd(&node_f.force[neighborindex[EK_LINK_UD0]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UD0]], -force / sqrt(2.0f));
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0U]] -
             ( ek_parameters_gpu.ext_force[0][species_index] +
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
             * boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_U0U]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[neighborindex[EK_LINK_U0U]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0U]], force / sqrt(2.0f));
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0D]] -
             ( ek_parameters_gpu.ext_force[0][species_index] -
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_U0D]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + index], -force / sqrt(2.0f));
//    atomicadd(&node_f.force[neighborindex[EK_LINK_U0D]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0D]], -force / sqrt(2.0f));
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UU]] -
             ( ek_parameters_gpu.ext_force[1][species_index] +
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] * boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_0UU]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force / sqrt(2.0f));
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UD]] -
             ( ek_parameters_gpu.ext_force[1][species_index] -
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ek_parameters_gpu.d[species_index] *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
//    flux *= (1 - lb_node.boundary[index]) * (1 - lb_node.boundary[neighborindex[EK_LINK_0UD]]);

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
               flux * ek_parameters_gpu.time_step );
    
//    force = flux / ek_parameters_gpu.d[species_index] * tune2 + (ek_parameters_gpu.rho[species_index][index] - ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]]) / (sqrt(2.0f) * ek_parameters_gpu.agrid) * tune3;
//    force *= 0.5f * powf(ek_parameters_gpu.agrid, 1) * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step * ek_parameters_gpu.T * tune1;
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + index], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + index], -force / sqrt(2.0f));
//    atomicadd(&node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], force / sqrt(2.0f));
//    atomicadd(&node_f.force[2 * ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], -force / sqrt(2.0f));
    
    
    /* advective contribution to flux */
    
    ek_displacement( dx, lb_node, index, ek_lbparameters_gpu );
    
    di[0] = 1 - signbit(dx[0]);
    di[1] = 1 - signbit(dx[1]);
    di[2] = 1 - signbit(dx[2]);

    dx[0] = fabs(dx[0]);
    dx[1] = fabs(dx[1]);
    dx[2] = fabs(dx[2]);
    
    //face in x
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        coord[2]
      );
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_U00 )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * (1.0 - dx[1]) * (1.0 - dx[2])
             );
    
    //face in y
    node =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] + di[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_0U0 )],
              (2 * di[1] - 1) * ek_parameters_gpu.rho[species_index][index] *
              (1.0 - dx[0]) * dx[1] * (1.0 - dx[2]) );
    
    //face in z
    node =
      rhoindex_cartesian2linear(
        coord[0],
        coord[1],
        (coord[2] + di[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_00U )],
               (2 * di[2] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0 - dx[0]) * (1.0 - dx[1]) * dx[2] );
    
    //edge in x
    node =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] + di[1] - 1                   + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] + (1 - di[1]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
        
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_0UU + (di[1] + di[2] == 1) )],
               (2 * di[1] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0 - dx[0]) * dx[1] * dx[2]
             );
    
    //edge in y
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        (coord[2] + (1 - di[0]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_U0U + (di[0] + di[2] == 1) )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * (1.0 - dx[1]) * dx[2] );
    
    //edge in z
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + (1 - di[0]) * (2*di[1] - 1) + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_UU0 + (di[0] + di[1] == 1) )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * dx[1] * (1.0 - dx[2]) );
    
    //corner
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + (1 - di[0]) * (2*di[1] - 1) + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] + (1 - di[0]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    atomicadd( &ek_parameters_gpu.j[
                jindex_getByRhoLinear( node, (1 - di[0]) *
                                             (EK_LINK_UUU + 2*di[1] + di[2]) +
                                             di[0] * (EK_LINK_UDD - 2*di[1] - di[2])
                                     ) ],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * dx[1] * dx[2] );
  }
}


__global__ void ek_propagate_densities( unsigned int species_index
                                      ) {
                                      
  unsigned int index = ek_getThreadIndex();
  
  if( index < ek_parameters_gpu.number_of_nodes ) {
  
    unsigned int neighborindex[13];
    unsigned int coord[3];
    
    rhoindex_linear2cartesian(index, coord);
    
    /* Indices of the neighbors storing the other half
       of the fluxes associated with this link */
    neighborindex[EK_LINK_D00-13] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        coord[2]
      );
      
    neighborindex[EK_LINK_0D0-13] =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    neighborindex[EK_LINK_00D-13] =
      rhoindex_cartesian2linear(
        coord[0],
        coord[1],
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_DD0-13] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    neighborindex[EK_LINK_DU0-13] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    neighborindex[EK_LINK_D0D-13] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_D0U-13] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_0DD-13] =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    neighborindex[EK_LINK_0DU-13] =
      rhoindex_cartesian2linear(
       coord[0],
       (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
       (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
     );
      
    
    neighborindex[EK_LINK_DDD-13] =
      rhoindex_cartesian2linear(
       (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
       (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
       (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
     );
      
    neighborindex[EK_LINK_DDU-13] =
      rhoindex_cartesian2linear(
       (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
       (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
       (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
     );
      
    neighborindex[EK_LINK_DUD-13] = 
      rhoindex_cartesian2linear(
       (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
       (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
       (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
     );
      
    neighborindex[EK_LINK_DUU-13] =
      rhoindex_cartesian2linear(
       (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
       (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
       (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
     );
      
    
    /* Calculate change of densities due to diffusive fluxes */
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_U00 ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_D00-13], EK_LINK_U00 ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_0U0 ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_0D0-13], EK_LINK_0U0 ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_00U ) ];
    ek_parameters_gpu.rho[species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_00D-13], EK_LINK_00U ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UU0 ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DD0-13], EK_LINK_UU0 ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UD0 ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DU0-13], EK_LINK_UD0 ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_U0U ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_D0D-13], EK_LINK_U0U ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_U0D ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_D0U-13], EK_LINK_U0D ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_0UU ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_0DD-13], EK_LINK_0UU ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_0UD ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_0DU-13], EK_LINK_0UD ) ];
    
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UUU ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DDD-13], EK_LINK_UUU ) ];
      
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UUD ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DDU-13], EK_LINK_UUD ) ];
      
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UDU ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DUD-13], EK_LINK_UDU ) ];
      
    ek_parameters_gpu.rho[ species_index ][index] -=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, EK_LINK_UDD ) ];
    ek_parameters_gpu.rho[ species_index ][index] +=
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[EK_LINK_DUU-13], EK_LINK_UDD ) ];
  }
}


__global__ void ek_apply_boundaries( unsigned int species_index,
                                     LB_nodes_gpu lbnode,
                                     LB_node_force_gpu node_f
                                   ) {

  unsigned int index = ek_getThreadIndex();
  unsigned int neighborindex[22];
  unsigned int coord[3];

  if( index < ek_parameters_gpu.number_of_nodes ) {
  
    if( lbnode.boundary[index] ) {
    
      rhoindex_linear2cartesian(index, coord);
      
      /* Indices of the neighbors */
      neighborindex[EK_LINK_D00-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          coord[1],
          coord[2]
        );
          
      neighborindex[EK_LINK_0D0-13] =
        rhoindex_cartesian2linear(
          coord[0], 
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          coord[2]
        );
        
      neighborindex[EK_LINK_00D-13] =
        rhoindex_cartesian2linear(
          coord[0],
          coord[1],
          (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
        );
        
      neighborindex[EK_LINK_DD0-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          coord[2]
        );
        
      neighborindex[EK_LINK_DU0-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
          coord[2]
        );
        
      neighborindex[EK_LINK_D0D-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          coord[1],
          (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
        );
        
      neighborindex[EK_LINK_D0U-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          coord[1],
          (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
        );
        
      neighborindex[EK_LINK_0DD-13] =
        rhoindex_cartesian2linear(
          coord[0],
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
        );
        
      neighborindex[EK_LINK_0DU-13] =
        rhoindex_cartesian2linear(
          coord[0],
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
        );
        
      neighborindex[EK_LINK_DDD-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
        );
      
      neighborindex[EK_LINK_DDU-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
          (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
        );
      
      neighborindex[EK_LINK_DUD-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
          (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
        );
      
      neighborindex[EK_LINK_DUU-13] =
        rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
          (coord[1] + 1                          ) % ek_parameters_gpu.dim_y,
          (coord[2] + 1                          ) % ek_parameters_gpu.dim_z
        );
      
      /* Clear fluxes on links connecting a boundary node */
      for( int i = 0; i < 13; i++ ) {
      
        ek_parameters_gpu.j[jindex_getByRhoLinear(index, i)] = 0.0f;
      }
        
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_D00-13 ], EK_LINK_U00 ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_0D0-13 ], EK_LINK_0U0 ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_00D-13 ], EK_LINK_00U ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DD0-13 ], EK_LINK_UU0 ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DU0-13 ], EK_LINK_UD0 ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_D0D-13 ], EK_LINK_U0U ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_D0U-13 ], EK_LINK_U0D ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_0DD-13 ], EK_LINK_0UU ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_0DU-13 ], EK_LINK_0UD ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DDD-13 ], EK_LINK_UUU ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DDU-13 ], EK_LINK_UUD ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DUD-13 ], EK_LINK_UDU ) ] = 0.0f;
      ek_parameters_gpu.j[ jindex_getByRhoLinear( neighborindex[ EK_LINK_DUU-13 ], EK_LINK_UDD ) ] = 0.0f;
    }
  }
}


//TODO maybe make this obsolete by a multiplication in the advective fluxes, just as it's done for the diffusive ones
__global__ void ek_clear_fluxes() {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) {
  
    for( int i = 0; i < 13; i++ ) {
    
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, i ) ] = 0.0f;
    }
  }
}


__global__ void ek_init_species_density_homogeneous() {

  unsigned int index = ek_getThreadIndex();

  if(index < ek_parameters_gpu.number_of_nodes) {
  
    for(int i = 0; i < ek_parameters_gpu.number_of_species; i++) {
    
      ek_parameters_gpu.rho[ i ][ index ] = ek_parameters_gpu.density[ i ] *
                                            ek_parameters_gpu.agrid *
                                            ek_parameters_gpu.agrid *
                                            ek_parameters_gpu.agrid;
    }
  }
}


__global__ void ek_multiply_greensfcn() {

  unsigned int index = ek_getThreadIndex();
  
  if( index < ek_parameters_gpu.dim_z *
              ek_parameters_gpu.dim_y *
              (ek_parameters_gpu.dim_x / 2 + 1) ) {
  
    ek_parameters_gpu.charge_potential[ index ].x *= ek_parameters_gpu.greensfcn[ index ];
    ek_parameters_gpu.charge_potential[ index ].y *= ek_parameters_gpu.greensfcn[ index ];
  }
}


__global__ void ek_gather_species_charge_density() {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) {
    ((cufftReal*) ek_parameters_gpu.charge_potential)[ index ] = 0.0f;
    
    for( int i = 0; i < ek_parameters_gpu.number_of_species; i++ ) {
    
      ((cufftReal*) ek_parameters_gpu.charge_potential)[ index ] +=
        ek_parameters_gpu.valency[ i ] * ek_parameters_gpu.rho[ i ][ index ] /
        powf( ek_parameters_gpu.agrid, 3 );
    }
  }
}


__global__ void ek_gather_particle_charge_density( CUDA_particle_data * particle_data,
                                                   LB_parameters_gpu * ek_lbparameters_gpu
                                                 ) {

  unsigned int index = ek_getThreadIndex();
  unsigned int lowernode[3];
  float cellpos[3];
  float gridpos;

  if( index < ek_lbparameters_gpu->number_of_particles ) {
  
    gridpos      = particle_data[ index ].p[0] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[0] = (int) floorf( gridpos );
    cellpos[0]   = gridpos - lowernode[0];
  
    gridpos      = particle_data[ index ].p[1] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[1] = (int) floorf( gridpos );
    cellpos[1]   = gridpos - lowernode[1];
  
    gridpos      = particle_data[ index ].p[2] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[2] = (int) floorf( gridpos );
    cellpos[2]   = gridpos - lowernode[2];
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( lowernode[0],
                                            lowernode[1],
                                            lowernode[2]  )
               ],
               particle_data[ index ].q *
               ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] )
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters_gpu.dim_x,
                                            lowernode[1],
                                            lowernode[2]                                    )
               ],
               particle_data[ index ].q *
               cellpos[0] * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] )
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( lowernode[0],
                                            ( lowernode[1] + 1 ) % ek_parameters_gpu.dim_y,
                                            lowernode[2]                                    )
               ],
               particle_data[ index ].q *
               ( 1 - cellpos[0] ) * cellpos[1] * ( 1 - cellpos[2] )
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( lowernode[0],
                                            lowernode[1],
                                            ( lowernode[2] + 1 ) % ek_parameters_gpu.dim_z  )
               ],
               particle_data[ index ].q *
               ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * cellpos[2]
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters_gpu.dim_x,
                                            ( lowernode[1] + 1 ) % ek_parameters_gpu.dim_y,
                                            lowernode[2]                                    )
               ],
               particle_data[ index ].q *
               cellpos[0] * cellpos[1] * ( 1 - cellpos[2] )
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters_gpu.dim_x,
                                            lowernode[1],
                                            ( lowernode[2] + 1 ) % ek_parameters_gpu.dim_z  )
               ],
               particle_data[ index ].q *
               cellpos[0] * ( 1 - cellpos[1] ) * cellpos[2]
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( lowernode[0],
                                            ( lowernode[1] + 1 ) % ek_parameters_gpu.dim_y,
                                            ( lowernode[2] + 1 ) % ek_parameters_gpu.dim_z  )
               ],
               particle_data[ index ].q *
               ( 1 - cellpos[0] ) * cellpos[1] * cellpos[2]
    );
    
    atomicadd( &((cufftReal*) ek_parameters_gpu.charge_potential)[
                 rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters_gpu.dim_x,
                                            ( lowernode[1] + 1 ) % ek_parameters_gpu.dim_y,
                                            ( lowernode[2] + 1 ) % ek_parameters_gpu.dim_z  )
               ],
               particle_data[ index ].q *
               cellpos[0] * cellpos[1] * cellpos[2]
    );
    
    //((cufftReal*) ek_parameters_gpu.charge_potential)[ index ] = 0.0f;
    printf("particle %d (%d):\n  charge %f\n  pos %f %f %f\n  lowernode %d %d %d\n  cellpos %f %f %f\n\n", index, ek_lbparameters_gpu->number_of_particles, particle_data[index].q, particle_data[index].p[0], particle_data[index].p[1], particle_data[index].p[2], lowernode[0], lowernode[1], lowernode[2], cellpos[0], cellpos[1], cellpos[2]); //TODO delete
  }
}


__global__ void ek_create_greensfcn() {

  unsigned int index = ek_getThreadIndex();
  unsigned int tmp;
  unsigned int coord[3];
  
  coord[0] = index % ( ek_parameters_gpu.dim_x / 2 + 1 );
  tmp      = index / ( ek_parameters_gpu.dim_x / 2 + 1 );
  coord[1] = tmp % ek_parameters_gpu.dim_y;
  coord[2] = tmp / ek_parameters_gpu.dim_y;
  
  if( index < ek_parameters_gpu.dim_z *
              ek_parameters_gpu.dim_y *
              ( ek_parameters_gpu.dim_x / 2 + 1 ) ) {
              
    if( index == 0 ) {
    
      //setting 0th fourier mode to 0 enforces charge neutrality
      ek_parameters_gpu.greensfcn[index] = 0.0f;
    }
    else {
    
      ek_parameters_gpu.greensfcn[ index ] =
        -4.0f * PI_FLOAT * ek_parameters_gpu.bjerrumlength *
        ek_parameters_gpu.T * ek_parameters_gpu.agrid * ek_parameters_gpu.agrid *
        0.5f /
        ( cos( 2.0f * PI_FLOAT * coord[0] / (cufftReal) ek_parameters_gpu.dim_x ) +
          cos( 2.0f * PI_FLOAT * coord[1] / (cufftReal) ek_parameters_gpu.dim_y ) +
          cos( 2.0f * PI_FLOAT * coord[2] / (cufftReal) ek_parameters_gpu.dim_z ) -
          3.0f
        ) /
        ( ek_parameters_gpu.dim_x *
          ek_parameters_gpu.dim_y *
          ek_parameters_gpu.dim_z
        );
    }
  }
}


__global__ void ek_clear_boundary_densities( LB_nodes_gpu lbnode ) {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) {
  
    if( lbnode.boundary[ index ] ) {
    
      for( int i = 0; i < ek_parameters_gpu.number_of_species; i++ ) {
      
        ek_parameters_gpu.rho[ i ][ index ] = 0.0f;
      }
    }
  }
}


//TODO delete
__global__ void ek_clear_node_force( LB_node_force_gpu node_f ) {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) {
  
    node_f.force[ index ]                                         = 0.0f;
    node_f.force[ ek_parameters_gpu.number_of_nodes + index ]     = 0.0f;
    node_f.force[ 2 * ek_parameters_gpu.number_of_nodes + index ] = 0.0f;
  }
}


#ifdef __cplusplus
extern "C" {
#endif

void ek_integrate_electrostatics() {

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 ) /
    ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_gather_species_charge_density, dim_grid, threads_per_block, () );
  
  if ( lbpar_gpu.number_of_particles != 0 ) { //TODO make it an if number_of_charged_particles != 0
  
    blocks_per_grid_x =
      ( lbpar_gpu.number_of_particles + threads_per_block * blocks_per_grid_y - 1 ) /
      ( threads_per_block * blocks_per_grid_y );
    dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
    
    particle_data_gpu = gpu_get_particle_pointer();
    
    KERNELCALL( ek_gather_particle_charge_density,
                dim_grid, threads_per_block,
                ( particle_data_gpu, ek_lbparameters_gpu ) );
  }
  
  if( cufftExecR2C( plan_fft,
                    (cufftReal*) ek_parameters.charge_potential,
                    ek_parameters.charge_potential               ) != CUFFT_SUCCESS ) {
                    
    fprintf(stderr, "ERROR: Unable to execute FFT plan\n");
  }
  
  blocks_per_grid_x =
    ( ek_parameters.dim_z * ek_parameters.dim_y * ( ek_parameters.dim_x / 2 + 1 ) +
      threads_per_block * blocks_per_grid_y - 1) / 
    ( threads_per_block * blocks_per_grid_y );
  dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_multiply_greensfcn, dim_grid, threads_per_block, () );
    
  if( cufftExecC2R( plan_ifft,
                    ek_parameters.charge_potential,
                    (cufftReal*) ek_parameters.charge_potential ) != CUFFT_SUCCESS ) {
                    
    fprintf(stderr, "ERROR: Unable to execute iFFT plan\n");
  }
}


void ek_integrate() {

  /** values for the kernel call */
  
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
    / (threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  //TODO delete
  KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, ( node_f ) );
  
  /* Integrate diffusion-advection */
  
  for( int i = 0; i < ek_parameters.number_of_species; i++ ) {
  
    KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
    KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                ( i, *current_nodes, node_f, ek_lbparameters_gpu )    );
              
#ifdef EK_BOUNDARIES
    KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                ( i, *current_nodes, node_f )                     );
#endif

    KERNELCALL( ek_propagate_densities, dim_grid, threads_per_block, ( i ) );
  }
  
  /* Integrate electrostatics */
  
  ek_integrate_electrostatics();
  
  /* Integrate Navier-Stokes */
  
  lb_integrate_GPU();
  
  //TODO delete - needed for printfs
  cudaDeviceSynchronize();
}


#ifdef EK_BOUNDARIES
void ek_init_species_density_wallcharge( float* wallcharge_species_density,
                                         int wallcharge_species             ) {
  
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
    / ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_init_species_density_homogeneous, dim_grid, threads_per_block, () );
  KERNELCALL( ek_clear_boundary_densities, dim_grid, threads_per_block, ( *current_nodes ) );
  
  if( wallcharge_species != -1 ) {
  
    cuda_safe_mem( cudaMemcpy( ek_parameters.rho[wallcharge_species], wallcharge_species_density,
                               ek_parameters.number_of_nodes * sizeof( float ),
                               cudaMemcpyHostToDevice                                             ) );
  }
}
#endif


void ek_init_species( int species ) {

  if( !initialized ) {
  
    ek_init();
  }
  
  if( ek_parameters.species_index[ species ] == -1 ) {
  
    ek_parameters.species_index[ species ] = ek_parameters.number_of_species;
    ek_parameters.number_of_species++;
    
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                               ek_parameters.number_of_nodes * sizeof( float )                        ) );
    
    ek_parameters.density[      ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.D[            ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.valency[      ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[0][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[1][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[2][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.d[            ek_parameters.species_index[ species ] ] =
      ek_parameters.D[          ek_parameters.species_index[ species ] ] / ( 1.0 + 2.0 * sqrt( 2.0 ) );
  }
}


int ek_init() {

  if( ek_parameters.agrid < 0.0 ||
      ek_parameters.viscosity < 0.0 ||
      ek_parameters.T < 0.0 ||
      ek_parameters.bjerrumlength < 0.0 ) {
      
    fprintf( stderr, "ERROR: invalid agrid, viscosity, T or bjerrum_length\n" );
    
    return 1;
  }
    
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x;
  dim3 dim_grid;
  
  if(!initialized) {
  
    if( cudaGetSymbolAddress( (void**) &ek_parameters_gpu_pointer, ek_parameters_gpu ) != cudaSuccess) {
    
      fprintf( stderr, "ERROR: Fetching constant memory pointer\n" );
      
      return 1;
    }
    
    for( int i = 0; i < MAX_NUMBER_OF_SPECIES; i++ ) {
    
      ek_parameters.species_index[i] = -1;
    }
      
    lattice_switch = LATTICE_LB_GPU;
    lbpar_gpu.agrid = ek_parameters.agrid;
    lbpar_gpu.viscosity = ek_parameters.viscosity;
    lbpar_gpu.bulk_viscosity = ek_parameters.bulk_viscosity;
    lbpar_gpu.friction = ek_parameters.friction;
    
    lbpar_gpu.rho = 1.0;
    lbpar_gpu.external_force = 0;
    lbpar_gpu.ext_force[0] = 0.0;
    lbpar_gpu.ext_force[1] = 0.0;
    lbpar_gpu.ext_force[2] = 0.0;
    
    lb_init_gpu();

    ek_parameters.dim_x = lbpar_gpu.dim_x;
    ek_parameters.dim_y = lbpar_gpu.dim_y;
    ek_parameters.dim_z = lbpar_gpu.dim_z;
    ek_parameters.time_step = lbpar_gpu.time_step;
    ek_parameters.number_of_nodes = ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z;

    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.j,
                             ek_parameters.number_of_nodes * 13 * sizeof( float ) ) );
    cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
    
    lb_get_para_pointer( &ek_lbparameters_gpu );
    lb_set_ek_pointer( ek_parameters_gpu_pointer );
    
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.charge_potential,
                             sizeof( cufftComplex ) *
                             ek_parameters.dim_z * ek_parameters.dim_y * ( ek_parameters.dim_x / 2 + 1 ) ) );
    
    if( cudaGetLastError() != cudaSuccess ) {
    
        fprintf(stderr, "ERROR: Failed to allocate\n");
        
        return 1;
    }
    
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.greensfcn,
                             sizeof( cufftReal ) * 
                ek_parameters.dim_z * ek_parameters.dim_y * ( ek_parameters.dim_x / 2 + 1 ) ) );
    
    if( cudaGetLastError() != cudaSuccess ) {
    
        fprintf(stderr, "ERROR: Failed to allocate\n");
        
        return 1;
    }
    
    cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
    
    blocks_per_grid_x =
      ( ek_parameters.dim_z * ek_parameters.dim_y * (ek_parameters.dim_x / 2 + 1) +
        threads_per_block * blocks_per_grid_y - 1
      ) / ( threads_per_block * blocks_per_grid_y );
    dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
    KERNELCALL( ek_create_greensfcn, dim_grid, threads_per_block, () );

    /* create 3D FFT plans */
    
    if( cufftPlan3d( &plan_fft,
                     ek_parameters.dim_z,
                     ek_parameters.dim_y,
                     ek_parameters.dim_x,
                     CUFFT_R2C            ) != CUFFT_SUCCESS ) {
    
        fprintf(stderr, "ERROR: Unable to create fft plan\n");
        return 1;
    }
    
    if( cufftSetCompatibilityMode( plan_fft, CUFFT_COMPATIBILITY_NATIVE ) != CUFFT_SUCCESS ) {
    
        fprintf(stderr, "ERROR: Unable to set fft compatibility mode to native\n");
        return 1;
    }
    
    if( cufftSetStream( plan_fft, stream[0]) != CUFFT_SUCCESS ) {
    
        fprintf(stderr, "ERROR: Unable to assign FFT to cuda stream\n");
        return 1;
    }

    if( cufftPlan3d( &plan_ifft,
                     ek_parameters.dim_z,
                     ek_parameters.dim_y,
                     ek_parameters.dim_x,
                     CUFFT_C2R            ) != CUFFT_SUCCESS ) {
    
        fprintf(stderr, "ERROR: Unable to create ifft plan\n");
        return 1;
    }
    
    if( cufftSetCompatibilityMode( plan_ifft, CUFFT_COMPATIBILITY_NATIVE ) != CUFFT_SUCCESS) {
    
        fprintf(stderr, "ERROR: Unable to set ifft compatibility mode to native\n");
        return 1;
    }
    
    if( cufftSetStream( plan_ifft, stream[0] ) != CUFFT_SUCCESS ) {
    
        fprintf(stderr, "ERROR: Unable to assign FFT to cuda stream\n");
        return 1;
    }
    
    initialized = true;
  }
  
  cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
  
#ifdef EK_BOUNDARIES
  lb_init_boundaries();
#else
  blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
    / (threads_per_block * blocks_per_grid_y );
  dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_init_species_density_homogeneous, dim_grid, threads_per_block, () );
#endif

  ek_integrate_electrostatics();
  
  //ek_print_parameters(); //TODO delete

  return 0;
}


int ek_lb_print_vtk_velocity( char* filename ) {
  FILE* fp = fopen( filename, "w" );
	
  if( fp == NULL ) {
  
  	return 1;
	}
  
  LB_values_gpu *host_values = (LB_values_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_values_gpu )     );
  lb_get_values_GPU( host_values );
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
velocity\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\nPOINT_DATA %u\n\
SCALARS velocity float 3\n\
LOOKUP_TABLE default\n",
           lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
           lbpar_gpu.agrid / 2, lbpar_gpu.agrid / 2, lbpar_gpu.agrid / 2,
           lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
           lbpar_gpu.number_of_nodes                                      );

  for( int i = 0; i < lbpar_gpu.number_of_nodes; i++ ) {
  
    fprintf( fp, "%f %f %f ", host_values[ i ].v[0],
                              host_values[ i ].v[1],
                              host_values[ i ].v[2]  );
  }
  
  free(host_values);	
  fclose(fp);
  
	return 0;
}


int ek_lb_print_vtk_density( char* filename ) {

  FILE* fp = fopen( filename, "w" );
	
  if( fp == NULL ) {
  
  	return 1;
	}
  
  LB_values_gpu *host_values = (LB_values_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_values_gpu )     );
  lb_get_values_GPU( host_values );
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
density_lb\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS density_lb float 1\n\
LOOKUP_TABLE default\n",
           lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
           lbpar_gpu.agrid / 2, lbpar_gpu.agrid / 2, lbpar_gpu.agrid / 2,
           lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
           lbpar_gpu.number_of_nodes                                      );

  for( int i = 0; i < lbpar_gpu.number_of_nodes; i++ ) {
  
    fprintf( fp, "%f ", host_values[ i ].rho );
  }
  
  free( host_values );	
  fclose( fp );
  
	return 0;
}


int ek_print_vtk_density( int species, char* filename ) {

  FILE* fp = fopen( filename, "w" );
	
  if( fp == NULL )
  	return 1;
  	
  float* densities = (float*) malloc( ek_parameters.number_of_nodes *
                                      sizeof( float )                 );
  
  if( ek_parameters.species_index[ species ] != -1 ) {
  
    cudaMemcpy( densities, ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                ek_parameters.number_of_nodes * sizeof( float ),
                cudaMemcpyDeviceToHost                                                  );
  }
  else
    return 1;
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
density_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS density_%d float 1\n\
LOOKUP_TABLE default\n",
           species,
           ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
           ek_parameters.agrid / 2, ek_parameters.agrid / 2, ek_parameters.agrid / 2,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes,
           species                                                                    );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) {
  
    fprintf( fp, "%f\n", densities[ i ] );
  }
  
  free( densities );
  fclose( fp );
  
	return 0;
}


int ek_print_vtk_potential( char* filename ) {

  FILE* fp = fopen( filename, "w" );
	
  if( fp == NULL ) {
  
  	return 1;
	}
  	
  float* potential = (float*) malloc( ek_parameters.number_of_nodes * sizeof( cufftReal ) );
  
  cudaMemcpy( potential, ek_parameters.charge_potential,
              ek_parameters.number_of_nodes * sizeof( cufftReal ),
              cudaMemcpyDeviceToHost                               );
  
  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
potential\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS potential float 1\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid / 2, ek_parameters.agrid / 2, ek_parameters.agrid / 2,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) {
  
    fprintf( fp, "%f\n", potential[ i ] );
  }
  
  free( potential );	
  fclose( fp );
  
	return 0;
}


int ek_print_vtk_lbforce( char* filename ) {

  FILE* fp = fopen( filename, "w" );
	
  if( fp == NULL ) {
  
  	return 1;
	}
  	
  float* lbforce = (float*) malloc( ek_parameters.number_of_nodes * 3 *sizeof( float ) );
  
  cudaMemcpy( lbforce, node_f.force,
              ek_parameters.number_of_nodes * 3 * sizeof( float ),
              cudaMemcpyDeviceToHost                               );
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
lbforce\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS lbforce float 3\n\
LOOKUP_TABLE default\n",
           ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
           ek_parameters.agrid / 2, ek_parameters.agrid / 2, ek_parameters.agrid / 2,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) {
  
    fprintf( fp, "%f %f %f\n", lbforce[ i ],
                              lbforce[ i + ek_parameters.number_of_nodes ],
                              lbforce[ i + 2 * ek_parameters.number_of_nodes ] );
  }
  
  free( lbforce );
  fclose( fp );
  
	return 0;
}


void ek_print_parameters() {

  printf( "ek_parameters {\n" );
  
  printf( "  float agrid = %f;\n",                    ek_parameters.agrid );
  printf( "  float time_step = %f;\n",                ek_parameters.time_step );
  printf( "  unsigned int dim_x = %d;\n",             ek_parameters.dim_x );
  printf( "  unsigned int dim_y = %d;\n",             ek_parameters.dim_y );
  printf( "  unsigned int dim_z = %d;\n",             ek_parameters.dim_z );
  printf( "  unsigned int number_of_nodes = %d;\n",   ek_parameters.number_of_nodes );
  printf( "  float viscosity = %f;\n",                ek_parameters.viscosity );
  printf( "  float bulk_viscosity = %f;\n",           ek_parameters.bulk_viscosity );
  printf( "  float gamma_odd = %f;\n",                ek_parameters.gamma_odd );
  printf( "  float gamma_even = %f;\n",               ek_parameters.gamma_even );
  printf( "  float friction = %f;\n",                 ek_parameters.friction );
  printf( "  float T = %f;\n",                        ek_parameters.T );
  printf( "  float bjerrumlength = %f;\n",            ek_parameters.bjerrumlength );
  printf( "  unsigned int number_of_species = %d;\n", ek_parameters.number_of_species);  
  printf( "  float* j = %p;\n",                       ek_parameters.j );
  
  printf( "  float* rho[] = {%p, %p, %p, %p, %p, %p, %p, %p, %p, %p};\n",
          ek_parameters.rho[0], ek_parameters.rho[1], ek_parameters.rho[2],
          ek_parameters.rho[3], ek_parameters.rho[4], ek_parameters.rho[5],
          ek_parameters.rho[6], ek_parameters.rho[7], ek_parameters.rho[8],
          ek_parameters.rho[9]                                              );
  
  printf( "  int species_index[] = {%d, %d, %d, %d, %d, %d, %d, %d, %d, %d};\n",
          ek_parameters.species_index[0], ek_parameters.species_index[1],
          ek_parameters.species_index[2], ek_parameters.species_index[3],
          ek_parameters.species_index[4], ek_parameters.species_index[5],
          ek_parameters.species_index[6], ek_parameters.species_index[7],
          ek_parameters.species_index[8], ek_parameters.species_index[9]         );
  
  printf( "  float density = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.density[0], ek_parameters.density[1],
          ek_parameters.density[2], ek_parameters.density[3],
          ek_parameters.density[4], ek_parameters.density[5],
          ek_parameters.density[6], ek_parameters.density[7],
          ek_parameters.density[8], ek_parameters.density[9]                );
  
  printf( "  float D[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.D[0], ek_parameters.D[1], ek_parameters.D[2],
          ek_parameters.D[3], ek_parameters.D[4], ek_parameters.D[5],
          ek_parameters.D[6], ek_parameters.D[7], ek_parameters.D[8],
          ek_parameters.D[9]                                           );
  
  printf( "  float d[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.d[0], ek_parameters.d[1], ek_parameters.d[2],
          ek_parameters.d[3], ek_parameters.d[4], ek_parameters.d[5],
          ek_parameters.d[6], ek_parameters.d[7], ek_parameters.d[8],
          ek_parameters.d[9]                                                   );
  
  printf( "  float valency[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.valency[0], ek_parameters.valency[1], ek_parameters.valency[2],
          ek_parameters.valency[3], ek_parameters.valency[4], ek_parameters.valency[5],
          ek_parameters.valency[6], ek_parameters.valency[7], ek_parameters.valency[8],
          ek_parameters.valency[9]                                                      );
  
  printf( "  float ext_force[0][] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.ext_force[0][0], ek_parameters.ext_force[0][1], ek_parameters.ext_force[0][2],
          ek_parameters.ext_force[0][3], ek_parameters.ext_force[0][4], ek_parameters.ext_force[0][5],
          ek_parameters.ext_force[0][6], ek_parameters.ext_force[0][7], ek_parameters.ext_force[0][8],
          ek_parameters.ext_force[0][9]                                                                );
  
  printf( "  float ext_force[1][] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.ext_force[1][0], ek_parameters.ext_force[1][1], ek_parameters.ext_force[1][2],
          ek_parameters.ext_force[1][3], ek_parameters.ext_force[1][4], ek_parameters.ext_force[1][5],
          ek_parameters.ext_force[1][6], ek_parameters.ext_force[1][7], ek_parameters.ext_force[1][8],
          ek_parameters.ext_force[1][9]                                                                );
  
  printf( "  float ext_force[2][] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
          ek_parameters.ext_force[2][0], ek_parameters.ext_force[2][1], ek_parameters.ext_force[2][2],
          ek_parameters.ext_force[2][3], ek_parameters.ext_force[2][4], ek_parameters.ext_force[2][5],
          ek_parameters.ext_force[2][6], ek_parameters.ext_force[2][7], ek_parameters.ext_force[2][8],
          ek_parameters.ext_force[2][9]                                                                );
  
  printf( "}\n" );
}


void ek_print_lbpar() {

  printf("lbpar_gpu {\n");
  
  printf("    float rho = %f;\n",                        lbpar_gpu.rho );
  printf("    float mu = %f;\n",                         lbpar_gpu.mu );
  printf("    float viscosity = %f;\n",                  lbpar_gpu.viscosity );
  printf("    float gamma_shear = %f;\n",                lbpar_gpu.gamma_shear );
  printf("    float gamma_bulk = %f;\n",                 lbpar_gpu.gamma_bulk );
  printf("    float gamma_odd = %f;\n",                  lbpar_gpu.gamma_odd );
  printf("    float gamma_even = %f;\n",                 lbpar_gpu.gamma_even );
  printf("    float agrid = %f;\n",                      lbpar_gpu.agrid );
  printf("    float tau = %f;\n",                        lbpar_gpu.tau );
  printf("    float friction = %f;\n",                   lbpar_gpu.friction );
  printf("    float time_step = %f;\n",                  lbpar_gpu.time_step );
  printf("    float lb_coupl_pref = %f;\n",              lbpar_gpu.lb_coupl_pref );
  printf("    float lb_coupl_pref2 = %f;\n",             lbpar_gpu.lb_coupl_pref2 );
  printf("    float bulk_viscosity = %f;\n",             lbpar_gpu.bulk_viscosity );
  printf("    unsigned int dim_x = %d;\n",               lbpar_gpu.dim_x );
  printf("    unsigned int dim_y = %d;\n",               lbpar_gpu.dim_y );
  printf("    unsigned int dim_z = %d;\n",               lbpar_gpu.dim_z );
  printf("    unsigned int number_of_nodes = %d;\n",     lbpar_gpu.number_of_nodes );
  printf("    unsigned int number_of_particles = %d;\n", lbpar_gpu.number_of_particles );
  printf("    int fluct = %d;\n",                        lbpar_gpu.fluct );
  printf("    int calc_val = %d;\n",                     lbpar_gpu.calc_val );
  printf("    int external_force = %d;\n",               lbpar_gpu.external_force );
  printf("    float ext_force[3] = {%f, %f, %f};\n",     lbpar_gpu.ext_force[0],
                                                         lbpar_gpu.ext_force[1],
                                                         lbpar_gpu.ext_force[2] );
  printf("    unsigned int your_seed = %d;\n",           lbpar_gpu.your_seed );
  printf("    unsigned int reinit = %d;\n",              lbpar_gpu.reinit );
  
  printf("}\n");
}


int ek_set_agrid( double agrid ) {  

  if( ek_parameters.agrid < 0.0 ) {
  
    ek_parameters.agrid = agrid;
    
    return 0;
  }
  else {
  
    printf("ERROR: electrokinetics agrid can not be changed\n");
    
    return 1;
  }
}


int ek_set_bjerrumlength( double bjerrumlength ) {

  if( ek_parameters.bjerrumlength < 0.0 ) {
  
    ek_parameters.bjerrumlength = bjerrumlength;
    
    return 0;
  }
  else {
  
    printf("ERROR: electrokinetics bjerrum_length can not be changed\n");
    
    return 1;
  }
}


int ek_set_viscosity( double viscosity ) {

  ek_parameters.viscosity = viscosity;
  
  return 0;
}


int ek_set_friction( double friction ) {

  ek_parameters.friction = friction;
  
  return 0;
}


int ek_set_bulk_viscosity( double bulk_viscosity ) {

  ek_parameters.bulk_viscosity = bulk_viscosity;
  
  return 0;
}


int ek_set_gamma_odd( double gamma_odd ) {

  ek_parameters.gamma_odd = gamma_odd;
  
  return 0;
}


int ek_set_gamma_even( double gamma_even ) {

  ek_parameters.gamma_even = gamma_even;
  
  return 0;
}


int ek_set_density( int species, double density ) {

  ek_init_species( species );
  ek_parameters.density[ ek_parameters.species_index[ species ] ] = density;
  
  lbpar_gpu.rho = 0.0;
  
  for( int i = 0; i < MAX_NUMBER_OF_SPECIES; i++ ) {
  
    lbpar_gpu.rho += ek_parameters.density[i];
  }
  
  lb_reinit_parameters_gpu();
  
  return 0;
}


int ek_set_D( int species, double D ) {

  ek_init_species( species );
  
  ek_parameters.D[ ek_parameters.species_index[ species ] ] = D;
  ek_parameters.d[ ek_parameters.species_index[ species ] ] = D / ( 1.0 + 2.0 * sqrt(2.0)) ;
  
  return 0;
}

int ek_set_T(double T) {

  ek_parameters.T = T;
  
  return 0;
}


int ek_set_valency( int species, double valency ) {

  ek_init_species( species );
  
  ek_parameters.valency[ ek_parameters.species_index[ species ] ] = valency;
  
  return 0;
}


int ek_set_ext_force( int species,
                      double ext_force_x,
                      double ext_force_y,
                      double ext_force_z
                    ) {
                    
  ek_init_species( species );
  
  ek_parameters.ext_force[0][ ek_parameters.species_index[ species ] ] = ext_force_x;
  ek_parameters.ext_force[1][ ek_parameters.species_index[ species ] ] = ext_force_y;
  ek_parameters.ext_force[2][ ek_parameters.species_index[ species ] ] = ext_force_z;
  
  return 0;
}

#ifdef __cplusplus
}
#endif

#endif /* ELECTROKINETICS */
