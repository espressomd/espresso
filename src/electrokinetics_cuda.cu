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
#include "config.hpp"
#ifdef CUDA /* Terminates at end of file */



#include <cuda.h>
#include <cufft.h>



#include <stdio.h>
#include "lb-boundaries.hpp"
#include "electrokinetics.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "lbgpu.hpp"
#include "constraint.hpp"



#ifdef ELECTROKINETICS /* Terminates at end of file */

  /* TODO: get rid of this code duplication with lb-boundaries.h by solving the
           cuda-mpi incompatibility */

#define LATTICE_OFF      0
#define LATTICE_LB_CPU   1
#define LATTICE_LB_GPU   2
extern int lattice_switch;
extern int ek_initialized;
extern EK_parameters* lb_ek_parameters_gpu;

// Used to limit register use for the pressure calculation
#define EK_LINK_U00_pressure 0
#define EK_LINK_0U0_pressure 1
#define EK_LINK_00U_pressure 2
#define EK_LINK_D00_pressure 3
#define EK_LINK_0D0_pressure 4
#define EK_LINK_00D_pressure 5
     
#ifdef EK_BOUNDARIES
  extern int n_lb_boundaries;
  extern LB_Boundary *lb_boundaries;

  void lb_init_boundaries();
#endif
  /* end of code duplication */

  extern cudaStream_t stream[1];
  extern cudaError_t _err;

  #define PI_FLOAT 3.14159265358979323846f

  EK_parameters ek_parameters = { -1.0, -1.0, -1.0,
                                     0,    0,    0,
                                     0,
                                  -1.0, -1.0,  0.0,
                                   0.0,  0.0, -1.0,
                                  -1.0,
                                     0,    0,
                                  -1.0, -1.0, -1.0, 
                                  {0.0,  0.0, 0.0},
                                  { -1,   -1,  -1},
                                  -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0,
                                  -1.0
                                };
                                
  static __device__ __constant__ EK_parameters ek_parameters_gpu;
  static __device__ float ek_accelerated_frame_boundary_force [3] = { 0.0f, 0.0f, 0.0f };
  static float* ek_accelerated_frame_boundary_velocity;
  EK_parameters *ek_parameters_gpu_pointer;
  LB_parameters_gpu *ek_lbparameters_gpu;
  CUDA_particle_data *particle_data_gpu;
  float *ek_lb_boundary_force;
  char *ek_node_is_catalyst;
  unsigned int old_number_of_species = 0;
  unsigned int old_number_of_boundaries = 0;

  cufftHandle plan_fft;
  cufftHandle plan_ifft;
  
  bool initialized = false;
  
  extern LB_parameters_gpu lbpar_gpu;
  extern LB_node_force_gpu node_f;
  extern LB_nodes_gpu *current_nodes;
  extern EK_parameters *lb_ek_parameters;
  
  LB_rho_v_gpu *ek_lb_device_values;



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
                                 
  float rho = ek_lbparameters_gpu->rho[0] *
              ek_lbparameters_gpu->agrid *
              ek_lbparameters_gpu->agrid *
              ek_lbparameters_gpu->agrid;

  float mode [19];

  for ( int i = 0; i < 19; i++ )
  {
    mode[i] = n.vd[  i * ek_lbparameters_gpu->number_of_nodes + node_index ];
  }
  
  rho += mode[  0 ] +
         mode[  1 ] +
         mode[  2 ] +
         mode[  3 ] +
         mode[  4 ] +
         mode[  5 ] +
         mode[  6 ] +
         mode[  7 ] +
         mode[  8 ] +
         mode[  9 ] +
         mode[ 10 ] +
         mode[ 11 ] +
         mode[ 12 ] +
         mode[ 13 ] +
         mode[ 14 ] +
         mode[ 15 ] +
         mode[ 16 ] +
         mode[ 17 ] +
         mode[ 18 ];

  dx[0] = ( mode[  1 ] - mode[  2 ] ) +
          ( mode[  7 ] - mode[  8 ] ) +
          ( mode[  9 ] - mode[ 10 ] ) +
          ( mode[ 11 ] - mode[ 12 ] ) +
          ( mode[ 13 ] - mode[ 14 ] );
                 
  dx[1] = ( mode[  3 ] - mode[  4 ] ) +
          ( mode[  7 ] - mode[  8 ] ) -
          ( mode[  9 ] - mode[ 10 ] ) +
          ( mode[ 15 ] - mode[ 16 ] ) +
          ( mode[ 17 ] - mode[ 18 ] );
          
  dx[2] = ( mode[  5 ] - mode[  6 ] ) +
          ( mode[ 11 ] - mode[ 12 ] ) -
          ( mode[ 13 ] - mode[ 14 ] ) +
          ( mode[ 15 ] - mode[ 16 ] ) -
          ( mode[ 17 ] - mode[ 18 ] );

  // Velocity requires half the force in the previous time step

  dx[0] += 0.5f * ek_parameters_gpu.lb_force_previous[ node_index ];
  dx[1] += 0.5f * ek_parameters_gpu.lb_force_previous[ ek_parameters_gpu.number_of_nodes + node_index ];
  dx[2] += 0.5f * ek_parameters_gpu.lb_force_previous[ 2 * ek_parameters_gpu.number_of_nodes + node_index ];

  dx[0] *= 1.0f / rho;
  dx[1] *= 1.0f / rho;
  dx[2] *= 1.0f / rho;
}

#ifdef EK_REACTION
__global__ void ek_pressure(
                             LB_nodes_gpu n_a,
                             LB_parameters_gpu *ek_lbparameters_gpu,
                             LB_rho_v_gpu *d_v
                           )
{
  unsigned int index = ek_getThreadIndex ();

  if( index < ek_parameters_gpu.number_of_nodes )
  {  
    ek_parameters_gpu.pressure[ index ] = 0.0f;
 
    // Add the ideal-gas contribution f from the EK
    // species, which is given by n_i * k. In MD units
    // the proper expression is n_i * T / ag^2, where 
    // there is a 1/ag^3 factor coming from converting the
    // internal EK particle number back to a density,
    // and an ag factor that is required to get the 
    // proper pressure difference

    for ( int i = 0; i < ek_parameters_gpu.number_of_species; i++ )
    {
      ek_parameters_gpu.pressure[ index ] += ek_parameters_gpu.rho[ i ][ index ] *
                                             ek_parameters_gpu.T /
                                             powf(ek_parameters_gpu.agrid, 2);
    }

    // Set pressure to zero inside boundary

    ek_parameters_gpu.pressure[ index ] *= (n_a.boundary[index] == 0);
  }
}

__global__ void ek_add_ideal_pressure_to_lb_force(
                                                   LB_nodes_gpu lb_node,
                                                   LB_node_force_gpu node_f,
                                                   LB_parameters_gpu *ek_lbparameters_gpu
                                                 ) 
{
  unsigned int coord[3];
  unsigned int neighborindex[6];
  unsigned int index = ek_getThreadIndex ();

  if(index < ek_parameters_gpu.number_of_nodes)
  {
    float pressure_gradient;
    
    rhoindex_linear2cartesian( index, coord );

    // Calculate the indices of the neighbours to which
    // the force is to be applied
       
    neighborindex[EK_LINK_U00_pressure] =
      rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu.dim_x,
         coord[1],
         coord[2]
      );
      
    neighborindex[EK_LINK_0U0_pressure] =
      rhoindex_cartesian2linear(
         coord[0],
        (coord[1] + 1) % ek_parameters_gpu.dim_y,
         coord[2]
      );
      
    neighborindex[EK_LINK_00U_pressure] =
      rhoindex_cartesian2linear(
         coord[0],
         coord[1],
        (coord[2] + 1) % ek_parameters_gpu.dim_z
      );

    neighborindex[EK_LINK_D00_pressure] =
      rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
         coord[1],
         coord[2]
      );
      
    neighborindex[EK_LINK_0D0_pressure] =
      rhoindex_cartesian2linear(
         coord[0],
        (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
         coord[2]
      );
      
    neighborindex[EK_LINK_00D_pressure] =
      rhoindex_cartesian2linear(
         coord[0],
         coord[1],
        (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );

    // Force in x direction (multiplicative factor
    // comes from converting MD force into LB force)

    pressure_gradient = (   ek_parameters_gpu.pressure[ neighborindex[EK_LINK_D00_pressure] ]
                          - ek_parameters_gpu.pressure[ neighborindex[EK_LINK_U00_pressure] ] )/
                        ( 2.0f * ek_parameters_gpu.agrid );

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 3) *
                         ek_parameters_gpu.time_step *
                         ek_parameters_gpu.time_step;

    pressure_gradient *= ( (   lb_node.boundary[ neighborindex[EK_LINK_U00_pressure] ]
                             + lb_node.boundary[ index ]
                             + lb_node.boundary[ neighborindex[EK_LINK_D00_pressure] ] ) == 0 );

    atomicadd( &node_f.force[index], pressure_gradient );
    
    // Force in y direction

    pressure_gradient = (   ek_parameters_gpu.pressure[ neighborindex[EK_LINK_0D0_pressure] ]
                          - ek_parameters_gpu.pressure[ neighborindex[EK_LINK_0U0_pressure] ] )/
                        ( 2.0f * ek_parameters_gpu.agrid );

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 3) *
                         ek_parameters_gpu.time_step *
                         ek_parameters_gpu.time_step;

    pressure_gradient *= ( (   lb_node.boundary[ neighborindex[EK_LINK_0U0_pressure] ]
                             + lb_node.boundary[ index ]
                             + lb_node.boundary[ neighborindex[EK_LINK_0D0_pressure] ] ) == 0 );

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], pressure_gradient );
              
    // Force in z direction

    pressure_gradient = (   ek_parameters_gpu.pressure[ neighborindex[EK_LINK_00D_pressure] ]
                          - ek_parameters_gpu.pressure[ neighborindex[EK_LINK_00U_pressure] ] )/
                        ( 2.0f * ek_parameters_gpu.agrid );

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 3) *
                         ek_parameters_gpu.time_step *
                         ek_parameters_gpu.time_step;

    pressure_gradient *= ( (   lb_node.boundary[ neighborindex[EK_LINK_00U_pressure] ]
                             + lb_node.boundary[ index ]
                             + lb_node.boundary[ neighborindex[EK_LINK_00D_pressure] ] ) == 0 );

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], pressure_gradient );
  }
}
#endif

__global__ void ek_accelerated_frame_transformation( LB_node_force_gpu node_f )
{
                                       
  unsigned int index = ek_getThreadIndex ();

  if(index < ek_parameters_gpu.number_of_nodes)
  {
    //face in x

    atomicadd( &node_f.force[index],
               ek_accelerated_frame_boundary_force[0] );
    
    //face in y

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index],
               ek_accelerated_frame_boundary_force[1] );
              
    //face in z

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index],
               ek_accelerated_frame_boundary_force[2] );
  }
}

__global__ void ek_calculate_quantities( unsigned int species_index,
                                         LB_nodes_gpu lb_node,
                                         LB_node_force_gpu node_f,
                                         LB_parameters_gpu *ek_lbparameters_gpu,
                                         LB_rho_v_gpu *d_v
                                       ) {
                                       
  unsigned int index = ek_getThreadIndex ();

  if(index < ek_parameters_gpu.number_of_nodes)
  {
  
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
 
    //face in x
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
             ek_parameters_gpu.ext_force[0][species_index] * ek_parameters_gpu.agrid
           )
         );
         
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;
           
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U00)],
               flux * ek_parameters_gpu.time_step );

    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*)ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
               ((cufftReal*)ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;

    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;
             
    atomicadd( &node_f.force[index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force / 2.0f +
                 ek_parameters_gpu.ext_force[0][species_index] *
                 (
                   powf(ek_parameters_gpu.agrid, 1) *
                   ek_parameters_gpu.time_step *
                   ek_parameters_gpu.time_step
                 )
               )
             );

    atomicadd( &node_f.force[neighborindex[EK_LINK_U00]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
                force / 2.0f );
    
    //face in y
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
             ek_parameters_gpu.ext_force[1][species_index] * ek_parameters_gpu.agrid
           )
         );
         
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;
           
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
               flux * ek_parameters_gpu.time_step );
              
    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;

    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force / 2.0f +
                 ek_parameters_gpu.ext_force[1][species_index] *
                 (
                   powf(ek_parameters_gpu.agrid, 1) *
                   ek_parameters_gpu.time_step *
                   ek_parameters_gpu.time_step
                 )
               )
             );

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
                force / 2.0f );
              
    //face in z
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ( (cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
               ek_parameters_gpu.ext_force[2][species_index] *
               ek_parameters_gpu.agrid
             )
           );
           
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
             boltzmannfactor_neighbor
           ) / ek_parameters_gpu.agrid;

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_00U)],
               flux * ek_parameters_gpu.time_step );

    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) / ek_parameters_gpu.agrid;

    force *= powf(ek_parameters_gpu.agrid, 1) *
             ek_parameters_gpu.time_step *
             ek_parameters_gpu.time_step;

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force / 2.0f +
                 ek_parameters_gpu.ext_force[2][species_index] *
                 (
                   powf(ek_parameters_gpu.agrid, 1) *
                   ek_parameters_gpu.time_step *
                   ek_parameters_gpu.time_step
                 )
               )
             );

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
                force / 2.0f );
    
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
             
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f/boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
           
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
               flux * ek_parameters_gpu.time_step
             );

    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UD0]] -
             ( ek_parameters_gpu.ext_force[0][species_index] -
               ek_parameters_gpu.ext_force[1][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
               flux * ek_parameters_gpu.time_step );
        
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0U]] -
             ( ek_parameters_gpu.ext_force[0][species_index] +
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
             * boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
               flux * ek_parameters_gpu.time_step );
        
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0D]] -
             ( ek_parameters_gpu.ext_force[0][species_index] -
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
               flux * ek_parameters_gpu.time_step );
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UU]] -
             ( ek_parameters_gpu.ext_force[1][species_index] +
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] * boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
               flux * ek_parameters_gpu.time_step );
    
    boltzmannfactor_neighbor =
      exp( 1.0f / ek_parameters_gpu.T *
           ( ek_parameters_gpu.valency[species_index] *
             ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UD]] -
             ( ek_parameters_gpu.ext_force[1][species_index] -
               ek_parameters_gpu.ext_force[2][species_index]
             ) * ek_parameters_gpu.agrid
           )
         );
    
    flux = ( ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid ) *
           ( 1.0f / boltzmannfactor_local +
             1.0f / boltzmannfactor_neighbor
           ) / 2.0f *
           ( ek_parameters_gpu.rho[species_index][index] *
             boltzmannfactor_local -
             ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]] *
             boltzmannfactor_neighbor
           ) /
           ( sqrt(2.0f) * ek_parameters_gpu.agrid );

    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
               flux * ek_parameters_gpu.time_step );    
    
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
               dx[0] * (1.0f - dx[1]) * (1.0f - dx[2])
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
              (1.0f - dx[0]) * dx[1] * (1.0f - dx[2]) );
    
    //face in z
    node =
      rhoindex_cartesian2linear(
        coord[0],
        coord[1],
        (coord[2] + di[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_00U )],
               (2 * di[2] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0f - dx[0]) * (1.0f - dx[1]) * dx[2] );
    
    //edge in x
    node =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] + di[1] - 1                   + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] + (1 - di[1]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
        
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_0UU + (di[1] + di[2] == 1) )],
               (2 * di[1] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0f - dx[0]) * dx[1] * dx[2]
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
               dx[0] * (1.0f - dx[1]) * dx[2] );
    
    //edge in z
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + (1 - di[0]) * (2*di[1] - 1) + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_UU0 + (di[0] + di[1] == 1) )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * dx[1] * (1.0f - dx[2]) );
    
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
  
  if( index < ek_parameters_gpu.number_of_nodes ) 
  {
  
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
    ek_parameters_gpu.rho[species_index ][index]  +=
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

  if( index < ek_parameters_gpu.number_of_nodes ) 
  {
    if( lbnode.boundary[index] ) 
    {
    
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
      for( int i = 0; i < 13; i++ )
        ek_parameters_gpu.j[jindex_getByRhoLinear(index, i)] = 0.0f;
        
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

  if( index < ek_parameters_gpu.number_of_nodes ) 
  {
    for( int i = 0; i < 13; i++ ) 
    {
      ek_parameters_gpu.j[ jindex_getByRhoLinear( index, i ) ] = 0.0f;
    }
  }
}


__global__ void ek_init_species_density_homogeneous() {

  unsigned int index = ek_getThreadIndex();

  if(index < ek_parameters_gpu.number_of_nodes) 
  {  
    for(int i = 0; i < ek_parameters_gpu.number_of_species; i++) 
    {
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
              (ek_parameters_gpu.dim_x / 2 + 1) ) 
  {
    ek_parameters_gpu.charge_potential[ index ].x *= ek_parameters_gpu.greensfcn[ index ];
    ek_parameters_gpu.charge_potential[ index ].y *= ek_parameters_gpu.greensfcn[ index ];
  }
}


__global__ void ek_gather_species_charge_density() {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) 
  {
    ((cufftReal*) ek_parameters_gpu.charge_potential)[ index ] = 0.0f;
    
    for( int i = 0; i < ek_parameters_gpu.number_of_species; i++ ) 
    {
    
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

  if( index < ek_lbparameters_gpu->number_of_particles ) 
  {  
    gridpos      = particle_data[ index ].p[0] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[0] = (int) floorf( gridpos );
    cellpos[0]   = gridpos - lowernode[0];
  
    gridpos      = particle_data[ index ].p[1] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[1] = (int) floorf( gridpos );
    cellpos[1]   = gridpos - lowernode[1];
  
    gridpos      = particle_data[ index ].p[2] / ek_parameters_gpu.agrid - 0.5f;
    lowernode[2] = (int) floorf( gridpos );
    cellpos[2]   = gridpos - lowernode[2];

    lowernode[0] = (lowernode[0] + ek_lbparameters_gpu->dim_x) % ek_lbparameters_gpu->dim_x;
    lowernode[1] = (lowernode[1] + ek_lbparameters_gpu->dim_y) % ek_lbparameters_gpu->dim_y;
    lowernode[2] = (lowernode[2] + ek_lbparameters_gpu->dim_z) % ek_lbparameters_gpu->dim_z;

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
    //printf("particle %d (%d):\n  charge %f\n  pos %f %f %f\n  lowernode %d %d %d\n  cellpos %f %f %f\n\n", index, ek_lbparameters_gpu->number_of_particles, particle_data[index].q, particle_data[index].p[0], particle_data[index].p[1], particle_data[index].p[2], lowernode[0], lowernode[1], lowernode[2], cellpos[0], cellpos[1], cellpos[2]); //TODO delete
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
              ( ek_parameters_gpu.dim_x / 2 + 1 ) ) 
  {
              
    if( index == 0 ) 
    {
    
      //setting 0th fourier mode to 0 enforces charge neutrality
      ek_parameters_gpu.greensfcn[index] = 0.0f;
    }
    else 
    {
    
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

  if( index < ek_parameters_gpu.number_of_nodes ) 
  {  
    if( lbnode.boundary[ index ] ) 
    {
    
      for( int i = 0; i < ek_parameters_gpu.number_of_species; i++ ) 
      {     
        ek_parameters_gpu.rho[ i ][ index ] = 0.0f;
      }
    }
  }
}


//TODO delete ?? (it has the previous step setting now)
__global__ void ek_clear_node_force( LB_node_force_gpu node_f ) {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes )
  {
    ek_parameters_gpu.lb_force_previous[ index ] = 
                           node_f.force[ index ];
    ek_parameters_gpu.lb_force_previous[ ek_parameters_gpu.number_of_nodes + index ] =
                           node_f.force[ ek_parameters_gpu.number_of_nodes + index ];
    ek_parameters_gpu.lb_force_previous[ 2 * ek_parameters_gpu.number_of_nodes + index ] = 
                           node_f.force[ 2 * ek_parameters_gpu.number_of_nodes + index ];

    node_f.force[ index ]                                         = 0.0f;
    node_f.force[ ek_parameters_gpu.number_of_nodes + index ]     = 0.0f;
    node_f.force[ 2 * ek_parameters_gpu.number_of_nodes + index ] = 0.0f;
  }
}


#ifdef EK_REACTION
__global__ void ek_reaction( ) {

  unsigned int index = ek_getThreadIndex();
  unsigned int coord[3];

  float* rho_reactant = &ek_parameters_gpu.rho[ek_parameters_gpu.reaction_species[0]][index];
  float* rho_product0 = &ek_parameters_gpu.rho[ek_parameters_gpu.reaction_species[1]][index];
  float* rho_product1 = &ek_parameters_gpu.rho[ek_parameters_gpu.reaction_species[2]][index];

  float dt = ek_parameters_gpu.time_step;
  float ct_rate = ek_parameters_gpu.reaction_ct_rate;
  float fraction_0 = ek_parameters_gpu.reaction_fraction_0;
  float fraction_1 = ek_parameters_gpu.reaction_fraction_1;

  float rho_change = *rho_reactant * ( 1.0f - expf(-dt*ct_rate) );

  rhoindex_linear2cartesian(index, coord);

  if ( index < ek_parameters_gpu.number_of_nodes )
  {
    if ( ek_parameters_gpu.node_is_catalyst[index] == 1 )
    {
      *rho_reactant -= rho_change;
      *rho_product0 += rho_change * fraction_0;
      *rho_product1 += rho_change * fraction_1;
    }
    else if ( ek_parameters_gpu.node_is_catalyst[index] == 2 )
    {
      *rho_reactant = ek_parameters_gpu.rho_reactant_reservoir * powf(ek_parameters_gpu.agrid,3);
      *rho_product0 = ek_parameters_gpu.rho_product0_reservoir * powf(ek_parameters_gpu.agrid,3);
      *rho_product1 = ek_parameters_gpu.rho_product1_reservoir * powf(ek_parameters_gpu.agrid,3); 
    } 
  }
}
#endif


__global__ void ek_calculate_boundary_forces( int n_lb_boundaries, 
                                              float* ek_lb_boundary_force, 
                                              float* ek_accelerated_frame_boundary_velocity, 
                                              LB_parameters_gpu *ek_lbparameters_gpu )
{
  // Set force to zero

  ek_accelerated_frame_boundary_force[0] = 0.0f;
  ek_accelerated_frame_boundary_force[1] = 0.0f;
  ek_accelerated_frame_boundary_force[2] = 0.0f;

  if ( ek_parameters_gpu.accelerated_frame_enabled == 1 )
  {
    // The number of nodes belonging to the fluid

    float number_of_fluid_nodes = ek_parameters_gpu.accelerated_frame_fluid_mass / 
                                  ek_lbparameters_gpu->rho[0] / 
                                  powf(ek_lbparameters_gpu->agrid,3);

    // Calculate the correct external force, on the Tcl level a volume density
    // is specified, which only acts on the fluid nodes!

    double ext_acc_force[3];
    ext_acc_force[0] = ek_parameters_gpu.ext_acceleration_force[0] * number_of_fluid_nodes *  
                       powf(ek_lbparameters_gpu->agrid,3);
    ext_acc_force[1] = ek_parameters_gpu.ext_acceleration_force[1] * number_of_fluid_nodes *  
                       powf(ek_lbparameters_gpu->agrid,3);
    ext_acc_force[2] = ek_parameters_gpu.ext_acceleration_force[2] * number_of_fluid_nodes *  
                       powf(ek_lbparameters_gpu->agrid,3);

    // Add/Start with the external force applied to the particle, which is
    // technically a boundary composite, this has a positive sign because it is
    // applied to the particle

    ek_accelerated_frame_boundary_force[0] += ext_acc_force[0];
    ek_accelerated_frame_boundary_force[1] += ext_acc_force[1];
    ek_accelerated_frame_boundary_force[2] += ext_acc_force[2];

// TODO : REMOVE
//printf("ext_force %f ", ek_accelerated_frame_boundary_force[2]);

    for ( int i = 0; i < n_lb_boundaries; i++)
    {
// TODO : REMOVE
//printf("bndry_force %f ", -ek_lb_boundary_force[3*i + 2]);

      // Sum over all the boundaries that make up the composite and add the total
      // friction force to the external force vector wise. The boundary force
      // coming from the LB has the incorrect sign, it points in the direction
      // of the fluid, not in the direction of the boundary, hence the minus.

      ek_accelerated_frame_boundary_force[0] += -ek_lb_boundary_force[3*i + 0];
      ek_accelerated_frame_boundary_force[1] += -ek_lb_boundary_force[3*i + 1];
      ek_accelerated_frame_boundary_force[2] += -ek_lb_boundary_force[3*i + 2];
    }

// TODO : REMOVE
//printf("ext_force+bndry_force %f ", ek_accelerated_frame_boundary_force[2]);

    // INTERMEZZO: Calculate the velocity of the frame, and thus obtain the 
    // velocity of the boundary with respect to the center of mass of the system.
    // In the finite box, the fluid's center of mass moves in the opposite
    // direction.

    ek_accelerated_frame_boundary_velocity[0] += ( ek_accelerated_frame_boundary_force[0] / 
                                                   ek_parameters_gpu.accelerated_frame_boundary_mass *
                                                   ek_parameters_gpu.time_step );
    ek_accelerated_frame_boundary_velocity[1] += ( ek_accelerated_frame_boundary_force[1] / 
                                                   ek_parameters_gpu.accelerated_frame_boundary_mass *
                                                   ek_parameters_gpu.time_step );
    ek_accelerated_frame_boundary_velocity[2] += ( ek_accelerated_frame_boundary_force[2] / 
                                                   ek_parameters_gpu.accelerated_frame_boundary_mass *
                                                   ek_parameters_gpu.time_step );

    // Now calculate the acceleration on the particle by dividing the total force
    // on the particle by the boundary mass. This acceleration is applied on the
    // fluid in the transformation (with a minus sign) and must be multiplied 
    // by the fluid mass to get the right force on the fluid
  
    ek_accelerated_frame_boundary_force[0] *= - (   ek_parameters_gpu.accelerated_frame_fluid_mass
                                                  / ek_parameters_gpu.accelerated_frame_boundary_mass );
    ek_accelerated_frame_boundary_force[1] *= - (   ek_parameters_gpu.accelerated_frame_fluid_mass
                                                  / ek_parameters_gpu.accelerated_frame_boundary_mass );
    ek_accelerated_frame_boundary_force[2] *= - (   ek_parameters_gpu.accelerated_frame_fluid_mass
                                                  / ek_parameters_gpu.accelerated_frame_boundary_mass );

// TODO : REMOVE
//printf("(ef+bf)*(mf/mp) %f ", ek_accelerated_frame_boundary_force[2]);

    // In a finite system there is also always the negative of the external force 
    // on the particle acting on the fluid. This force ensures that the there is
    // no total force acting on the system, hence momentum is conserved. That is,
    // for a moving particle there would otherwise be no stationary state.

    ek_accelerated_frame_boundary_force[0] -= ext_acc_force[0];
    ek_accelerated_frame_boundary_force[1] -= ext_acc_force[1];
    ek_accelerated_frame_boundary_force[2] -= ext_acc_force[2];

// TODO : REMOVE
//printf("(ef+bf)*(mf/mp) + ef %f ", ek_accelerated_frame_boundary_force[2]);

    // Do the unit conversion from LB units (boundary force) to units that 
    // can be added back into the LB fluid via the LB external force (MD units),
    // apparently this requires scaling by the fluid node number

    ek_accelerated_frame_boundary_force[0] *= ( ek_parameters_gpu.agrid * powf(ek_parameters_gpu.time_step, 2) / 
                                                number_of_fluid_nodes );
    ek_accelerated_frame_boundary_force[1] *= ( ek_parameters_gpu.agrid * powf(ek_parameters_gpu.time_step, 2) / 
                                                number_of_fluid_nodes );
    ek_accelerated_frame_boundary_force[2] *= ( ek_parameters_gpu.agrid * powf(ek_parameters_gpu.time_step, 2) / 
                                                number_of_fluid_nodes );
  }
}



void ek_integrate_electrostatics() {

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 ) /
    ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_gather_species_charge_density, dim_grid, threads_per_block, () );
  
  if ( lbpar_gpu.number_of_particles != 0 ) //TODO make it an if number_of_charged_particles != 0
  { 
  
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
                    ek_parameters.charge_potential               ) != CUFFT_SUCCESS ) 
  {
                    
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
                    (cufftReal*) ek_parameters.charge_potential ) != CUFFT_SUCCESS )
  {
                    
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



  /* Clears the force on the nodes and must be called before fluxes are calculated,
     since in the reaction set up the previous-step LB force is added to the flux
     (in ek_calculate_quantities / ek_displacement), which is copied in this routine;
     it also needs to be called before ek_accelerated_frame_transformation */

  KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, ( node_f ) );



#ifdef EK_REACTION
  if ( ek_parameters.reaction_species[0] != -1 &&
       ek_parameters.reaction_species[1] != -1 &&
       ek_parameters.reaction_species[2] != -1 )
  {
    /* Performs the catalytic reaction and sets the reservoir densities at 
       the boundary of the simulation box */

    KERNELCALL( ek_reaction, dim_grid, threads_per_block, ());

    /* Determines the excess pressure that follows from the creation of 
       species by the reaction */

    KERNELCALL( ek_pressure, dim_grid, threads_per_block, ( *current_nodes, 
                                                            ek_lbparameters_gpu, 
                                                            ek_lb_device_values ) );
  }
#endif


#ifdef EK_BOUNDARIES
  if ( ek_parameters.accelerated_frame_enabled == 1 && n_lb_boundaries > 0 )
  {
    /* Adds the force required to perform the accelerated frame tranformation,
       must be done before the integration of the LB is called; force is 
       technically from the previous step and is calculated in the kernel
       ek_calculate_boundary_forces */

    KERNELCALL( ek_accelerated_frame_transformation, dim_grid, threads_per_block, ( node_f ) );
  }
#endif


  /* Integrate diffusion-advection */
  
  for( int i = 0; i < ek_parameters.number_of_species; i++ )
  {
  
    KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
    KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                ( i, *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values ) );
              
#ifdef EK_BOUNDARIES
    KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                ( i, *current_nodes, node_f ) );
#endif

    KERNELCALL( ek_propagate_densities, dim_grid, threads_per_block, ( i ) );
  }



#ifdef EK_REACTION
  if ( ek_parameters.reaction_species[0] != -1 &&
       ek_parameters.reaction_species[1] != -1 &&
       ek_parameters.reaction_species[2] != -1 )
  {
    /* Add pressure force to LB must be done outside of loop,
       otherwise the force gets added several times */

    KERNELCALL( ek_add_ideal_pressure_to_lb_force, dim_grid, threads_per_block,
                  ( *current_nodes, node_f, ek_lbparameters_gpu ) );
  }
#endif



  /* Integrate electrostatics */
  
  ek_integrate_electrostatics();
  
  /* Integrate Navier-Stokes */
  
  lb_integrate_GPU();


  
#ifdef EK_BOUNDARIES
  if ( ek_parameters.accelerated_frame_enabled == 1 && n_lb_boundaries > 0 )
  {
    /* Calculate the total force on the boundaries for the accelerated frame transformation,
       can only be done after the LB integration is called */

    ek_calculate_boundary_forces<<<1,1>>>( n_lb_boundaries, ek_lb_boundary_force, ek_accelerated_frame_boundary_velocity, ek_lbparameters_gpu );
  }
#endif


  
  //TODO delete - needed for printfs
  cudaDeviceSynchronize();

// TODO : REMOVE
/*
LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_rho_v_pi_gpu ) );
lb_get_values_GPU( host_values ); 
printf( "ve %e %e %e\n", host_values[ 0 ].v[0], host_values[ 0 ].v[1], host_values[ 0 ].v[2] );
free(host_values);
*/
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
  
  if( wallcharge_species != -1 ) 
  {  
    cuda_safe_mem( cudaMemcpy( ek_parameters.rho[wallcharge_species], 
                               wallcharge_species_density,
                               ek_parameters.number_of_nodes * sizeof( float ),
                               cudaMemcpyHostToDevice )
                 );
  }
}
#endif


void ek_init_species( int species ) {

  if( !initialized ) 
  {  
    ek_init();
  }
  
  if( ek_parameters.species_index[ species ] == -1 ) 
  {  
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
      ek_parameters.bjerrumlength < 0.0 ) 
  {
      
    fprintf( stderr, "ERROR: invalid agrid, viscosity, T or bjerrum_length\n" );
    
    return 1;
  }
    
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x;
  dim3 dim_grid;
  
  if(!initialized) 
  {
    if( cudaGetSymbolAddress( (void**) &ek_parameters_gpu_pointer, ek_parameters_gpu ) != cudaSuccess) 
    {
      fprintf( stderr, "ERROR: Fetching constant memory pointer\n" );

      return 1;
    }
    
    for( int i = 0; i < MAX_NUMBER_OF_SPECIES; i++ ) 
    {    
      ek_parameters.species_index[i] = -1;
    }

    if ( lattice_switch != LATTICE_OFF ) 
    {
      fprintf( stderr, "ERROR: Electrokinetics automatically intializes the LB on the GPU and can therefore not be used in conjunction with LB.\n");
      fprintf( stderr, "ERROR: Please run either electrokinetics or LB.\n");
      
      return 1;
    }

    lattice_switch = LATTICE_LB_GPU;
    ek_initialized = 1;         

    lbpar_gpu.agrid = ek_parameters.agrid;
    lbpar_gpu.viscosity[0] = ek_parameters.viscosity;
    lbpar_gpu.bulk_viscosity[0] = ek_parameters.bulk_viscosity;
    lbpar_gpu.friction[0] = ek_parameters.friction;

    lbpar_gpu.rho[0] = ( ek_parameters.lb_density < 0.0 ? 1.0 : ek_parameters.lb_density );
    lb_reinit_parameters_gpu();

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

    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.lb_force_previous,
                             ek_parameters.number_of_nodes * 3 * sizeof( float ) ) );

#ifdef EK_REACTION
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.pressure,
                             ek_parameters.number_of_nodes * sizeof( float ) ) );
    ek_node_is_catalyst = (char*) calloc( ek_parameters.number_of_nodes , sizeof( char ) );
#endif

    lb_get_device_values_pointer( &ek_lb_device_values );
    
    if( cudaGetLastError() != cudaSuccess ) 
    {
      fprintf(stderr, "ERROR: Failed to allocate\n");
      return 1;
    }
    
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.greensfcn,
                             sizeof( cufftReal ) * 
                ek_parameters.dim_z * ek_parameters.dim_y * ( ek_parameters.dim_x / 2 + 1 ) ) );
    
    if( cudaGetLastError() != cudaSuccess ) 
    {
      fprintf(stderr, "ERROR: Failed to allocate\n");
      return 1;
    }

    cudaMallocHost((void**) &ek_parameters.node_is_catalyst,
                             sizeof( char ) * 
                ek_parameters.dim_z*ek_parameters.dim_y*ek_parameters.dim_x );
    
    if(cudaGetLastError() != cudaSuccess) 
    {
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
                     CUFFT_R2C            ) != CUFFT_SUCCESS ) 
    {
      fprintf(stderr, "ERROR: Unable to create fft plan\n");
      return 1;
    }
    
    if( cufftSetCompatibilityMode( plan_fft, CUFFT_COMPATIBILITY_NATIVE ) != CUFFT_SUCCESS ) 
    {    
      fprintf(stderr, "ERROR: Unable to set fft compatibility mode to native\n");
      return 1;
    }
    
    if( cufftSetStream( plan_fft, stream[0]) != CUFFT_SUCCESS ) 
    {
        fprintf(stderr, "ERROR: Unable to assign FFT to cuda stream\n");
        return 1;
    }

    if( cufftPlan3d( &plan_ifft,
                     ek_parameters.dim_z,
                     ek_parameters.dim_y,
                     ek_parameters.dim_x,
                     CUFFT_C2R            ) != CUFFT_SUCCESS ) 
    {   
      fprintf(stderr, "ERROR: Unable to create ifft plan\n");
      return 1;
    }
    
    if( cufftSetCompatibilityMode( plan_ifft, CUFFT_COMPATIBILITY_NATIVE ) != CUFFT_SUCCESS) 
    {   
      fprintf(stderr, "ERROR: Unable to set ifft compatibility mode to native\n");
      return 1;
    }
    
    if( cufftSetStream( plan_ifft, stream[0] ) != CUFFT_SUCCESS )
    {    
      fprintf(stderr, "ERROR: Unable to assign FFT to cuda stream\n");
      return 1;
    }
    
    blocks_per_grid_x =
      ( ek_parameters.dim_z * ek_parameters.dim_y * (ek_parameters.dim_x ) +
        threads_per_block * blocks_per_grid_y - 1
      ) / ( threads_per_block * blocks_per_grid_y );
    dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
    KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, ( node_f ) );
    
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.charge_potential,
                             sizeof( cufftComplex ) *
                             ek_parameters.dim_z * ek_parameters.dim_y * ( ek_parameters.dim_x / 2 + 1 ) ) );

    cuda_safe_mem( cudaMalloc( (void**) &ek_accelerated_frame_boundary_velocity,
                               3 * sizeof( float ) ) );

    cuda_safe_mem( cudaMemset( ek_accelerated_frame_boundary_velocity,
                               0,
                               3 * sizeof( float ) ) );

    initialized = true;

    cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
  }
  else
  {
    if ( lbpar_gpu.agrid != ek_parameters.agrid ||
         lbpar_gpu.viscosity[0] != ek_parameters.viscosity ||
         lbpar_gpu.bulk_viscosity[0] != ek_parameters.bulk_viscosity ||
         lbpar_gpu.friction[0] != ek_parameters.friction ||
         ( ( lbpar_gpu.rho[0] != 1.0 ) && ( lbpar_gpu.rho[0] != ek_parameters.lb_density ) )
       )
    {
      fprintf( stderr, "ERROR: The LB parameters on the GPU cannot be reinitialized.\n");
      
      return 1;
    }
    else
    {
      cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );

#ifdef EK_BOUNDARIES
      if ( old_number_of_boundaries != n_lb_boundaries )
      {
        lb_init_boundaries();
        lb_get_boundary_force_pointer( &ek_lb_boundary_force );
        old_number_of_boundaries = n_lb_boundaries;
      }

      // Determine the total boundary mass and the fluid mass

      unsigned int number_of_boundary_nodes = ek_calculate_boundary_mass( );
      unsigned int number_of_fluid_nodes = ek_parameters.number_of_nodes - number_of_boundary_nodes;
      
      ek_parameters.accelerated_frame_boundary_mass =   static_cast<float>( number_of_boundary_nodes )
                                                      * ek_parameters.accelerated_frame_boundary_mass_density
                                                      * powf(lbpar_gpu.agrid,3);

      ek_parameters.accelerated_frame_fluid_mass =   static_cast<float>( number_of_fluid_nodes )
                                                   * lbpar_gpu.rho[0] * powf(lbpar_gpu.agrid,3);
        
      cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
#else
      blocks_per_grid_x =
        ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
        / (threads_per_block * blocks_per_grid_y );
      dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
      

      if ( old_number_of_species != ek_parameters.number_of_species )
      {
        KERNELCALL( ek_init_species_density_homogeneous, dim_grid, threads_per_block, () );
        old_number_of_species = ek_parameters.number_of_species;
      }
#endif

#ifdef EK_REACTION
      // added to ensure that the pressure is set to the proper value in the first time step
      blocks_per_grid_x = (ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) / (threads_per_block * blocks_per_grid_y);
      dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
      KERNELCALL( ek_pressure, dim_grid, threads_per_block, ( *current_nodes, ek_lbparameters_gpu, ek_lb_device_values ) );
#endif

      ek_integrate_electrostatics();
    }
  }

  //ek_print_parameters(); //TODO delete
      
  return 0; 
}


void lb_set_ek_pointer(EK_parameters* pointeradress) {
  lb_ek_parameters_gpu = pointeradress;
}


unsigned int ek_calculate_boundary_mass( )
{
  unsigned int* bound_array = (unsigned int*) malloc( lbpar_gpu.number_of_nodes*sizeof(unsigned int) );

  lb_get_boundary_flags_GPU(bound_array);

  unsigned int boundary_node_number = 0;

  for( int j=0; j<ek_parameters.number_of_nodes; j++)
    if( bound_array[j] != 0 ) boundary_node_number++;

  free(bound_array);

  return boundary_node_number;
}


void rhoindex_linear2cartesian_host( unsigned int index,
                                     unsigned int * coord
                                   ) {

  coord[0]  = index % ek_parameters.dim_x;
  index    /= ek_parameters.dim_x;
  coord[1]  = index % ek_parameters.dim_y;
  coord[2]  = index / ek_parameters.dim_y;
}


unsigned int jindex_cartesian2linear_host( unsigned int x,
                                           unsigned int y,
                                           unsigned int z,
                                           unsigned int c
                                         ) {
                                    
  x = ( x + ek_parameters.dim_x ) % ek_parameters.dim_x; //this does not happen in the GPU version of this function
  y = ( y + ek_parameters.dim_y ) % ek_parameters.dim_y;
  z = ( z + ek_parameters.dim_z ) % ek_parameters.dim_z;
  
  return c * ek_parameters.number_of_nodes + 
         z * ek_parameters.dim_y * ek_parameters.dim_x +
         y * ek_parameters.dim_x +
         x;
}


unsigned int jindex_getByRhoLinear_host( unsigned int rho_index,
                                         unsigned int c
                                       ) {
                                               
  return c * ek_parameters.number_of_nodes +
         rho_index;
}


unsigned int rhoindex_cartesian2linear_host( unsigned int x,
                                             unsigned int y,
                                             unsigned int z
                                           ) {

  return z * ek_parameters.dim_y * ek_parameters.dim_x +
         y * ek_parameters.dim_x +
         x;
}


int ek_lb_print_vtk_velocity( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {  
    return 1;
  }
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_rho_v_pi_gpu ) );
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
           lbpar_gpu.agrid*0.5f, lbpar_gpu.agrid*0.5f, lbpar_gpu.agrid*0.5f,
           lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
           lbpar_gpu.number_of_nodes                                      );

  for( int i = 0; i < lbpar_gpu.number_of_nodes; i++ ) 
  {  
    fprintf( fp, "%e %e %e ", host_values[ i ].v[0],
                              host_values[ i ].v[1],
                              host_values[ i ].v[2]  );
  }
  
  free(host_values);
  fclose(fp);
  
  return 0;
}


int ek_node_print_velocity( int x, int y, int z, double* velocity ) { //TODO only calculate single node velocity
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_rho_v_pi_gpu ) );
  lb_get_values_GPU( host_values );
  
  int i = z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;
  
  velocity[0] = host_values[i].v[0];
  velocity[1] = host_values[i].v[1];
  velocity[2] = host_values[i].v[2];
  
  free(host_values);
  
  return 0;
}


int ek_node_print_mass_flux( int x, int y, int z, double* mass_flux ) { //TODO only calculate single node mass_flux

  float current_mass_fraction;
  float total_mass =   ek_parameters.rho_reactant_reservoir*ek_parameters.mass_reactant
                     + ek_parameters.rho_product0_reservoir*ek_parameters.mass_product0
                     + ek_parameters.rho_product1_reservoir*ek_parameters.mass_product1;

  unsigned int coord[3];
  float flux_local_cartesian[3]; // temporary variable for converting fluxes into cartesian coordinates for output
  float mass_flux_local_cartesian[3] = {0.0,0.0,0.0}; // total mass flux

  float* fluxes = (float*) malloc( ek_parameters.number_of_nodes * 13 * sizeof( float ) );
  
  for ( int species_number = 0; species_number < 3; species_number++ )
  {
    int species = ek_parameters.reaction_species[species_number];

    if ( species_number == 0 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_reactant/total_mass;
    else if ( species_number == 1 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_product0/total_mass;
    else if ( species_number == 2 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_product1/total_mass;
    else 
      return 1;

    if( ek_parameters.species_index[ species ] != -1 ) {
    
      int threads_per_block = 64;
      int blocks_per_grid_y = 4;
      int blocks_per_grid_x =
        ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
        / (threads_per_block * blocks_per_grid_y );
      dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
      
      KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
      KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                  ( ek_parameters.species_index[ species ], *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values ) );
                
#ifdef EK_BOUNDARIES
      KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                  ( ek_parameters.species_index[ species ], *current_nodes, node_f ) );
#endif
    
      cuda_safe_mem( cudaMemcpy( fluxes, 
                                 ek_parameters.j,
                                 ek_parameters.number_of_nodes * 13*sizeof( float ),
                                 cudaMemcpyDeviceToHost )
                   );
    }
    else
      return 1;

    int i = z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;

    rhoindex_linear2cartesian_host(i, coord);
     
    flux_local_cartesian[0]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U00) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2], EK_LINK_D00-13) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


    flux_local_cartesian[1]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0U0) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2], EK_LINK_0D0-13) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


    flux_local_cartesian[2]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_00U) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1], coord[2]-1, EK_LINK_00D-13) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];

    mass_flux_local_cartesian[0] += current_mass_fraction * flux_local_cartesian[0] 
                                    / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
    mass_flux_local_cartesian[1] += current_mass_fraction * flux_local_cartesian[1] 
                                    / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
    mass_flux_local_cartesian[2] += current_mass_fraction * flux_local_cartesian[2] 
                                    / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
  }

  mass_flux[0] = mass_flux_local_cartesian[0];
  mass_flux[1] = mass_flux_local_cartesian[1];
  mass_flux[2] = mass_flux_local_cartesian[2];

  free( fluxes );

  return 0;
}


int ek_lb_print_vtk_density( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_rho_v_pi_gpu ) );
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
           lbpar_gpu.agrid*0.5f, lbpar_gpu.agrid*0.5f, lbpar_gpu.agrid*0.5f,
           lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
           lbpar_gpu.number_of_nodes                                      );

  for( int i = 0; i < lbpar_gpu.number_of_nodes; i++ ) 
  {  
    fprintf( fp, "%e ", host_values[ i ].rho[ 0 ] );
  }
  
  free( host_values );
  fclose( fp );
  
  return 0;
}


int ek_print_vtk_density( int species, char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ){
    return 1;
  }

  float* densities = (float*) malloc( ek_parameters.number_of_nodes *
                                      sizeof( float )                 );
  
  if( ek_parameters.species_index[ species ] != -1 ) 
  {  
    cuda_safe_mem( cudaMemcpy( densities, 
                               ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                               ek_parameters.number_of_nodes * sizeof( float ),
                               cudaMemcpyDeviceToHost )
                 );
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
           ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes,
           species                                                                    );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {  
    fprintf( fp, "%e\n", densities[ i ] / (ek_parameters.agrid*ek_parameters.agrid*ek_parameters.agrid) );
  }
  
  free( densities );
  fclose( fp );
  
  return 0;
}


int ek_node_print_density( int species, int x, int y, int z, double* density ) {

  float* densities = (float*) malloc( ek_parameters.number_of_nodes *
                                      sizeof( float )                 );
  
  if( ek_parameters.species_index[ species ] != -1 ) 
  {  
    cuda_safe_mem( cudaMemcpy( densities, 
                               ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                               ek_parameters.number_of_nodes * sizeof( float ),
                               cudaMemcpyDeviceToHost )
                 );
  }
  else
    return 1;
  
  *density = densities[z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x] / (ek_parameters.agrid*ek_parameters.agrid*ek_parameters.agrid);
  
  free( densities );
  
  return 0;
}


int ek_print_vtk_flux( int species, char* filename ) {

  FILE* fp = fopen( filename, "w" );
  float flux_local_cartesian[3]; //temporary variable for converting fluxes into cartesian coordinates for output

  unsigned int coord[3];

  if( fp == NULL ){
    return 1;
  }

  float* fluxes = (float*) malloc( ek_parameters.number_of_nodes * 13 * sizeof( float ) );
  
  if( ek_parameters.species_index[ species ] != -1 ) 
  {  
    int threads_per_block = 64;
    int blocks_per_grid_y = 4;
    int blocks_per_grid_x =
      ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
      / (threads_per_block * blocks_per_grid_y );
    dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
    
    KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
    KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                ( ek_parameters.species_index[ species ], *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values )    );
              
#ifdef EK_BOUNDARIES
    KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                ( ek_parameters.species_index[ species ], *current_nodes, node_f )                     );
#endif
  
    cuda_safe_mem( cudaMemcpy( fluxes, 
                               ek_parameters.j,
                               ek_parameters.number_of_nodes * 13*sizeof( float ),
                               cudaMemcpyDeviceToHost )
                 );
  }
  else
    return 1;
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
flux_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS flux_%d float 3\n\
LOOKUP_TABLE default\n",
           species,
           ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
           ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes,
           species                                                                    );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {    
    rhoindex_linear2cartesian_host(i, coord);
     
    flux_local_cartesian[0]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U00) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2], EK_LINK_D00-13) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];

    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


    flux_local_cartesian[1]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0U0) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2], EK_LINK_0D0-13) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


    flux_local_cartesian[2]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_00U) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1], coord[2]-1, EK_LINK_00D-13) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
    flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
    flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


    fprintf( fp, "%e %e %e\n",
             flux_local_cartesian[0] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid ),
             flux_local_cartesian[1] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid ),
             flux_local_cartesian[2] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid ) );
  }
  
  free( fluxes );
  fclose( fp );
  
  return 0;
}


int ek_print_vtk_potential( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {  
    return 1;
  }

  float* potential = (float*) malloc( ek_parameters.number_of_nodes * sizeof( cufftReal ) );
  
  cuda_safe_mem( cudaMemcpy( potential, 
                             ek_parameters.charge_potential,
                             ek_parameters.number_of_nodes * sizeof( cufftReal ),
                             cudaMemcpyDeviceToHost )                          
               );
  
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
          ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {  
    fprintf( fp, "%e\n", potential[ i ] );
  }
  
  free( potential );
  fclose( fp );
  
  return 0;
}


int ek_print_vtk_lbforce( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }

  float* lbforce = (float*) malloc( ek_parameters.number_of_nodes * 3 *sizeof( float ) );
  
  cuda_safe_mem( cudaMemcpy( lbforce, 
                             node_f.force,
                             ek_parameters.number_of_nodes * 3 * sizeof( float ),
                             cudaMemcpyDeviceToHost )
               );
  
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
           ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {
    fprintf( fp, "%e %e %e\n", lbforce[ i ] / 
                                 ( 2.0 * lbpar_gpu.rho[0] * powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ),
                               lbforce[ i + ek_parameters.number_of_nodes ] /
                                 ( 2.0 * lbpar_gpu.rho[0] * powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ),
                               lbforce[ i + 2 * ek_parameters.number_of_nodes ] /
                                 ( 2.0 * lbpar_gpu.rho[0] * powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ) );
  }
  
  free( lbforce );
  fclose( fp );
  
  return 0;
}


#ifdef EK_REACTION
int ek_print_vtk_pressure( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }

  float* pressure = (float*) malloc( ek_parameters.number_of_nodes * sizeof( float ) );
  
  cuda_safe_mem( cudaMemcpy( pressure, 
                             ek_parameters.pressure,
                             ek_parameters.number_of_nodes * sizeof( float ),
                             cudaMemcpyDeviceToHost )
               );
  
  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
pressure\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS pressure float 1\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  { 
    fprintf( fp, "%e\n", pressure[ i ] / ek_parameters.agrid );
  }
  
  free( pressure );
  fclose( fp );
  
  return 0;
}


int ek_print_vtk_reaction_tags( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
rection_tags\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS reaction_tags int 1\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes                                              );

  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {  
    fprintf( fp, "%d\n", ek_node_is_catalyst[ i ] );
  }
  
  fclose( fp );
  
  return 0;
}


int ek_print_vtk_mass_flux( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ){
    return 1;
  }

  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
mass_flux\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS mass_flux float 3\n\
LOOKUP_TABLE default\n",
           ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
           ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
           ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
           ek_parameters.number_of_nodes
         );

  float current_mass_fraction;
  float total_mass =   ek_parameters.rho_reactant_reservoir*ek_parameters.mass_reactant
                     + ek_parameters.rho_product0_reservoir*ek_parameters.mass_product0
                     + ek_parameters.rho_product1_reservoir*ek_parameters.mass_product1;

  unsigned int coord[3];
  float flux_local_cartesian[3]; // temporary variable for converting fluxes into cartesian coordinates for output

  float* fluxes = (float*) malloc( ek_parameters.number_of_nodes * 13 * sizeof( float ) );
  float* mass_flux_local_cartesian = (float*) calloc( ek_parameters.number_of_nodes * 3 , sizeof( float ) );
  
  for ( int species_number = 0; species_number < 3; species_number++ )
  {
    int species = ek_parameters.reaction_species[species_number];

    if ( species_number == 0 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_reactant/total_mass;
    else if ( species_number == 1 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_product0/total_mass;
    else if ( species_number == 2 )
      current_mass_fraction = ek_parameters.lb_density*ek_parameters.mass_product1/total_mass;
    else 
      return 1;

    if( ek_parameters.species_index[ species ] != -1 ) 
    {    
      int threads_per_block = 64;
      int blocks_per_grid_y = 4;
      int blocks_per_grid_x =
        ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
        / (threads_per_block * blocks_per_grid_y );
      dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
      
      KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
      KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                  ( ek_parameters.species_index[ species ], *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values ) );
                
#ifdef EK_BOUNDARIES
      KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                  ( ek_parameters.species_index[ species ], *current_nodes, node_f ) );
#endif
    
      cuda_safe_mem( cudaMemcpy( fluxes, 
                                 ek_parameters.j,
                                 ek_parameters.number_of_nodes * 13*sizeof( float ),
                                 cudaMemcpyDeviceToHost )
                   );
    }
    else
      return 1;

    for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
    {      
      rhoindex_linear2cartesian_host(i, coord);
       
      flux_local_cartesian[0]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U00) ];

      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];

      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2], EK_LINK_D00-13) ];

      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];

      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
      flux_local_cartesian[0] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


      flux_local_cartesian[1]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0U0) ];

      flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UU0) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UD0) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

      flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2], EK_LINK_0D0-13) ];

      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2], EK_LINK_DD0-13) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2], EK_LINK_DU0-13) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
      flux_local_cartesian[1] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
      flux_local_cartesian[1] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];


      flux_local_cartesian[2]  = 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_00U) ];

      flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0U) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_U0D) ];
      flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UU) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_0UD) ];

      flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUU) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UUD) ];
      flux_local_cartesian[2] += 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDU) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_getByRhoLinear_host(i, EK_LINK_UDD) ];

      flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1], coord[2]-1, EK_LINK_00D-13) ];

      flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]-1, EK_LINK_D0D-13) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1], coord[2]+1, EK_LINK_D0U-13) ];
      flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]-1, EK_LINK_0DD-13) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0], coord[1]-1, coord[2]+1, EK_LINK_0DU-13) ];

      flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]-1, EK_LINK_DDD-13) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]-1, coord[2]+1, EK_LINK_DDU-13) ];
      flux_local_cartesian[2] += 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]-1, EK_LINK_DUD-13) ];
      flux_local_cartesian[2] -= 0.5*fluxes[ jindex_cartesian2linear_host(coord[0]-1, coord[1]+1, coord[2]+1, EK_LINK_DUU-13) ];

      mass_flux_local_cartesian[3*i + 0] += current_mass_fraction * flux_local_cartesian[0]
                                            / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
      mass_flux_local_cartesian[3*i + 1] += current_mass_fraction * flux_local_cartesian[1]
                                            / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
      mass_flux_local_cartesian[3*i + 2] += current_mass_fraction * flux_local_cartesian[2]
                                            / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );

    }
  }
  
  for( int i = 0; i < ek_parameters.number_of_nodes; i++ ) 
  {
    fprintf( fp, "%e %e %e\n",
             mass_flux_local_cartesian[3*i + 0],
             mass_flux_local_cartesian[3*i + 1],
             mass_flux_local_cartesian[3*i + 2] );
  }

  free( fluxes );
  fclose( fp );
  
  return 0;
}
#endif


void ek_print_parameters() {

  printf( "ek_parameters {\n" );
  
  printf( "  float agrid = %f;\n",                      ek_parameters.agrid );
  printf( "  float time_step = %f;\n",                  ek_parameters.time_step );
  printf( "  float lb_density = %f;\n",                 ek_parameters.lb_density );
  printf( "  unsigned int dim_x = %d;\n",               ek_parameters.dim_x );
  printf( "  unsigned int dim_y = %d;\n",               ek_parameters.dim_y );
  printf( "  unsigned int dim_z = %d;\n",               ek_parameters.dim_z );
  printf( "  unsigned int number_of_nodes = %d;\n",     ek_parameters.number_of_nodes );
  printf( "  float viscosity = %f;\n",                  ek_parameters.viscosity );
  printf( "  float bulk_viscosity = %f;\n",             ek_parameters.bulk_viscosity );
  printf( "  float gamma_odd = %f;\n",                  ek_parameters.gamma_odd );
  printf( "  float gamma_even = %f;\n",                 ek_parameters.gamma_even );
  printf( "  float friction = %f;\n",                   ek_parameters.friction );
  printf( "  float T = %f;\n",                          ek_parameters.T );
  printf( "  float bjerrumlength = %f;\n",              ek_parameters.bjerrumlength );
  printf( "  unsigned int number_of_species = %d;\n",   ek_parameters.number_of_species);
  printf( "  unsigned int accelerated_frame_enabled = %d;\n",             ek_parameters.accelerated_frame_enabled);
  printf( "  float accelerated_frame_boundary_mass_density = %f;\n",      ek_parameters.accelerated_frame_boundary_mass_density);
  printf( "  float accelerated_frame_boundary_mass = %f;\n",              ek_parameters.accelerated_frame_boundary_mass);
  printf( "  float accelerated_frame_fluid_mass = %f;\n",                 ek_parameters.accelerated_frame_fluid_mass);
  printf( "  float ext_acceleration_force[] = {%f, %f, %f};\n",           ek_parameters.ext_acceleration_force[0], 
                                                                          ek_parameters.ext_acceleration_force[1], 
                                                                          ek_parameters.ext_acceleration_force[2] );
  printf( "  int reaction_species[] = {%d, %d, %d};\n", ek_parameters.reaction_species[0], 
                                                        ek_parameters.reaction_species[1], 
                                                        ek_parameters.reaction_species[2] );
  printf( "  float rho_reactant_reservoir = %f;\n",     ek_parameters.rho_reactant_reservoir);
  printf( "  float rho_product0_reservoir = %f;\n",     ek_parameters.rho_product0_reservoir);
  printf( "  float rho_product1_reservoir = %f;\n",     ek_parameters.rho_product1_reservoir);
  printf( "  float reaction_ct_rate = %f;\n",           ek_parameters.reaction_ct_rate); 
  printf( "  float reaction_fraction_0 = %f;\n",        ek_parameters.reaction_fraction_0);
  printf( "  float reaction_fraction_1 = %f;\n",        ek_parameters.reaction_fraction_0);
  printf( "  float reset_mode_0 = %f;\n",               ek_parameters.reset_mode_0);
  printf( "  float* j = %p;\n",                         ek_parameters.j );
  
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
  
  printf("    float rho = %f;\n",                        lbpar_gpu.rho[0] );
  printf("    float mu = %f;\n",                         lbpar_gpu.mu[0] );
  printf("    float viscosity = %f;\n",                  lbpar_gpu.viscosity[0] );
  printf("    float gamma_shear = %f;\n",                lbpar_gpu.gamma_shear[0] );
  printf("    float gamma_bulk = %f;\n",                 lbpar_gpu.gamma_bulk[0] );
  printf("    float gamma_odd = %f;\n",                  lbpar_gpu.gamma_odd[0] );
  printf("    float gamma_even = %f;\n",                 lbpar_gpu.gamma_even[0] );
  printf("    float agrid = %f;\n",                      lbpar_gpu.agrid );
  printf("    float tau = %f;\n",                        lbpar_gpu.tau );
  printf("    float friction = %f;\n",                   lbpar_gpu.friction[0] );
  printf("    float time_step = %f;\n",                  lbpar_gpu.time_step );
  printf("    float lb_coupl_pref = %f;\n",              lbpar_gpu.lb_coupl_pref[0] );
  printf("    float lb_coupl_pref2 = %f;\n",             lbpar_gpu.lb_coupl_pref2[0] );
  printf("    float bulk_viscosity = %f;\n",             lbpar_gpu.bulk_viscosity[0] );
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

  ek_parameters.agrid = agrid;    
  return 0;
}


int ek_set_lb_density( double lb_density ) {  

  ek_parameters.lb_density = lb_density;    
  return 0;
}


int ek_set_bjerrumlength( double bjerrumlength ) {

  ek_parameters.bjerrumlength = bjerrumlength;
  return 0;
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

int ek_set_accelerated_frame( int enabled, double boundary_mass_density, double* ext_acceleration_force ) {

#ifdef EK_BOUNDARIES
  ek_parameters.accelerated_frame_enabled = enabled;
  ek_parameters.accelerated_frame_boundary_mass_density = boundary_mass_density;

  ek_parameters.ext_acceleration_force[0] = ext_acceleration_force[0];
  ek_parameters.ext_acceleration_force[1] = ext_acceleration_force[1];
  ek_parameters.ext_acceleration_force[2] = ext_acceleration_force[2];

  return 0;
#else 
  printf("ERROR: Need boundaries (EK_BOUNDARIES) for the accelerated frame.\n");
  return 1;
#endif

}

int ek_accelerated_frame_print_boundary_velocity( double* accelerated_boundary_velocity ) {

#ifdef EK_BOUNDARIES
  float* temp_boundary_velocity = (float*) malloc( 3 * sizeof( float ) );

  cuda_safe_mem( cudaMemcpy( temp_boundary_velocity, 
                             ek_accelerated_frame_boundary_velocity,
                             3 * sizeof( float ),
                             cudaMemcpyDeviceToHost )
               );

  accelerated_boundary_velocity[0] = static_cast<double>(temp_boundary_velocity[0]);
  accelerated_boundary_velocity[1] = static_cast<double>(temp_boundary_velocity[1]);
  accelerated_boundary_velocity[2] = static_cast<double>(temp_boundary_velocity[2]);

  return 0;
#else 
  printf("ERROR: Need boundaries (EK_BOUNDARIES) for the accelerated frame.\n");
  return 1;
#endif

}


#ifdef EK_REACTION
int ek_set_reaction( int reactant, int product0, int product1, 
                     float rho_reactant_reservoir, float rho_product0_reservoir, float rho_product1_reservoir, 
                     float reaction_ct_rate, float reaction_fraction_0, float reaction_fraction_1,
                     float mass_reactant, float mass_product0, float mass_product1 ) 
{
  if ( ek_parameters.species_index[reactant] == -1 ||
       ek_parameters.species_index[product0] == -1 ||
       ek_parameters.species_index[product1] == -1 ) 
    return 1;

  ek_parameters.reaction_species[0] = reactant;
  ek_parameters.reaction_species[1] = product0;
  ek_parameters.reaction_species[2] = product1;

  ek_parameters.rho_reactant_reservoir = rho_reactant_reservoir;
  ek_parameters.rho_product0_reservoir = rho_product0_reservoir;
  ek_parameters.rho_product1_reservoir = rho_product1_reservoir;

  ek_parameters.reaction_ct_rate = reaction_ct_rate;

  ek_parameters.mass_reactant = mass_reactant;
  ek_parameters.mass_product0 = mass_product0;
  ek_parameters.mass_product1 = mass_product1;

  ek_parameters.reaction_fraction_0 = reaction_fraction_0;
  ek_parameters.reaction_fraction_1 = reaction_fraction_1;  

  return 0;
}

int ek_reset_mode_zero( double reset_mode_0 ) {

  if ( ek_parameters.reaction_species[0] == -1 ||
       ek_parameters.reaction_species[1] == -1 ||
       ek_parameters.reaction_species[2] == -1 ) 
    return 1;
                    
  ek_parameters.reset_mode_0 = reset_mode_0;
  
  return 0;
}

int ek_tag_reaction_nodes( LB_Boundary *boundary, char reaction_type )
{

#ifdef EK_BOUNDARIES
  char *errtxt;
  double pos[3], dist, dist_vec[3];

  for(int z=0; z<int(ek_parameters.dim_z); z++) {
  for(int y=0; y<int(ek_parameters.dim_y); y++) {
  for(int x=0; x<int(ek_parameters.dim_x); x++) {	 

    pos[0] = (x + 0.5)*lbpar_gpu.agrid;
    pos[1] = (y + 0.5)*lbpar_gpu.agrid;
    pos[2] = (z + 0.5)*lbpar_gpu.agrid;

    switch (boundary->type)
    {
      case LB_BOUNDARY_WAL:
        calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.wal, &dist, dist_vec);
        break;
                
      case LB_BOUNDARY_SPH:
        calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.sph, &dist, dist_vec);
        break;
                
      case LB_BOUNDARY_CYL:
        calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.cyl, &dist, dist_vec);
        break;
                
      case LB_BOUNDARY_RHOMBOID:
        calculate_rhomboid_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.rhomboid, &dist, dist_vec);
        break;
                
      case LB_BOUNDARY_POR:
        calculate_pore_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.pore, &dist, dist_vec);
        break;
                
      case LB_BOUNDARY_STOMATOCYTE:
        calculate_stomatocyte_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.stomatocyte, &dist, dist_vec);
        break;

      case LB_BOUNDARY_BOX:
        dist = -1.0;
        break;
                
      case LB_BOUNDARY_HOLLOW_CONE:
        calculate_hollow_cone_dist((Particle*) NULL, pos, (Particle*) NULL, &boundary->c.hollow_cone, &dist, dist_vec);
        break;
                
      default:
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{109 lbboundary type %d not implemented in ek_tag_reaction_nodes()\n", boundary->type);
    }

    if( dist <= 0.0 )
    {
      ek_node_is_catalyst[
                           z * ek_parameters.dim_y * ek_parameters.dim_x +
                           y * ek_parameters.dim_x +
                           x
                         ] = reaction_type;
    }

  }}}

  cuda_safe_mem( cudaMemcpy( ek_parameters.node_is_catalyst, 
                             ek_node_is_catalyst, 
                             ek_parameters.number_of_nodes * sizeof( char ), 
                             cudaMemcpyHostToDevice ) 
               );

  return 0;
#else 
  printf("ERROR: Need boundaries (EK_BOUNDARIES) for the catalytic reaction tagging.\n");
  return 1;
#endif

}
#endif


#endif /* ELECTROKINETICS */

#endif /* CUDA */
