/* 
   Copyright (C) 2010,2011,2012,2014,2015,2016 The ESPResSo project

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
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include "constraint.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "electrokinetics.hpp"
#include "errorhandling.hpp"
#include "fd-electrostatics.hpp"
#include "lb-boundaries.hpp"
#include "lbgpu.hpp"

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

  #define PI_FLOAT 3.14159265358979323846f

  EK_parameters ek_parameters = { -1.0, -1.0, -1.0,
                                     0,    0,    0,
                                     0,
                                  -1.0, -1.0,  0.0,
                                   0.0,  0.0, -1.0,
                                  -1.0,
                                  {0.0,  0.0,  0.0},
                                     0,
                                  { -1,   -1,  -1},
                                  -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0,
                                  0, -1,
                                  true,
                                  true,
#ifdef EK_ELECTROSTATIC_COUPLING
                                  false
#endif                                 
                                };
                                
  static __device__ __constant__ EK_parameters ek_parameters_gpu;
  __device__ float charge_gpu = 0.0f;
  EK_parameters *ek_parameters_gpu_pointer;
  LB_parameters_gpu *ek_lbparameters_gpu;
  CUDA_particle_data *particle_data_gpu;
  float *ek_lb_boundary_force;
  char *ek_node_is_catalyst;
  unsigned int old_number_of_species = 0;
  unsigned int old_number_of_boundaries = 0;

  FdElectrostatics* electrostatics = NULL;

  bool initialized = false;
  
  extern LB_parameters_gpu lbpar_gpu;
  extern LB_node_force_gpu node_f, node_f_buf;
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

__device__ inline void atomicadd (double* address, double value) {
  unsigned long long oldval, newval, readback;
  oldval = __double_as_longlong(*address);
  newval = __double_as_longlong(__longlong_as_double(oldval) + value);
  while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
  {
    oldval = readback;
    newval = __double_as_longlong(__longlong_as_double(oldval) + value);
  }
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

    // Calculate the indices of the neighbors to which
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

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 2) * //TODO CHANGE THIS SCALING AND THE FOLLOWING TWO?
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

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 2) *
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

    pressure_gradient *= powf(ek_parameters_gpu.agrid, 2) *
                         ek_parameters_gpu.time_step *
                         ek_parameters_gpu.time_step;

    pressure_gradient *= ( (   lb_node.boundary[ neighborindex[EK_LINK_00U_pressure] ]
                             + lb_node.boundary[ index ]
                             + lb_node.boundary[ neighborindex[EK_LINK_00D_pressure] ] ) == 0 );

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], pressure_gradient );
  }
}
#endif

__device__ void ek_diffusion_migration_lbforce_nonlinear_stencil(unsigned int index, unsigned int *neighborindex, unsigned int species_index, LB_node_force_gpu node_f) {
  ekfloat flux, force;
  float boltzmannfactor_local, boltzmannfactor_neighbor;

  float agrid_inv = 1.0f / ek_parameters_gpu.agrid;
  float sqrt2agrid_inv = 1.0f / (sqrtf(2.0f) * ek_parameters_gpu.agrid);
  float T_inv =  1.0f / ek_parameters_gpu.T;
  float force_conv = agrid_inv * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step;
  
  boltzmannfactor_local = 
    exp( T_inv *
         ek_parameters_gpu.valency[species_index] *
         ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
       );

  //face in x
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
           ek_parameters_gpu.ext_force[0][species_index] * ek_parameters_gpu.agrid
         )
       );
       
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
           boltzmannfactor_neighbor
         ) * agrid_inv;
         
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U00)],
             flux * ek_parameters_gpu.time_step );

  force  = -1.0f * ek_parameters_gpu.valency[species_index] *
           ( ((cufftReal*)ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
             ((cufftReal*)ek_parameters_gpu.charge_potential)[index]
           ) * agrid_inv;

  force *= force_conv;
           
  atomicadd( &node_f.force[index],
             ek_parameters_gpu.rho[species_index][index] *
             (
               force * 0.5f +
               ek_parameters_gpu.ext_force[0][species_index] * force_conv
             )
           );

  atomicadd( &node_f.force[neighborindex[EK_LINK_U00]],
              ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]] *
              force * 0.5f );
  
  //face in y
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
           ek_parameters_gpu.ext_force[1][species_index] * ek_parameters_gpu.agrid
         )
       );
       
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
           boltzmannfactor_neighbor
         ) * agrid_inv;
         
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
             flux * ek_parameters_gpu.time_step );
            
  force  = -1.0f * ek_parameters_gpu.valency[species_index] *
           ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
             ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
           ) * agrid_inv;

  force *= force_conv;

  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index],
             ek_parameters_gpu.rho[species_index][index] *
             (
               force * 0.5f +
               ek_parameters_gpu.ext_force[1][species_index] * force_conv
             )
           );

  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]],
              ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
              force * 0.5f );

  //face in z
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ( (cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
             ek_parameters_gpu.ext_force[2][species_index] *
             ek_parameters_gpu.agrid
           )
         );
         
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
           boltzmannfactor_neighbor
         ) * agrid_inv;

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_00U)],
             flux * ek_parameters_gpu.time_step );

  force  = -1.0f * ek_parameters_gpu.valency[species_index] *
           ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
             ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
           ) * agrid_inv;

  force *= force_conv;

  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index],
             ek_parameters_gpu.rho[species_index][index] *
             (
               force * 0.5f +
               ek_parameters_gpu.ext_force[2][species_index] * force_conv
             )
           );

  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]],
              ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
              force * 0.5f );

  //edge in z
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UU0]] -
           ( ek_parameters_gpu.ext_force[0][species_index] +
             ek_parameters_gpu.ext_force[1][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
           
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]] *
           boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;
        
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
             flux * ek_parameters_gpu.time_step
           );

  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UD0]] -
           ( ek_parameters_gpu.ext_force[0][species_index] -
             ek_parameters_gpu.ext_force[1][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
  
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]] *
           boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
             flux * ek_parameters_gpu.time_step );

  //edge in y
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0U]] -
           ( ek_parameters_gpu.ext_force[0][species_index] +
             ek_parameters_gpu.ext_force[2][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
  
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
           * boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;
  
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
             flux * ek_parameters_gpu.time_step );

  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0D]] -
           ( ek_parameters_gpu.ext_force[0][species_index] -
             ek_parameters_gpu.ext_force[2][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
  
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]] *
           boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
             flux * ek_parameters_gpu.time_step );

  //edge in x
  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UU]] -
           ( ek_parameters_gpu.ext_force[1][species_index] +
             ek_parameters_gpu.ext_force[2][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
  
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] * boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]] *
           boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
             flux * ek_parameters_gpu.time_step );

  boltzmannfactor_neighbor =
    exp( T_inv *
         ( ek_parameters_gpu.valency[species_index] *
           ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UD]] -
           ( ek_parameters_gpu.ext_force[1][species_index] -
             ek_parameters_gpu.ext_force[2][species_index]
           ) * ek_parameters_gpu.agrid
         )
       );
  
  flux = ( ek_parameters_gpu.d[species_index] * agrid_inv ) *
         ( 1.0f / boltzmannfactor_local +
           1.0f / boltzmannfactor_neighbor
         ) * 0.5f *
         ( ek_parameters_gpu.rho[species_index][index] *
           boltzmannfactor_local -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]] *
           boltzmannfactor_neighbor
         ) * sqrt2agrid_inv;
       
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
             flux * ek_parameters_gpu.time_step );
}

__device__ void ek_diffusion_migration_lbforce_linkcentered_stencil(unsigned int index, unsigned int *neighborindex, unsigned int species_index, LB_node_force_gpu node_f, LB_nodes_gpu lb_node) {
  ekfloat flux, force;

  float agrid_inv = 1.0f / ek_parameters_gpu.agrid;
  float sqrt2agrid_inv = 1.0f / (sqrtf(2.0f) * ek_parameters_gpu.agrid);
  float sqrt2_inv = 1.0f / sqrt(2.0f);
  float twoT_inv =  1.0f / (2.0f * ek_parameters_gpu.T);
  float D_inv = 1.0f / ek_parameters_gpu.D[species_index];
  float force_conv = agrid_inv * ek_parameters_gpu.time_step * ek_parameters_gpu.time_step;

  //face in x
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]]
         ) * agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]]
            ) * agrid_inv
            +
            ek_parameters_gpu.ext_force[0][species_index]
          );
         
  flux += force * 
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U00]]);
         
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U00)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[index], force * 0.5f);
    atomicadd( &node_f.force[neighborindex[EK_LINK_U00]], force * 0.5f);
  }
  else
  {
    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*)ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]] -
               ((cufftReal*)ek_parameters_gpu.charge_potential)[index]
             ) * agrid_inv;

    force *= force_conv;
             
    atomicadd( &node_f.force[index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force * 0.5f +
                 ek_parameters_gpu.ext_force[0][species_index] * force_conv
               )
             );
  }

  //face in y
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]]
         ) * agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]]
            ) * agrid_inv
            +
            ek_parameters_gpu.ext_force[1][species_index]
          );             
         
  flux += force * 
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;
         
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0U0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
             flux * ek_parameters_gpu.time_step );
            
  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]], force * 0.5f);
  }
  else
  {
    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) * agrid_inv;

    force *= force_conv;

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force * 0.5f +
                 ek_parameters_gpu.ext_force[1][species_index] * force_conv
               )
             );

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]] *
                force * 0.5f );
  }

  //face in z
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]]
         ) * agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]]
            ) * agrid_inv
            +
            ek_parameters_gpu.ext_force[2][species_index]
          );             
         
  flux += force * 
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_00U]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_00U)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]], force * 0.5f);
  }
  else
  {
    force  = -1.0f * ek_parameters_gpu.valency[species_index] *
             ( ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]] -
               ((cufftReal*) ek_parameters_gpu.charge_potential)[index]
             ) * agrid_inv;

    force *= force_conv;

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index],
               ek_parameters_gpu.rho[species_index][index] *
               (
                 force * 0.5f +
                 ek_parameters_gpu.ext_force[2][species_index] * force_conv
               )
             );

    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]],
                ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]] *
                force * 0.5f );
  }

  //edge in z
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UU0]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[0][species_index] +
              ek_parameters_gpu.ext_force[1][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;
         
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_UU0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
             flux * ek_parameters_gpu.time_step
           );
  
  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[index], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[neighborindex[EK_LINK_UU0]], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UU0]], force * 0.5f);
  }

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UD0]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[0][species_index] -
              ek_parameters_gpu.ext_force[1][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_UD0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[index], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], -force * 0.5f);
    atomicadd( &node_f.force[neighborindex[EK_LINK_UD0]], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UD0]], -force * 0.5f);
  }

  //edge in y
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0U]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[0][species_index] +
              ek_parameters_gpu.ext_force[2][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U0U]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[index], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[neighborindex[EK_LINK_U0U]], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0U]], force * 0.5f);
  }

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0D]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[0][species_index] -
              ek_parameters_gpu.ext_force[2][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;
  
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U0D]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[index], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], -force * 0.5f);
    atomicadd( &node_f.force[neighborindex[EK_LINK_U0D]], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0D]], -force * 0.5f);
  }

  //edge in x
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UU]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[1][species_index] +
              ek_parameters_gpu.ext_force[2][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0UU]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force * 0.5f);
  }

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]]
         ) * sqrt2agrid_inv;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UD]]
            ) * sqrt2agrid_inv +
            ( ek_parameters_gpu.ext_force[1][species_index] -
              ek_parameters_gpu.ext_force[2][species_index]
            ) * sqrt2_inv
          );
         
  flux += force *
          ( ek_parameters_gpu.rho[species_index][index] +
            ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]]
          ) * twoT_inv;

  flux *= ek_parameters_gpu.d[species_index] * agrid_inv;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0UD]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
             flux * ek_parameters_gpu.time_step );

  if(ek_parameters_gpu.fluidcoupling_ideal_contribution)
  {
    force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid * D_inv;
    force *= force_conv;

    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], -force * 0.5f);
    atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], force * 0.5f);
    atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], -force * 0.5f);
  }
}

__device__ void ek_diffusion_migration_lbforce_nodecentered_stencil(unsigned int index, unsigned int *neighborindex, unsigned int species_index, LB_node_force_gpu node_f, LB_nodes_gpu lb_node) {
  ekfloat flux, force;

  //face in x
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]]
         ) / ek_parameters_gpu.agrid;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U00]]
            ) / ek_parameters_gpu.agrid
            +
            ek_parameters_gpu.ext_force[0][species_index]
          );
         
  flux += force * 
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U00]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U00]]);
         
  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U00)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[index], force / 2.0f);
  atomicadd( &node_f.force[neighborindex[EK_LINK_U00]], force / 2.0f);

  //face in y
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]]
         ) / ek_parameters_gpu.agrid;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0U0]]
            ) / ek_parameters_gpu.agrid
            +
            ek_parameters_gpu.ext_force[1][species_index]
          );             
         
  flux += force * 
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0U0]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;
         
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0U0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
             flux * ek_parameters_gpu.time_step );
            
  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0U0]], force / 2.0f);

  //face in z
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]]
         ) / ek_parameters_gpu.agrid;

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_00U]]
            ) / ek_parameters_gpu.agrid
            +
            ek_parameters_gpu.ext_force[2][species_index]
          );             
         
  flux += force * 
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_00U]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_00U]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_00U)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_00U]], force / 2.0f);

  //edge in z
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UU0]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[0][species_index] +
              ek_parameters_gpu.ext_force[1][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UU0]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;
         
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_UU0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
             flux * ek_parameters_gpu.time_step
           );
  
  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[index], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[neighborindex[EK_LINK_UU0]], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UU0]], force / 2.0f);

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_UD0]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[0][species_index] -
              ek_parameters_gpu.ext_force[1][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_UD0]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_UD0]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[index], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], -force / 2.0f);
  atomicadd( &node_f.force[neighborindex[EK_LINK_UD0]], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_UD0]], -force / 2.0f);

  //edge in y
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0U]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[0][species_index] +
              ek_parameters_gpu.ext_force[2][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0U]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U0U]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[index], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[neighborindex[EK_LINK_U0U]], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0U]], force / 2.0f);

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_U0D]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[0][species_index] -
              ek_parameters_gpu.ext_force[2][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_U0D]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;
  
  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_U0D]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[index], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], -force / 2.0f);
  atomicadd( &node_f.force[neighborindex[EK_LINK_U0D]], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_U0D]], -force / 2.0f);

  //edge in x
  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UU]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[1][species_index] +
              ek_parameters_gpu.ext_force[2][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UU]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0UU]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UU]], force / 2.0f);

  flux = ( ek_parameters_gpu.rho[species_index][index] -
           ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]]
         ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid);

  force = ( ek_parameters_gpu.valency[species_index] *
            ( ((cufftReal*) ek_parameters_gpu.charge_potential)[index] -
              ((cufftReal*) ek_parameters_gpu.charge_potential)[neighborindex[EK_LINK_0UD]]
            ) / (sqrtf(2.0f) * ek_parameters_gpu.agrid) +
            ( ek_parameters_gpu.ext_force[1][species_index] -
              ek_parameters_gpu.ext_force[2][species_index]
            ) / sqrtf(2.0f)
          );
         
  flux += force *
          ( (force >= 0.0f) * ek_parameters_gpu.rho[species_index][index] +
            (force <  0.0f) * ek_parameters_gpu.rho[species_index][neighborindex[EK_LINK_0UD]]
          ) / ek_parameters_gpu.T;

  flux *= ek_parameters_gpu.d[species_index] / ek_parameters_gpu.agrid;

  flux *= !(lb_node.boundary[index] || lb_node.boundary[neighborindex[EK_LINK_0UD]]);

  atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
             flux * ek_parameters_gpu.time_step );

  force = flux * ek_parameters_gpu.T * ek_parameters_gpu.agrid / ek_parameters_gpu.D[species_index];

  force *= powf(ek_parameters_gpu.agrid, -1) *
           ek_parameters_gpu.time_step *
           ek_parameters_gpu.time_step;

  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + index], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + index], -force / 2.0f);
  atomicadd( &node_f.force[ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], force / 2.0f);
  atomicadd( &node_f.force[2*ek_parameters_gpu.number_of_nodes + neighborindex[EK_LINK_0UD]], -force / 2.0f);
}

__device__ void ek_add_advection_to_flux(unsigned int index, unsigned int *neighborindex, unsigned int *coord, unsigned int species_index, LB_node_force_gpu node_f, LB_nodes_gpu lb_node, LB_parameters_gpu *ek_lbparameters_gpu) {
    float dx[3];
    int di[3];
    int node;

    ek_displacement( dx, lb_node, index, ek_lbparameters_gpu );

    di[0] = 1 - signbit(dx[0]);
    di[1] = 1 - signbit(dx[1]);
    di[2] = 1 - signbit(dx[2]);

    dx[0] = fabs(dx[0]);
    dx[1] = fabs(dx[1]);
    dx[2] = fabs(dx[2]);

    int target_node[3];
    int target_node_index;
    int not_boundary;

    //face in x
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        coord[2]
      );
    
    target_node[0] = (coord[0] + 2*di[0]-1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x;
    target_node[1] = coord[1];
    target_node[2] = coord[2];
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

    //if(!not_boundary)
    //  printf("[%d,%d,%d]=%d\n", coord[0], coord[1], coord[2], not_boundary); //TODO delete
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_U00 )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * (1.0f - dx[1]) * (1.0f - dx[2]) * not_boundary );
    
    //face in y
    node =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] + di[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
    
    target_node[0] = coord[0];
    target_node[1] = (coord[1] + 2*di[1]-1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y;
    target_node[2] = coord[2];
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_0U0 )],
              (2 * di[1] - 1) * ek_parameters_gpu.rho[species_index][index] *
              (1.0f - dx[0]) * dx[1] * (1.0f - dx[2]) * not_boundary );

    //face in z
    node =
      rhoindex_cartesian2linear(
        coord[0],
        coord[1],
        (coord[2] + di[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
    
    target_node[0] = coord[0];
    target_node[1] = coord[1];
    target_node[2] = (coord[2] + 2*di[2]-1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z;
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
    
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_00U )],
               (2 * di[2] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0f - dx[0]) * (1.0f - dx[1]) * dx[2] * not_boundary );
    
    //edge in x
    node =
      rhoindex_cartesian2linear(
        coord[0],
        (coord[1] + di[1] - 1                   + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] + (1 - di[1]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
    
    target_node[0] = coord[0];
    target_node[1] = (coord[1] + 2*di[1]-1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y;
    target_node[2] = (coord[2] + 2*di[2]-1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z;
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
           
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_0UU + (di[1] + di[2] == 1) )],
               (2 * di[1] - 1) * ek_parameters_gpu.rho[species_index][index] *
               (1.0f - dx[0]) * dx[1] * dx[2] * not_boundary );
    
    //edge in y
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        coord[1],
        (coord[2] + (1 - di[0]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
    
    target_node[0] = (coord[0] + 2*di[0]-1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x;
    target_node[1] = coord[1];
    target_node[2] = (coord[2] + 2*di[2]-1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z;
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_U0U + (di[0] + di[2] == 1) )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * (1.0f - dx[1]) * dx[2] * not_boundary );
    
    //edge in z
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + (1 - di[0]) * (2*di[1] - 1) + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        coord[2]
      );
    
    target_node[0] = (coord[0] + 2*di[0]-1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x;
    target_node[1] = (coord[1] + 2*di[1]-1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y;
    target_node[2] = coord[2];
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
      
    atomicadd( &ek_parameters_gpu.j[jindex_getByRhoLinear( node, EK_LINK_UU0 + (di[0] + di[1] == 1) )],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * dx[1] * (1.0f - dx[2]) * not_boundary );
    
    //corner
    node =
      rhoindex_cartesian2linear(
        (coord[0] + di[0] - 1                   + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x,
        (coord[1] + (1 - di[0]) * (2*di[1] - 1) + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y,
        (coord[2] + (1 - di[0]) * (2*di[2] - 1) + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z
      );
    
    target_node[0] = (coord[0] + 2*di[0]-1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x;
    target_node[1] = (coord[1] + 2*di[1]-1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y;
    target_node[2] = (coord[2] + 2*di[2]-1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z;
    target_node_index = rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
    not_boundary = (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;
      
    atomicadd( &ek_parameters_gpu.j[
                jindex_getByRhoLinear( node, (1 - di[0]) *
                                             (EK_LINK_UUU + 2*di[1] + di[2]) +
                                             di[0] * (EK_LINK_UDD - 2*di[1] - di[2])
                                     ) ],
               (2 * di[0] - 1) * ek_parameters_gpu.rho[species_index][index] *
               dx[0] * dx[1] * dx[2] * not_boundary );
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
    if(ek_parameters_gpu.stencil == 1) //nonlinear
      ek_diffusion_migration_lbforce_nonlinear_stencil(index, neighborindex, species_index, node_f);
    else if(ek_parameters_gpu.stencil == 0) //link centered
      ek_diffusion_migration_lbforce_linkcentered_stencil(index, neighborindex, species_index, node_f, lb_node);
    else if(ek_parameters_gpu.stencil == 2) //node centered
      ek_diffusion_migration_lbforce_nodecentered_stencil(index, neighborindex, species_index, node_f, lb_node);

    /* advective contribution to flux */
    if(ek_parameters_gpu.advection)
      ek_add_advection_to_flux(index, neighborindex, coord, species_index, node_f, lb_node, ek_lbparameters_gpu);
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
  unsigned int coord[3];

  rhoindex_linear2cartesian(index, coord);

  if(index < ek_parameters_gpu.number_of_nodes) 
  {  
    for(int i = 0; i < ek_parameters_gpu.number_of_species; i++) 
    {
//      if(coord[0] == ek_parameters_gpu.dim_x/2 && coord[1] == ek_parameters_gpu.dim_y/2 && coord[2] == ek_parameters_gpu.dim_z/2)
        ek_parameters_gpu.rho[ i ][ index ] = ek_parameters_gpu.density[ i ] *
                                              ek_parameters_gpu.agrid *
                                              ek_parameters_gpu.agrid *
                                              ek_parameters_gpu.agrid;
//      else
//        ek_parameters_gpu.rho[ i ][ index ] = 0.0f;
    }
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
  int lowernode[3];
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
    
    //printf("particle %d (%d):\n  charge %f\n  pos %f %f %f\n  lowernode %d %d %d\n  cellpos %f %f %f\n\n", index, ek_lbparameters_gpu->number_of_particles, particle_data[index].q, particle_data[index].p[0], particle_data[index].p[1], particle_data[index].p[2], lowernode[0], lowernode[1], lowernode[2], cellpos[0], cellpos[1], cellpos[2]); //TODO delete
  }
}
#ifdef EK_ELECTROSTATIC_COUPLING
__global__ void ek_spread_particle_force( CUDA_particle_data * particle_data,
                                          float *particle_forces,
                                          LB_parameters_gpu * ek_lbparameters_gpu ) {

  unsigned int index = ek_getThreadIndex();
  int lowernode[3];
  float cellpos[3];
  float gridpos;

  if( index < ek_lbparameters_gpu->number_of_particles ) 
    {  
      gridpos      = particle_data[ index ].p[0] / ek_parameters_gpu.agrid - 0.5f;
      lowernode[0] = (int) floorf( gridpos );
      cellpos[0]   = gridpos - (float)(lowernode[0]);
  
      gridpos      = particle_data[ index ].p[1] / ek_parameters_gpu.agrid - 0.5f;
      lowernode[1] = (int) floorf( gridpos );
      cellpos[1]   = gridpos - (float)(lowernode[1]);
  
      gridpos      = particle_data[ index ].p[2] / ek_parameters_gpu.agrid - 0.5f;
      lowernode[2] = (int) floorf( gridpos );
      cellpos[2]   = gridpos - (float)(lowernode[2]);

      lowernode[0] = (lowernode[0] + ek_lbparameters_gpu->dim_x) % ek_lbparameters_gpu->dim_x;
      lowernode[1] = (lowernode[1] + ek_lbparameters_gpu->dim_y) % ek_lbparameters_gpu->dim_y;
      lowernode[2] = (lowernode[2] + ek_lbparameters_gpu->dim_z) % ek_lbparameters_gpu->dim_z;

      float efield[3] = { 0., 0., 0. };
#pragma unroll 3
      for(unsigned int dim = 0; dim < 3; ++dim) {
        // 0 0 0
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( lowernode[0],
                                                                                  lowernode[1],
                                                                                    lowernode[2]) + dim]
          * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

        // 0 0 1
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( lowernode[0],
                                                                                  lowernode[1],
                                                                                  (lowernode[2] + 1 ) % ek_lbparameters_gpu->dim_z ) + dim]
          * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * cellpos[2];

        // 0 1 0
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( lowernode[0],
                                                                                  (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                                                                                  lowernode[2]  ) + dim]
          * ( 1 - cellpos[0] ) * cellpos[1] * ( 1 - cellpos[2] );

        // 0 1 1
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( lowernode[0],
                                                                                  (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                                                                                  (lowernode[2] + 1 ) % ek_lbparameters_gpu->dim_z ) + dim]
          * ( 1 - cellpos[0] ) * cellpos[1] * cellpos[2];

        // 1 0 0
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                                                                                  lowernode[1],
                                                                                  lowernode[2]  ) + dim]
          * cellpos[0] * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

        // 1 0 1
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                                                                                  lowernode[1],
                                                                                  (lowernode[2] + 1 ) % ek_lbparameters_gpu->dim_z ) + dim]
          * cellpos[0] * ( 1 - cellpos[1] ) * cellpos[2];

        // 1 1 0
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                                                                                  (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                                                                                  lowernode[2]  ) + dim]
          * cellpos[0] * cellpos[1] * ( 1 - cellpos[2] );

        // 1 1 1
        efield[dim] += ek_parameters_gpu.electric_field[3*rhoindex_cartesian2linear( (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                                                                                  (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                                                                                  (lowernode[2] + 1 ) % ek_lbparameters_gpu->dim_z ) + dim]
          * cellpos[0] * cellpos[1] * cellpos[2];
      }
      particle_forces[3*index + 0] += particle_data[ index ].q * efield[0];
      particle_forces[3*index + 1] += particle_data[ index ].q * efield[1];
      particle_forces[3*index + 2] += particle_data[ index ].q * efield[2];
    }  
}

__global__ void ek_calc_electric_field(const float *potential) {
  unsigned int coord[3];
  const unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) {
    rhoindex_linear2cartesian(index, coord);
    const float agrid_inv = 1.0f / ek_parameters_gpu.agrid;

    ek_parameters_gpu.electric_field[3*index + 0] = -0.5f * agrid_inv *
      (
         potential[rhoindex_cartesian2linear((coord[0] + 1) % ek_parameters_gpu.dim_x, coord[1], coord[2])]
       - potential[rhoindex_cartesian2linear((coord[0] - 1 + ek_parameters_gpu.dim_x) % ek_parameters_gpu.dim_x, coord[1], coord[2])]
       );
    ek_parameters_gpu.electric_field[3*index + 1] = -0.5f * agrid_inv *
      (
         potential[rhoindex_cartesian2linear(coord[0], (coord[1] + 1) % ek_parameters_gpu.dim_y, coord[2])]
       - potential[rhoindex_cartesian2linear(coord[0], (coord[1] - 1 + ek_parameters_gpu.dim_y) % ek_parameters_gpu.dim_y, coord[2])]
       );
    ek_parameters_gpu.electric_field[3*index + 2] = -0.5f * agrid_inv *
      (
         potential[rhoindex_cartesian2linear(coord[0], coord[1], (coord[2] + 1) % ek_parameters_gpu.dim_z)]
       - potential[rhoindex_cartesian2linear(coord[0], coord[1], (coord[2] - 1 + ek_parameters_gpu.dim_z) % ek_parameters_gpu.dim_z)]
      );
  }
}
#endif

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


__global__ void ek_calculate_system_charge() {

  unsigned int index = ek_getThreadIndex();

  if( index < ek_parameters_gpu.number_of_nodes ) 
  {
    for(int i = 0; i < ek_parameters_gpu.number_of_species; i++)
    {
      atomicadd(&charge_gpu, ek_parameters_gpu.rho[i][index] * ek_parameters_gpu.valency[i]);
    }
  }
}


//TODO delete ?? (it has the previous step setting now)
//This is not compatible with external LB forces!
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

#ifdef EK_ELECTROSTATIC_COUPLING
void ek_calculate_electrostatic_coupling() {
  int blocks_per_grid_x;
  int blocks_per_grid_y = 4;
  int threads_per_block = 64;
  dim3 dim_grid;

if((!ek_parameters.es_coupling) || (!initialized))
    return;

  blocks_per_grid_x =
    ( lbpar_gpu.number_of_particles + threads_per_block * blocks_per_grid_y - 1 ) /
    ( threads_per_block * blocks_per_grid_y );
  dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );

  KERNELCALL( ek_spread_particle_force, dim_grid, threads_per_block, 
		(gpu_get_particle_pointer(), gpu_get_particle_force_pointer(), ek_lbparameters_gpu ));
}
#endif

void ek_integrate_electrostatics() {

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 ) /
    ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( ek_gather_species_charge_density, dim_grid, threads_per_block, () );

#ifdef EK_ELECTROSTATIC_COUPLING
    if(ek_parameters.es_coupling) {
      cuda_safe_mem( cudaMemcpy(ek_parameters.charge_potential_buffer, ek_parameters.charge_potential, ek_parameters.number_of_nodes * sizeof(cufftReal), cudaMemcpyDeviceToDevice));
      electrostatics->calculatePotential((cufftComplex *)ek_parameters.charge_potential_buffer);
      KERNELCALL( ek_calc_electric_field, dim_grid, threads_per_block, (ek_parameters.charge_potential_buffer));
    }
#endif

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
  
  electrostatics->calculatePotential();
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
     (in ek_calculate_quantities / ek_displacement), which is copied in this routine */

  //KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, ( node_f ) );



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


  /* Integrate diffusion-advection */
  
  for( int i = 0; i < ek_parameters.number_of_species; i++ )
  {
  
    KERNELCALL( ek_clear_fluxes, dim_grid, threads_per_block, () );
    KERNELCALL( ek_calculate_quantities, dim_grid, threads_per_block,
                ( i, *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values ) );
              
#ifdef EK_BOUNDARIES
    if( ek_parameters.stencil == 1)
    {
      KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                  ( i, *current_nodes, node_f ) );
    }
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

  
  //TODO delete - needed for printfs
  //cudaDeviceSynchronize();
}


#ifdef EK_BOUNDARIES
void ek_init_species_density_wallcharge( ekfloat* wallcharge_species_density,
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
                               ek_parameters.number_of_nodes * sizeof( ekfloat ),
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
                               ek_parameters.number_of_nodes * sizeof( ekfloat )                        ) );
    
    ek_parameters.density[      ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.D[            ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.valency[      ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[0][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[1][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.ext_force[2][ ek_parameters.species_index[ species ] ] = 0.0;
    ek_parameters.d[            ek_parameters.species_index[ species ] ] =
    ek_parameters.D[            ek_parameters.species_index[ species ] ] / ( 1.0 + 2.0 * sqrt( 2.0 ) );
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

    lbpar_gpu.is_TRT = true;

    lb_reinit_parameters_gpu();
    
    lb_init_gpu();

    if (ek_parameters.lb_force[0] != 0 || ek_parameters.lb_force[1] != 0 || ek_parameters.lb_force[2] != 0)
    {
      lbpar_gpu.external_force = 1;
      lbpar_gpu.ext_force[0] = ek_parameters.lb_force[0];
      lbpar_gpu.ext_force[1] = ek_parameters.lb_force[1];
      lbpar_gpu.ext_force[2] = ek_parameters.lb_force[2];
      lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    }
    else
    {
      lbpar_gpu.external_force = 0;
      lbpar_gpu.ext_force[0] = 0;
      lbpar_gpu.ext_force[1] = 0;
      lbpar_gpu.ext_force[2] = 0;
    }

    ek_parameters.dim_x = lbpar_gpu.dim_x;
    ek_parameters.dim_y = lbpar_gpu.dim_y;
    ek_parameters.dim_z = lbpar_gpu.dim_z;
    ek_parameters.time_step = lbpar_gpu.time_step;
    ek_parameters.number_of_nodes = ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z;

    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.j,
                             ek_parameters.number_of_nodes * 13 * sizeof( ekfloat ) ) );
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

#ifdef EK_ELECTROSTATIC_COUPLING
    if(ek_parameters.es_coupling) {
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.charge_potential_buffer,
                               ek_parameters.number_of_nodes * sizeof( cufftComplex ) ) );
    cuda_safe_mem( cudaMalloc( (void**) &ek_parameters.electric_field,
                               ek_parameters.number_of_nodes * 3 * sizeof( float ) ) );
  }
#endif

    lb_get_device_values_pointer( &ek_lb_device_values );
    
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
   
    //initialize electrostatics
    if(electrostatics != NULL)
      delete electrostatics;

    FdElectrostatics::InputParameters es_parameters = {ek_parameters.bjerrumlength, ek_parameters.T, int(ek_parameters.dim_x), int(ek_parameters.dim_y), int(ek_parameters.dim_z), ek_parameters.agrid};
    try {
      electrostatics = new FdElectrostatics(es_parameters, stream[0]);
    }
    catch(std::string e) {
      std::cout << "Error in initialization of electrokinetics electrostatics solver: " << e << std::endl;
      return 1;
    }

    ek_parameters.charge_potential = electrostatics->getGrid().grid;
    cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );

    //clear initial LB force and finish up
    blocks_per_grid_x =
      ( ek_parameters.dim_z * ek_parameters.dim_y * (ek_parameters.dim_x ) +
        threads_per_block * blocks_per_grid_y - 1
      ) / ( threads_per_block * blocks_per_grid_y );
    dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
    KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, ( node_f ) );

    initialized = true;
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
      lb_init_boundaries();
      lb_get_boundary_force_pointer( &ek_lb_boundary_force );
      
      cuda_safe_mem( cudaMemcpyToSymbol( ek_parameters_gpu, &ek_parameters, sizeof( EK_parameters ) ) );
#else
      blocks_per_grid_x =
        ( ek_parameters.number_of_nodes + threads_per_block * blocks_per_grid_y - 1 )
        / (threads_per_block * blocks_per_grid_y );
      dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );

      KERNELCALL( ek_init_species_density_homogeneous, dim_grid, threads_per_block, () );
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
  unsigned int* bound_array = (unsigned int*) Utils::malloc( lbpar_gpu.number_of_nodes*sizeof(unsigned int) );

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
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) Utils::malloc( lbpar_gpu.number_of_nodes *
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
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) Utils::malloc( lbpar_gpu.number_of_nodes *
                                                        sizeof( LB_rho_v_pi_gpu ) );
  lb_get_values_GPU( host_values );
  
  int i = z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;
  
  velocity[0] = host_values[i].v[0];
  velocity[1] = host_values[i].v[1];
  velocity[2] = host_values[i].v[2];
  
  free(host_values);
  
  return 0;
}


int ek_lb_print_vtk_density( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }
  
  LB_rho_v_pi_gpu *host_values = (LB_rho_v_pi_gpu*) Utils::malloc( lbpar_gpu.number_of_nodes *
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

  ekfloat* densities = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes *
                                      sizeof( ekfloat )                 );
  
  if( ek_parameters.species_index[ species ] != -1 ) 
  {  
    cuda_safe_mem( cudaMemcpy( densities, 
                               ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                               ek_parameters.number_of_nodes * sizeof( ekfloat ),
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

  ekfloat* densities = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes *
                                      sizeof( ekfloat )                 );
  
  if( ek_parameters.species_index[ species ] != -1 ) 
  {  
    cuda_safe_mem( cudaMemcpy( densities, 
                               ek_parameters.rho[ ek_parameters.species_index[ species ] ],
                               ek_parameters.number_of_nodes * sizeof( ekfloat ),
                               cudaMemcpyDeviceToHost )
                 );
  }
  else
    return 1;
  
  *density = densities[z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x] / (ek_parameters.agrid*ek_parameters.agrid*ek_parameters.agrid);
  
  free( densities );
  
  return 0;
}


int ek_node_print_flux( int species, int x, int y, int z, double* flux ) {

  ekfloat flux_local_cartesian[3]; //temporary variable for converting fluxes into cartesian coordinates for output
  unsigned int coord[3];

  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  ekfloat* fluxes = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes * 13 * sizeof( ekfloat ) );
  
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
    reset_LB_forces_GPU(false);
              
#ifdef EK_BOUNDARIES
    KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                ( ek_parameters.species_index[ species ], *current_nodes, node_f )                     );
#endif
  
    cuda_safe_mem( cudaMemcpy( fluxes, 
                               ek_parameters.j,
                               ek_parameters.number_of_nodes * 13*sizeof( ekfloat ),
                               cudaMemcpyDeviceToHost )
                 );
  }
  else
    return 1;
  
  int i = rhoindex_cartesian2linear_host(coord[0], coord[1], coord[2]);
   
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

  flux[0] = flux_local_cartesian[0] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
  flux[1] = flux_local_cartesian[1] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
  flux[2] = flux_local_cartesian[2] / ( ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid );
  
  free( fluxes );
  
  return 0;
}


int ek_node_set_density(int species, int x, int y, int z, double density) {
  if(ek_parameters.species_index[species] != -1) 
  {  
    int index = z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;
    ekfloat num_particles = density * ek_parameters.agrid*ek_parameters.agrid*ek_parameters.agrid;

    cuda_safe_mem( cudaMemcpy( &ek_parameters.rho[ek_parameters.species_index[species]][index],
                               &num_particles,
                               sizeof(ekfloat),
                               cudaMemcpyHostToDevice )
                 );
  }
  else
    return 1;
  
  return 0;
}


int ek_print_vtk_flux( int species, char* filename ) {

  FILE* fp = fopen( filename, "w" );
  ekfloat flux_local_cartesian[3]; //temporary variable for converting fluxes into cartesian coordinates for output

  unsigned int coord[3];

  if( fp == NULL ){
    return 1;
  }

  ekfloat* fluxes = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes * 13 * sizeof( ekfloat ) );
  
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
    reset_LB_forces_GPU(false);
              
#ifdef EK_BOUNDARIES
    KERNELCALL( ek_apply_boundaries, dim_grid, threads_per_block,
                ( ek_parameters.species_index[ species ], *current_nodes, node_f )                     );
#endif
  
    cuda_safe_mem( cudaMemcpy( fluxes, 
                               ek_parameters.j,
                               ek_parameters.number_of_nodes * 13*sizeof( ekfloat ),
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

  float* potential = (float*) Utils::malloc( ek_parameters.number_of_nodes * sizeof( cufftReal ) );
  
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

#ifdef EK_ELECTROSTATIC_COUPLING
int ek_print_vtk_particle_potential( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {  
    return 1;
  }

  float* potential = (float*) Utils::malloc( ek_parameters.number_of_nodes * sizeof( cufftReal ) );
  
  cuda_safe_mem( cudaMemcpy( potential, 
                             ek_parameters.charge_potential_buffer,
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
#endif

int ek_print_vtk_lbforce( char* filename ) {
#ifndef EK_DEBUG
  return 1;
#else

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }

  lbForceFloat* lbforce = (lbForceFloat*) Utils::malloc( ek_parameters.number_of_nodes * 3 *sizeof( lbForceFloat ) );
  
  cuda_safe_mem( cudaMemcpy( lbforce, 
                             node_f.force_buf,
                             ek_parameters.number_of_nodes * 3 * sizeof( lbForceFloat ),
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
                                 ( powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ),
                               lbforce[ i + ek_parameters.number_of_nodes ] /
                                 ( powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ),
                               lbforce[ i + 2 * ek_parameters.number_of_nodes ] /
                                 ( powf( ek_parameters.time_step , 2.0 ) * powf( ek_parameters.agrid, 4.0 ) ) );
  }
  
  free( lbforce );
  fclose( fp );
  
  return 0;
#endif
}


#ifdef EK_REACTION
int ek_print_vtk_pressure( char* filename ) {

  FILE* fp = fopen( filename, "w" );

  if( fp == NULL ) 
  {
    return 1;
  }

  float* pressure = (float*) Utils::malloc( ek_parameters.number_of_nodes * sizeof( float ) );
  
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
  printf( "  float lb_force[] = {%f, %f, %f};\n",       ek_parameters.lb_force[0], 
                                                        ek_parameters.lb_force[1], 
                                                        ek_parameters.lb_force[2] );
  printf( "  unsigned int number_of_species = %d;\n",   ek_parameters.number_of_species);
  printf( "  int reaction_species[] = {%d, %d, %d};\n", ek_parameters.reaction_species[0], 
                                                        ek_parameters.reaction_species[1], 
                                                        ek_parameters.reaction_species[2] );
  printf( "  float rho_reactant_reservoir = %f;\n",     ek_parameters.rho_reactant_reservoir);
  printf( "  float rho_product0_reservoir = %f;\n",     ek_parameters.rho_product0_reservoir);
  printf( "  float rho_product1_reservoir = %f;\n",     ek_parameters.rho_product1_reservoir);
  printf( "  float reaction_ct_rate = %f;\n",           ek_parameters.reaction_ct_rate); 
  printf( "  float reaction_fraction_0 = %f;\n",        ek_parameters.reaction_fraction_0);
  printf( "  float reaction_fraction_1 = %f;\n",        ek_parameters.reaction_fraction_0);
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

int ek_set_lb_force(double* ext_force) {
  for (int i = 0; i < 3; i++)
    ek_parameters.lb_force[i] = ext_force[i];
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
#ifdef EK_ELECTROSTATIC_COUPLING
int ek_set_electrostatics_coupling( bool electrostatics_coupling ) {
  ek_parameters.es_coupling = electrostatics_coupling;
  return 0;
}
#endif

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


int ek_set_stencil( int stencil ) {
  if(!ek_parameters.fluidcoupling_ideal_contribution)
    return 1; //combination not implemented

  ek_parameters.stencil = stencil;
  return 0;
}


int ek_set_advection( bool advection ) {
  ek_parameters.advection = advection;
  return 0;
}


int ek_set_fluidcoupling( bool ideal_contribution ) {
  if(ek_parameters.stencil != 0)
    return 1; //combination not implemented

  ek_parameters.fluidcoupling_ideal_contribution = ideal_contribution;
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

float ek_calculate_net_charge() {
  float charge = 0.0f;
  cuda_safe_mem( cudaMemcpyToSymbol(charge_gpu, &charge, sizeof(float), 0, cudaMemcpyHostToDevice) );

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( ek_parameters.number_of_nodes +
      threads_per_block * blocks_per_grid_y - 1
    ) / ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );

  KERNELCALL( ek_calculate_system_charge, dim_grid, threads_per_block, () );

  cuda_safe_mem( cudaMemcpyFromSymbol(&charge, charge_gpu, sizeof(float), 0, cudaMemcpyDeviceToHost ) );

  return charge;
}

int ek_neutralize_system(int species) {
  int species_index = ek_parameters.species_index[species];

  if(species_index == -1)
    return 1;

  if(ek_parameters.valency[species_index] == 0.0f)
    return 2;

  float compensating_species_density = 0.0f;

#ifndef EK_BOUNDARIES
  for(int i = 0; i < ek_parameters.number_of_species; i++)
    compensating_species_density += ek_parameters.density[i] * ek_parameters.valency[i];

  compensating_species_density = ek_parameters.density[species_index] - compensating_species_density / ek_parameters.valency[species_index];
#else
  float charge = ek_calculate_net_charge();

  compensating_species_density = ek_parameters.density[species_index] - (charge / ek_parameters.valency[species_index]) / (ek_parameters.agrid * ek_parameters.agrid * ek_parameters.agrid * double(ek_parameters.number_of_nodes-ek_parameters.number_of_boundary_nodes));
#endif

  if(compensating_species_density < 0.0f)
    return 3;

  ek_parameters.density[species_index] = compensating_species_density;

  return 0;
}

int ek_save_checkpoint(char* filename) {
  std::string fname(filename);
  std::ofstream fout((const char *) (fname + ".ek").c_str(), std::ofstream::binary);
  ekfloat* densities = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes *
                                      sizeof( ekfloat )                 );

  for(int i = 0; i < ek_parameters.number_of_species; i++)
  {
    cuda_safe_mem( cudaMemcpy( densities, 
                               ek_parameters.rho[i],
                               ek_parameters.number_of_nodes * sizeof(ekfloat),
                               cudaMemcpyDeviceToHost )
                 );

    if(!fout.write((char*) densities, sizeof(ekfloat)*ek_parameters.number_of_nodes))
    {
      free(densities);
      fout.close();
      return 1;
    }
  }

  free(densities);
  fout.close();
 
  lb_lbfluid_save_checkpoint_wrapper((char*) (fname + ".lb").c_str(), 1);

  return 0;
}

int ek_load_checkpoint(char* filename) {
  std::string fname(filename);
  std::ifstream fin((const char *) (fname + ".ek").c_str(), std::ifstream::binary);
  ekfloat* densities = (ekfloat*) Utils::malloc( ek_parameters.number_of_nodes *
                                      sizeof( ekfloat )                 );

  for(int i = 0; i < ek_parameters.number_of_species; i++)
  {
    if(!fin.read((char*) densities, sizeof(ekfloat)*ek_parameters.number_of_nodes))
    {
      free(densities);
      fin.close();
      return 1;
    }

    cuda_safe_mem( cudaMemcpy( ek_parameters.rho[i],
                               densities, 
                               ek_parameters.number_of_nodes * sizeof(ekfloat),
                               cudaMemcpyHostToDevice )
                 );
  }

  free(densities);
  fin.close();

  lb_lbfluid_load_checkpoint_wrapper((char*) (fname + ".lb").c_str(), 1);

  ek_integrate_electrostatics();
  
  return 0;
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

int ek_tag_reaction_nodes( LB_Boundary *boundary, char reaction_type )
{

#ifdef EK_BOUNDARIES
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
        std::ostringstream msg;
        msg << "lbboundary type " << boundary->type << " not implemented";
        runtimeError(msg);
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
