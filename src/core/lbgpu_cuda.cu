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

/** \file lbgpu_cuda.cu
 *
 * Cuda (.cu) file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.hpp.
 */

#include "config.hpp"

#ifdef LB_GPU
#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <vector>
#include <cassert>

#include "electrokinetics.hpp"
#include "electrokinetics_pdb_parse.hpp"
#include "lbgpu.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "debug.hpp"

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif

int extended_values_flag=0; /* TODO: this has to be set to one by
                               appropriate functions if there is 
                               the need to compute pi at every 
                               step (e.g. moving boundaries)*/

/**defining structures residing in global memory */

/** device_rho_v: struct for hydrodynamic fields: this is for internal use 
    (i.e. stores values in LB units) and should not used for 
    printing values  */
static LB_rho_v_gpu *device_rho_v= nullptr;

/** device_rho_v_pi: extended struct for hydrodynamic fields: this is the interface
    and stores values in MD units. It should not be used
    as an input for any LB calculations. TODO: This structure is not yet 
    used, and it is here to allow access to the stress tensor at any
    timestep, e.g. for future implementations of moving boundary codes */
static LB_rho_v_pi_gpu *device_rho_v_pi= nullptr;

/** print_rho_v_pi: struct for hydrodynamic fields: this is the interface
    and stores values in MD units. It should not used
    as an input for any LB calculations. TODO: in the future,
    one might want to have several structures for printing 
    separately rho, v, pi without having to compute/store 
    the complete set. */
static LB_rho_v_pi_gpu *print_rho_v_pi= nullptr;

/** structs for velocity densities */
static LB_nodes_gpu nodes_a = { nullptr, nullptr, nullptr};
static LB_nodes_gpu nodes_b = { nullptr, nullptr, nullptr};;
/** struct for node force density*/

LB_node_force_density_gpu node_f = {nullptr, nullptr} ;

static LB_extern_nodeforcedensity_gpu *extern_node_force_densities = nullptr;

#ifdef LB_BOUNDARIES_GPU
static float* lb_boundary_force = nullptr;

static float* lb_boundary_velocity = nullptr;

/** pointer for bound index array*/
static int *boundary_node_list;
static int *boundary_index_list;
static size_t size_of_boundindex;
#endif

EK_parameters* lb_ek_parameters_gpu;

/** pointers for additional cuda check flag*/
static int *gpu_check = nullptr;
static int *h_gpu_check = nullptr;

static unsigned int intflag = 1;
LB_nodes_gpu *current_nodes = nullptr;
/**defining size values for allocating global memory */
static size_t size_of_rho_v;
static size_t size_of_rho_v_pi;
static size_t size_of_extern_node_force_densities;

/**parameters residing in constant memory */
static __device__ __constant__ LB_parameters_gpu para;
static const float c_sound_sq = 1.0f/3.0f;

/*-------------------------------------------------------*/
/*********************************************************/
/** \name device functions called by kernel functions */
/*********************************************************/
/*-------------------------------------------------------*/

/*-------------------------------------------------------*/

/** atomic add function for sveral cuda architectures 
*/
__device__ inline void atomicadd(float* address, float value){
#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110
#warning Using slower atomicAdd emulation
// float-atomic-add from 
// [url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]
  float old = value;
  while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
#else
#error I need at least compute capability 1.1
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

/**randomgenerator which generates numbers [0,1]
 * @param *rn Pointer to randomnumber array of the local node or particle 
*/
__device__ void random_01(LB_randomnr_gpu *rn){

  const float mxi = 1.0f/(float)(1ul<<31);
  unsigned int curr = rn->seed;

  curr = 1103515245 * curr + 12345;
  rn->randomnr[0] = (float)(curr & ((1ul<<31)-1))*mxi;
  curr = 1103515245 * curr + 12345;
  rn->randomnr[1] = (float)(curr & ((1ul<<31)-1))*mxi;
  rn->seed = curr;

}

/**randomgenerator which generates numbers between -2 sigma and 2 sigma in the form of a Gaussian with standard deviation sigma=1.118591404 resulting in 
 * an actual standard deviation of 1.
 * @param *rn Pointer to randomnumber array of the local node or particle 
*/
__device__ void gaussian_random_cut(LB_randomnr_gpu *rn){

  float x1, x2;
  float r2, fac;
  /** On every second call two gaussian random numbers are calculated
   via the Box-Muller transformation.*/
  /** draw two uniform random numbers in the unit circle */
  do {
    random_01(rn);
    x1 = 2.0f*rn->randomnr[0] - 1.0f;
    x2 = 2.0f*rn->randomnr[1] - 1.0f;
    r2 = x1*x1 + x2*x2;
  } while (r2 >= 1.0f || r2 == 0.0f);

  /** perform Box-Muller transformation and cutoff the ends and replace with flat noise */
  /*
  fac = sqrtf(-2.0f*__logf(r2)/r2)*1.118591404f;
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  random_01(rn);
  if ( fabs(rn->randomnr[0]) > 2.0f*1.118591404f) {
    rn->randomnr[0] = (2.0f*rn->randomnr[0]-1.0f)*2.0f*1.118591404f;
  }
  if ( fabs(rn->randomnr[1]) > 2.0f*1.118591404f ) {
    rn->randomnr[0] = (2.0f*rn->randomnr[1]-1.0f)*2.0f*1.118591404f;
  }
  */
  
  fac = sqrtf(-2.0f*__logf(r2)/r2)*1.042267973f;
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  if ( fabs(rn->randomnr[0]) > 2.0f*1.042267973f) {
    if ( rn->randomnr[0] > 0 ) rn->randomnr[0] = 2.0f*1.042267973f;
    else rn->randomnr[0] = -2.0f*1.042267973f;
  }
  if ( fabs(rn->randomnr[1]) > 2.0f*1.042267973f ) {
    if ( rn->randomnr[1] > 0 ) rn->randomnr[1] = 2.0f*1.042267973f;
    else rn->randomnr[1] = -2.0f*1.042267973f;
  }
}

/** gaussian random nummber generator for thermalisation
 * @param *rn Pointer to randomnumber array of the local node node or particle 
*/
__device__ void gaussian_random(LB_randomnr_gpu *rn){

  float x1, x2;
  float r2, fac;
  /** On every second call two gaussian random numbers are calculated
   via the Box-Muller transformation.*/
  /** draw two uniform random numbers in the unit circle */
  do {
    random_01(rn);
    x1 = 2.0f*rn->randomnr[0]-1.0f;
    x2 = 2.0f*rn->randomnr[1]-1.0f;
    r2 = x1*x1 + x2*x2;
  } while (r2 >= 1.0f || r2 == 0.0f);

  /** perform Box-Muller transformation */
  fac = sqrtf(-2.0f*__logf(r2)/r2);
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  
}
/* wrapper */
__device__ void random_wrapper(LB_randomnr_gpu *rn) { 

#if defined(FLATNOISE)
#define sqrt12 3.46410161514f
  random_01(rn);
  rn->randomnr[0]-=0.5f;
  rn->randomnr[0]*=sqrt12;
  rn->randomnr[1]-=0.5f;
  rn->randomnr[1]*=sqrt12;
#elif defined(GAUSSRANDOMCUT)
  gaussian_random_cut(rn);
#elif defined(GAUSSRANDOM)
  gaussian_random(rn);
#else
#error No noise type defined for the GPU LB
#endif  
  
}


/**tranformation from 1d array-index to xyz
 * @param index   node index / thread index (Input)
 * @param xyz     Pointer to calculated xyz array (Output)
 */
template < typename T >
__device__ void index_to_xyz(T index, T* xyz){
  xyz[0] = index%para.dim_x;
  index /= para.dim_x;
  xyz[1] = index%para.dim_y;
  index /= para.dim_y;
  xyz[2] = index;
}


/**tranformation from xyz to 1d array-index
 * @param xyz     Pointer xyz array (Input)
 * @param index   Calculated node index / thread index (Output)
 */
template < typename T >
__device__ T xyz_to_index(T* xyz){
  T x = (xyz[0] + para.dim_x) % para.dim_x;
  T y = (xyz[1] + para.dim_y) % para.dim_y;
  T z = (xyz[2] + para.dim_z) % para.dim_z;
  return x + para.dim_x*(y + para.dim_y*z);
}


/**calculation of the modes from the velocity densities (space-transform.)
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param index   Node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Output)
*/
__device__ void calc_m_from_n(LB_nodes_gpu n_a, unsigned int index, float *mode){

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {
    // The following convention is used:
    // The $\hat{c}_i$ form B. Duenweg's paper are given by:

    /* c_0  = { 0, 0, 0}
       c_1  = { 1, 0, 0}
       c_2  = {-1, 0, 0}
       c_3  = { 0, 1, 0}
       c_4  = { 0,-1, 0}
       c_5  = { 0, 0, 1}
       c_6  = { 0, 0,-1}
       c_7  = { 1, 1, 0}
       c_8  = {-1,-1, 0}
       c_9  = { 1,-1, 0}
       c_10 = {-1, 1, 0}
       c_11 = { 1, 0, 1}
       c_12 = {-1, 0,-1}
       c_13 = { 1, 0,-1}
       c_14 = {-1, 0, 1}
       c_15 = { 0, 1, 1}
       c_16 = { 0,-1,-1}
       c_17 = { 0, 1,-1}
       c_18 = { 0,-1, 1} */

    // The basis vectors (modes) are constructed as follows
    // $m_k = \sum_{i} e_{ki} n_{i}$, where the $e_{ki}$ form a 
    // linear transformation (matrix) that is given by

    /* $e{ 0,i} = 1$
       $e{ 1,i} = c_{i,x}$
       $e{ 2,i} = c_{i,y}$
       $e{ 3,i} = c_{i,z}$
       $e{ 4,i} = c_{i}^2 - 1$
       $e{ 5,i} = c_{i,x}^2 - c_{i,y}^2$
       $e{ 6,i} = c_{i}^2 - 3*c_{i,z}^2$
       $e{ 7,i} = c_{i,x}*c_{i,y}$
       $e{ 8,i} = c_{i,x}*c_{i,z}$
       $e{ 9,i} = c_{i,y}*c_{i,z}$
       $e{10,i} = (3*c_{i}^2 - 5)*c_{i,x}$
       $e{11,i} = (3*c_{i}^2 - 5)*c_{i,y}$
       $e{12,i} = (3*c_{i}^2 - 5)*c_{i,z}$
       $e{13,i} = (c_{i,y}^2 - c_{i,z}^2)*c_{i,x}$
       $e{14,i} = (c_{i,x}^2 - c_{i,z}^2)*c_{i,y}$
       $e{15,i} = (c_{i,x}^2 - c_{i,y}^2)*c_{i,z}$
       $e{16,i} = 3*c_{i}^2^2 - 6*c_{i}^2 + 1$
       $e{17,i} = (2*c_{i}^2 - 3)*(c_{i,x}^2 - c_{i,y}^2)$
       $e{18,i} = (2*c_{i}^2 - 3)*(c_{i}^2 - 3*c_{i,z}^2)$ */

    // Such that the transformation matrix is given by

    /* {{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1}, 
        { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1}, 
        {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1}, 
        { 0, 1, 1, 1, 1,-2,-2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1}, 
        { 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1}, 
        { 0,-2, 2, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0,-2, 2, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1}, 
        { 0, 0, 0, 0, 0,-2, 2, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1}, 
        { 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1, 1,-1}, 
        { 1,-2,-2,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1}, 
        { 0,-1,-1,-1,-1, 2, 2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1}} */

    // With weights 

    /* q^{c_{i}} = { 1/3, 1/18, 1/18, 1/18,
                    1/18, 1/18, 1/18, 1/36,
                    1/36, 1/36, 1/36, 1/36, 
                    1/36, 1/36, 1/36, 1/36, 
                    1/36, 1/36, 1/36 } */

    // Which makes the transformation satisfy the following
    // orthogonality condition:
    // \sum_{i} q^{c_{i}} e_{ki} e_{li} = w_{k} \delta_{kl},
    // where the weights are:

    /* w_{i} = {  1, 1/3, 1/3, 1/3,
                2/3, 4/9, 4/3, 1/9,
                1/9, 1/9, 2/3, 2/3,
                2/3, 2/9, 2/9, 2/9, 
                  2, 4/9, 4/3 } */

    // mass mode

    mode[0 + ii * LBQ] =   n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index];

    // momentum modes

    mode[1 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[2 + ii * LBQ] =   (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[3 + ii * LBQ] =   (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    // stress modes
    mode[4 + ii * LBQ] = - n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index];

    mode[5 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[6 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index])
                         - 2.0f*( 
                                    (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                  - (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                                  - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                                );

    mode[7 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[8 + ii * LBQ] =   (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[9 + ii * LBQ] =   (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    // kinetic modes

    mode[10 + ii * LBQ] = - 2.0f*(n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[11 + ii * LBQ] = - 2.0f*(n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[12 + ii * LBQ] = - 2.0f*(n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[13 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[14 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[15 + ii * LBQ] =   (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[16 + ii * LBQ] =   n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]
                          - 2.0f*(
                                     (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                 );

    mode[17 + ii * LBQ] = - (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[18 + ii * LBQ] = - (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index])
                          + 2.0f*(
                                     (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                                 );
  }
}

__device__ void reset_LB_force_densities(unsigned int index, LB_node_force_density_gpu node_f, bool buffer = true) {

  float force_factor=powf(para.agrid,2)*para.tau*para.tau;
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {

#if defined(VIRTUAL_SITES_INERTIALESS_TRACERS) || defined(EK_DEBUG)
// Store backup of the node forces
    if (buffer)
    {
      node_f.force_density_buf[(0 + ii*3 ) * para.number_of_nodes + index] = node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index];
      node_f.force_density_buf[(1 + ii*3 ) * para.number_of_nodes + index] = node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index];
      node_f.force_density_buf[(2 + ii*3 ) * para.number_of_nodes + index] = node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index];
    }
#endif

#ifdef EXTERNAL_FORCES
      if(para.external_force_density)
      {
        node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] = para.ext_force_density[0 + ii*3 ]*force_factor;
        node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index] = para.ext_force_density[1 + ii*3 ]*force_factor;
        node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index] = para.ext_force_density[2 + ii*3 ]*force_factor;
      }
      else
      {
        node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      }
#else
      /** reset force */
      node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
#endif
  }
}

__global__ void reset_LB_force_densities_kernel(LB_node_force_density_gpu node_f, bool buffer = true) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if( index < para.number_of_nodes )
    reset_LB_force_densities(index, node_f, buffer);
}

void reset_LB_force_densities_GPU(bool buffer) {
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_LB_force_densities_kernel, dim_grid, threads_per_block, (node_f, buffer));
}


__device__ void update_rho_v(float *mode, unsigned int index, LB_node_force_density_gpu node_f, LB_rho_v_gpu *d_v){

  float Rho_tot=0.0f;
  float u_tot[3]={0.0f,0.0f,0.0f};
  
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      /** re-construct the real density
      * remember that the populations are stored as differences to their
      * equilibrium value */

      d_v[index].rho[ii] = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      Rho_tot  += mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      u_tot[0] += mode[1 + ii * LBQ];
      u_tot[1] += mode[2 + ii * LBQ];
      u_tot[2] += mode[3 + ii * LBQ];

      /** if forces are present, the momentum density is redefined to
      * inlcude one half-step of the force action.  See the
      * Chapman-Enskog expansion in [Ladd & Verberg]. */

      u_tot[0] += 0.5f*node_f.force_density[(0+ii*3)*para.number_of_nodes + index];
      u_tot[1] += 0.5f*node_f.force_density[(1+ii*3)*para.number_of_nodes + index];
      u_tot[2] += 0.5f*node_f.force_density[(2+ii*3)*para.number_of_nodes + index];
  }

  u_tot[0]/=Rho_tot;
  u_tot[1]/=Rho_tot;
  u_tot[2]/=Rho_tot;

  d_v[index].v[0]=u_tot[0]; 
  d_v[index].v[1]=u_tot[1]; 
  d_v[index].v[2]=u_tot[2]; 
}

/**lb_relax_modes, means collision update of the modes
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__device__ void relax_modes(float *mode, unsigned int index, LB_node_force_density_gpu node_f, LB_rho_v_gpu *d_v){
  float u_tot[3]={0.0f,0.0f,0.0f};

  update_rho_v(mode, index, node_f, d_v);

  u_tot[0]=d_v[index].v[0];  
  u_tot[1]=d_v[index].v[1];  
  u_tot[2]=d_v[index].v[2];  
 
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      float Rho; float j[3]; float modes_from_pi_eq[6];

      Rho = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid ;
      j[0] = Rho * u_tot[0];
      j[1] = Rho * u_tot[1];
      j[2] = Rho * u_tot[2];

      /** equilibrium part of the stress modes (eq13 schiller)*/

      modes_from_pi_eq[0] = ((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2]))/Rho;
      modes_from_pi_eq[1] = ((j[0]*j[0])-(j[1]*j[1]))/Rho;
      modes_from_pi_eq[2] = (((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2])) - 3.0f*(j[2]*j[2]))/Rho;
      modes_from_pi_eq[3] = j[0]*j[1]/Rho;
      modes_from_pi_eq[4] = j[0]*j[2]/Rho;
      modes_from_pi_eq[5] = j[1]*j[2]/Rho;
 
      /** in Shan-Chen we have to relax the momentum modes as well using the mobility, but
          the total momentum is conserved */  

#ifdef SHANCHEN
      mode[1 + ii * LBQ] = j[0] + para.gamma_mobility[0]*(mode[1 + ii * LBQ] - j[0]);
      mode[2 + ii * LBQ] = j[1] + para.gamma_mobility[0]*(mode[2 + ii * LBQ] - j[1]);
      mode[3 + ii * LBQ] = j[2] + para.gamma_mobility[0]*(mode[3 + ii * LBQ] - j[2]);
#endif
 
      /** relax the stress modes (eq14 schiller)*/

      mode[4 + ii * LBQ] = modes_from_pi_eq[0] +  para.gamma_bulk[ii]*(mode[4 + ii * LBQ] - modes_from_pi_eq[0]);
      mode[5 + ii * LBQ] = modes_from_pi_eq[1] + para.gamma_shear[ii]*(mode[5 + ii * LBQ] - modes_from_pi_eq[1]);
      mode[6 + ii * LBQ] = modes_from_pi_eq[2] + para.gamma_shear[ii]*(mode[6 + ii * LBQ] - modes_from_pi_eq[2]);
      mode[7 + ii * LBQ] = modes_from_pi_eq[3] + para.gamma_shear[ii]*(mode[7 + ii * LBQ] - modes_from_pi_eq[3]);
      mode[8 + ii * LBQ] = modes_from_pi_eq[4] + para.gamma_shear[ii]*(mode[8 + ii * LBQ] - modes_from_pi_eq[4]);
      mode[9 + ii * LBQ] = modes_from_pi_eq[5] + para.gamma_shear[ii]*(mode[9 + ii * LBQ] - modes_from_pi_eq[5]);
    
      /** relax the ghost modes (project them out) */
      /** ghost modes have no equilibrium part due to orthogonality */

      mode[10 + ii * LBQ] =  para.gamma_odd[ii]*mode[10 + ii * LBQ];
      mode[11 + ii * LBQ] =  para.gamma_odd[ii]*mode[11 + ii * LBQ];
      mode[12 + ii * LBQ] =  para.gamma_odd[ii]*mode[12 + ii * LBQ];
      mode[13 + ii * LBQ] =  para.gamma_odd[ii]*mode[13 + ii * LBQ];
      mode[14 + ii * LBQ] =  para.gamma_odd[ii]*mode[14 + ii * LBQ];
      mode[15 + ii * LBQ] =  para.gamma_odd[ii]*mode[15 + ii * LBQ];
      mode[16 + ii * LBQ] = para.gamma_even[ii]*mode[16 + ii * LBQ];
      mode[17 + ii * LBQ] = para.gamma_even[ii]*mode[17 + ii * LBQ];
      mode[18 + ii * LBQ] = para.gamma_even[ii]*mode[18 + ii * LBQ];
  }
}


/**thermalization of the modes with gaussian random numbers
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param *rn     Pointer to randomnumber array of the local node
*/
__device__ void thermalize_modes(float *mode, unsigned int index, LB_randomnr_gpu *rn){
  float Rho;
#ifdef SHANCHEN
  float Rho_tot=0.0,c;
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii) { 
      Rho_tot  += mode[0 + ii * LBQ]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
  }
  c = (mode[0 + 0 * LBQ]+ para.rho[0]*para.agrid*para.agrid*para.agrid ) / Rho_tot;
  random_wrapper(rn);
  for(int ii=0;ii<LB_COMPONENTS;++ii) { 
      mode[1 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[0];
      mode[2 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[1];
  }                                      
  random_wrapper(rn);                    
  for(int ii=0;ii<LB_COMPONENTS;++ii)    
      mode[3 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-( para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[0];
#endif
  
  
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  

    /** mass mode */  
    Rho = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;

    /** momentum modes */

    /** stress modes */
    random_wrapper(rn);
    mode[4 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-( para.gamma_bulk[ii]* para.gamma_bulk[ii])))) * rn->randomnr[0];
    mode[5 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[6 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/3.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[0];
    mode[7 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[8 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[0];
    mode[9 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    /** ghost modes */
    random_wrapper(rn);
    mode[10 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[11 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[12 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[13 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[14 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[15 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[16 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii]))))     * rn->randomnr[0];
    mode[17 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/9.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[18 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/3.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii])))) * rn->randomnr[0];
  }
}


/*-------------------------------------------------------*/
/**normalization of the modes need befor backtransformation into velocity space
 * @param mode    Pointer to the local register values mode (Input/Output)
*/
__device__ void normalize_modes(float* mode){
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
    /** normalization factors enter in the back transformation */
    mode[ 0 + ii * LBQ] *= 1.0f;
    mode[ 1 + ii * LBQ] *= 3.0f;
    mode[ 2 + ii * LBQ] *= 3.0f;
    mode[ 3 + ii * LBQ] *= 3.0f;
    mode[ 4 + ii * LBQ] *= 3.0f/2.0f;
    mode[ 5 + ii * LBQ] *= 9.0f/4.0f;
    mode[ 6 + ii * LBQ] *= 3.0f/4.0f;
    mode[ 7 + ii * LBQ] *= 9.0f;
    mode[ 8 + ii * LBQ] *= 9.0f;
    mode[ 9 + ii * LBQ] *= 9.0f;
    mode[10 + ii * LBQ] *= 3.0f/2.0f;
    mode[11 + ii * LBQ] *= 3.0f/2.0f;
    mode[12 + ii * LBQ] *= 3.0f/2.0f;
    mode[13 + ii * LBQ] *= 9.0f/2.0f;
    mode[14 + ii * LBQ] *= 9.0f/2.0f;
    mode[15 + ii * LBQ] *= 9.0f/2.0f;
    mode[16 + ii * LBQ] *= 1.0f/2.0f;
    mode[17 + ii * LBQ] *= 9.0f/4.0f;
    mode[18 + ii * LBQ] *= 3.0f/4.0f;
  }
}



/*-------------------------------------------------------*/
/**backtransformation from modespace to desityspace and streaming with the push method using pbc
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param *n_b    Pointer to local node residing in array b (Output)
*/
__device__ void calc_n_from_modes_push(LB_nodes_gpu n_b, float *mode, unsigned int index){

  unsigned int xyz[3];
  index_to_xyz(index, xyz);
  unsigned int x = xyz[0];
  unsigned int y = xyz[1];
  unsigned int z = xyz[2];

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {
 
    n_b.vd[(0 + ii*LBQ ) * para.number_of_nodes + x 
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*z] = 
      1.0f/3.0f * (mode[0 + ii * LBQ] - mode[4 + ii * LBQ] + mode[16 + ii * LBQ]);

    n_b.vd[(1 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*y 
                                                + para.dim_x*para.dim_y*z] = 
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     - 2.0f*(mode[10 + ii * LBQ] + mode[16 + ii * LBQ])
                   );

    n_b.vd[(2 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     + 2.0f*(mode[10 + ii * LBQ] - mode[16 + ii * LBQ])
                   );

    n_b.vd[(3 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*((y+1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     - 2.0f*(mode[11 + ii * LBQ] + mode[16 + ii * LBQ])
                   );

    n_b.vd[(4 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     + 2.0f*(mode[11 + ii * LBQ] - mode[16 + ii * LBQ])
                   );

    n_b.vd[(5 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/18.0f * (
                       mode[0 + ii * LBQ] + mode[3 + ii * LBQ]
                     - 2.0f*(   mode[ 6 + ii * LBQ] + mode[12 + ii * LBQ]
                              + mode[16 + ii * LBQ] - mode[18 + ii * LBQ])
                   );

    n_b.vd[(6 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/18.0f * (
                       mode[0 + ii * LBQ] - mode[3 + ii * LBQ]
                     - 2.0f*(   mode[6 + ii * LBQ] - mode[12 + ii * LBQ]
                              + mode[16 + ii * LBQ] - mode[18 + ii * LBQ])
                   );

    n_b.vd[(7 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*((y+1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] + mode[ 7 + ii * LBQ]
                     + mode[10 + ii * LBQ] + mode[11 + ii * LBQ]
                     + mode[13 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(8 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     - mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] + mode[ 7 + ii * LBQ]
                     - mode[10 + ii * LBQ] - mode[11 + ii * LBQ]
                     - mode[13 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(9 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     - mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] - mode[ 7 + ii * LBQ]
                     + mode[10 + ii * LBQ] - mode[11 + ii * LBQ]
                     + mode[13 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(10 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*z] = 
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] - mode[ 7 + ii * LBQ]
                     - mode[10 + ii * LBQ] + mode[11 + ii * LBQ]
                     - mode[13 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(11 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 8 + ii * LBQ] + mode[10 + ii * LBQ]
                     + mode[12 + ii * LBQ] - mode[13 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(12 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 8 + ii * LBQ] - mode[10 + ii * LBQ]
                     - mode[12 + ii * LBQ] + mode[13 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(13 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 8 + ii * LBQ] + mode[10 + ii * LBQ]
                     - mode[12 + ii * LBQ] - mode[13 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(14 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 8 + ii * LBQ] - mode[10 + ii * LBQ]
                     + mode[12 + ii * LBQ] + mode[13 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(15 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 9 + ii * LBQ] + mode[11 + ii * LBQ]
                     + mode[12 + ii * LBQ] - mode[14 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(16 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 9 + ii * LBQ] - mode[11 + ii * LBQ]
                     - mode[12 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(17 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 9 + ii * LBQ] + mode[11 + ii * LBQ]
                     - mode[12 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(18 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 9 + ii * LBQ] - mode[11 + ii * LBQ]
                     + mode[12 + ii * LBQ] + mode[14 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );
  }
}



/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 * @param index   node index / thread index (Input)
 * @param n_curr  Pointer to local node which receives the current node field (Input)
 * @param lb_boundary_velocity    The constant velocity at the boundary, set by the user (Input)
 * @param lb_boundary_force       The force on the boundary nodes (Output)
*/
__device__ void bounce_back_boundaries(LB_nodes_gpu n_curr, unsigned int index, \
    float* lb_boundary_velocity, float* lb_boundary_force){
    
  unsigned int xyz[3];
  int c[3];
  float v[3];
  float shift, weight, pop_to_bounce_back;
  float boundary_force[3] = {0.0f,0.0f,0.0f};
  size_t to_index, to_index_x, to_index_y, to_index_z;
  int population, inverse;
  int boundary_index;

  boundary_index= n_curr.boundary[index];
  if(boundary_index != 0)
  {
    
    v[0]=lb_boundary_velocity[3*(boundary_index-1)+0];
    v[1]=lb_boundary_velocity[3*(boundary_index-1)+1];
    v[2]=lb_boundary_velocity[3*(boundary_index-1)+2];

    index_to_xyz(index, xyz);

    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];

    /* CPU analog of shift:
       lbpar.agrid*lbpar.agrid*lbpar.agrid*lbpar.rho*2*lbmodel.c[i][l]*lb_boundaries[lbfields[k].boundary-1].velocity[l] */
  
    /** store vd temporary in second lattice to avoid race conditions */

// TODO : PUT IN EQUILIBRIUM CONTRIBUTION TO THE BOUNCE-BACK DENSITY FOR THE BOUNDARY FORCE
// TODO : INITIALIZE BOUNDARY FORCE PROPERLY, HAS NONZERO ELEMENTS IN FIRST STEP
// TODO : SET INTERNAL BOUNDARY NODE VALUES TO ZERO

#ifndef SHANCHEN

#define BOUNCEBACK()  \
  shift = 2.0f*para.agrid*para.agrid*para.rho[0]*3.0f*weight*para.tau*(v[0]*c[0] + v[1]*c[1] + v[2]*c[2]); \
  pop_to_bounce_back =  n_curr.vd[population*para.number_of_nodes + index ]; \
  to_index_x = (x+c[0]+para.dim_x)%para.dim_x; \
  to_index_y = (y+c[1]+para.dim_y)%para.dim_y; \
  to_index_z = (z+c[2]+para.dim_z)%para.dim_z; \
  to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z; \
  if ( n_curr.boundary[to_index] == 0) \
  { \
    boundary_force[0] += (2.0f*pop_to_bounce_back+shift)*c[0]/para.tau/para.tau/para.agrid; \
    boundary_force[1] += (2.0f*pop_to_bounce_back+shift)*c[1]/para.tau/para.tau/para.agrid; \
    boundary_force[2] += (2.0f*pop_to_bounce_back+shift)*c[2]/para.tau/para.tau/para.agrid; \
     n_curr.vd[inverse*para.number_of_nodes + to_index ] = pop_to_bounce_back + shift; \
  }
#else 

#define BOUNCEBACK()  \
  for(int component=0; component<LB_COMPONENTS;component++){\
     shift = 2.0f*para.agrid*para.agrid*para.rho[component]*3.0f*weight*para.tau*(v[0]*c[0] + v[1]*c[1] + v[2]*c[2]); \
     pop_to_bounce_back =  n_curr.vd[(population+component*LBQ)*para.number_of_nodes + index ]; \
     to_index_x = (x+c[0]+para.dim_x)%para.dim_x; \
     to_index_y = (y+c[1]+para.dim_y)%para.dim_y; \
     to_index_z = (z+c[2]+para.dim_z)%para.dim_z; \
     to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z; \
     if ( n_curr.boundary[to_index] == 0) \
     { \
       boundary_force[0] += (2.0f*pop_to_bounce_back+shift)*c[0]/para.tau/para.tau/para.agrid; \
       boundary_force[1] += (2.0f*pop_to_bounce_back+shift)*c[1]/para.tau/para.tau/para.agrid; \
       boundary_force[2] += (2.0f*pop_to_bounce_back+shift)*c[2]/para.tau/para.tau/para.agrid; \
       n_curr.vd[(inverse+component*LBQ)*para.number_of_nodes + to_index ] = pop_to_bounce_back + shift; \
       n_curr.vd[(inverse+component*LBQ)*para.number_of_nodes + to_index ] = pop_to_bounce_back + shift; \
     } \
  }
#endif

    // the resting population does nothing, i.e., population 0.
    c[0]= 1;c[1]= 0;c[2]= 0; weight=1./18.; population= 2; inverse= 1; 
    BOUNCEBACK();
    
    c[0]=-1;c[1]= 0;c[2]= 0; weight=1./18.; population= 1; inverse= 2; 
    BOUNCEBACK();
    
    c[0]= 0;c[1]= 1;c[2]= 0; weight=1./18.; population= 4; inverse= 3; 
    BOUNCEBACK();

    c[0]= 0;c[1]=-1;c[2]= 0; weight=1./18.; population= 3; inverse= 4; 
    BOUNCEBACK();
    
    c[0]= 0;c[1]= 0;c[2]= 1; weight=1./18.; population= 6; inverse= 5; 
    BOUNCEBACK();

    c[0]= 0;c[1]= 0;c[2]=-1; weight=1./18.; population= 5; inverse= 6; 
    BOUNCEBACK(); 
    
    c[0]= 1;c[1]= 1;c[2]= 0; weight=1./36.; population= 8; inverse= 7; 
    BOUNCEBACK();
    
    c[0]=-1;c[1]=-1;c[2]= 0; weight=1./36.; population= 7; inverse= 8; 
    BOUNCEBACK();
    
    c[0]= 1;c[1]=-1;c[2]= 0; weight=1./36.; population=10; inverse= 9; 
    BOUNCEBACK();

    c[0]=-1;c[1]= 1;c[2]= 0; weight=1./36.; population= 9; inverse=10; 
    BOUNCEBACK();
    
    c[0]= 1;c[1]= 0;c[2]= 1; weight=1./36.; population=12; inverse=11; 
    BOUNCEBACK();
    
    c[0]=-1;c[1]= 0;c[2]=-1; weight=1./36.; population=11; inverse=12; 
    BOUNCEBACK();

    c[0]= 1;c[1]= 0;c[2]=-1; weight=1./36.; population=14; inverse=13; 
    BOUNCEBACK();
    
    c[0]=-1;c[1]= 0;c[2]= 1; weight=1./36.; population=13; inverse=14; 
    BOUNCEBACK();

    c[0]= 0;c[1]= 1;c[2]= 1; weight=1./36.; population=16; inverse=15; 
    BOUNCEBACK();
    
    c[0]= 0;c[1]=-1;c[2]=-1; weight=1./36.; population=15; inverse=16; 
    BOUNCEBACK();
    
    c[0]= 0;c[1]= 1;c[2]=-1; weight=1./36.; population=18; inverse=17; 
    BOUNCEBACK();
    
    c[0]= 0;c[1]=-1;c[2]= 1; weight=1./36.; population=17; inverse=18; 
    BOUNCEBACK();  
    
    atomicadd(&lb_boundary_force[3*( n_curr.boundary[index]-1)+0], boundary_force[0]);
    atomicadd(&lb_boundary_force[3*( n_curr.boundary[index]-1)+1], boundary_force[1]);
    atomicadd(&lb_boundary_force[3*( n_curr.boundary[index]-1)+2], boundary_force[2]);
  }
}

/** add of (external) forces within the modespace, needed for particle-interaction
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__device__ void apply_forces(unsigned int index, float *mode, LB_node_force_density_gpu node_f, LB_rho_v_gpu *d_v) {
  
  float u[3]={0.0f,0.0f,0.0f},
        C[6]={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  /* Note: the values d_v were calculated in relax_modes() */

  u[0]=d_v[index].v[0]; 
  u[1]=d_v[index].v[1]; 
  u[2]=d_v[index].v[2]; 

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
       C[0] += (1.0f + para.gamma_bulk[ii])*u[0]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[1]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[2] += (1.0f + para.gamma_bulk[ii])*u[1]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[1]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[5] += (1.0f + para.gamma_bulk[ii])*u[2]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] 
                                                                       + u[1]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[1] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[0]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[1]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index]
                                                     );

       C[3] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[0]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[2]*node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index]
                                                     );

       C[4] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[1]*node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[2]*node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index]
                                                     );
  }

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
      /** update momentum modes */
#ifdef SHANCHEN
      float mobility_factor=1.0f/2.0f*(1.0f+para.gamma_mobility[0]);
#else
      float mobility_factor=1.0f;
#endif 
 /** update momentum modes */
      mode[1 + ii * LBQ] += mobility_factor * node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index];
      mode[2 + ii * LBQ] += mobility_factor * node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index];
      mode[3 + ii * LBQ] += mobility_factor * node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index];

      /** update stress modes */
      mode[4 + ii * LBQ] += C[0] + C[2] + C[5];
      mode[5 + ii * LBQ] += C[0] - C[2];
      mode[6 + ii * LBQ] += C[0] + C[2] - 2.0f*C[5];
      mode[7 + ii * LBQ] += C[1];
      mode[8 + ii * LBQ] += C[3];
      mode[9 + ii * LBQ] += C[4];
    
  }

//#if !defined(IMMERSED_BOUNDARY)
  // This must not be done here since we need the forces after LB update for the velocity interpolation
  // It is done by calling IBM_ResetLBForces_GPU from integrate_vv
  reset_LB_force_densities(index, node_f);
//#endif

#ifdef SHANCHEN
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
     node_f.force_density[(0 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce_density[(0+ii*3)*para.number_of_nodes + index];
     node_f.force_density[(1 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce_density[(1+ii*3)*para.number_of_nodes + index];
     node_f.force_density[(2 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce_density[(2+ii*3)*para.number_of_nodes + index];
  }
#endif
}

/**function used to calculate hydrodynamic fields in MD units.
 * @param n_a     Pointer to local node residing in array a for boundary flag(Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param d_p_v   Pointer to local print values (Output)
 * @param d_v     Pointer to local device values (Input)
 * @param node_f  Pointer to local node force (Input)
 * @param index   node index / thread index (Input)
 * @param print_index   node index / thread index (Output)
*/
__device__ void calc_values_in_MD_units(LB_nodes_gpu n_a, float *mode, LB_rho_v_pi_gpu *d_p_v, LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f, unsigned int index, unsigned int print_index) {
  
  float j[3]; 
  float modes_from_pi_eq[6]; 
  float pi[6]={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};

  if(n_a.boundary[index] == 0)
  {
    /* Ensure we are working with the current values of d_v */

    update_rho_v(mode, index, node_f, d_v);

    for(int ii= 0; ii < LB_COMPONENTS; ii++)
    {
      d_p_v[print_index].rho[ii] = d_v[index].rho[ii] / para.agrid / para.agrid / para.agrid;
    }
      
    d_p_v[print_index].v[0] = d_v[index].v[0] * para.agrid / para.tau;
    d_p_v[print_index].v[1] = d_v[index].v[1] * para.agrid / para.tau;
    d_p_v[print_index].v[2] = d_v[index].v[2] * para.agrid / para.tau;

    /* stress calculation */ 
    for(int ii = 0; ii < LB_COMPONENTS; ii++)
    {
      float Rho = d_v[index].rho[ii];
      
      /* note that d_v[index].v[] already includes the 1/2 f term, accounting for the pre- and post-collisional average */

      j[0] = Rho * d_v[index].v[0];
      j[1] = Rho * d_v[index].v[1];
      j[2] = Rho * d_v[index].v[2];

      // equilibrium part of the stress modes, which comes from 
      // the equality between modes and stress tensor components

      /* m4 = trace(pi) - rho
         m5 = pi_xx - pi_yy
         m6 = trace(pi) - 3 pi_zz
         m7 = pi_xy
         m8 = pi_xz
         m9 = pi_yz */

      // and pluggin in the Euler stress for the equilibrium:
      // pi_eq = rho_0*c_s^2*I3 + (j \otimes j)/rho
      // with I3 the 3D identity matrix and
      // rho = \trace(rho_0*c_s^2*I3), which yields

      /* m4_from_pi_eq = j.j
         m5_from_pi_eq = j_x*j_x - j_y*j_y
         m6_from_pi_eq = j.j - 3*j_z*j_z
         m7_from_pi_eq = j_x*j_y
         m8_from_pi_eq = j_x*j_z
         m9_from_pi_eq = j_y*j_z */

      // where the / Rho term has been dropped. We thus obtain: 

      modes_from_pi_eq[0] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] ) / Rho;
      modes_from_pi_eq[1] = ( j[0]*j[0] - j[1]*j[1] ) / Rho;
      modes_from_pi_eq[2] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] - 3.0f*j[2]*j[2] ) / Rho;
      modes_from_pi_eq[3] = j[0]*j[1] / Rho;
      modes_from_pi_eq[4] = j[0]*j[2] / Rho;
      modes_from_pi_eq[5] = j[1]*j[2] / Rho;
     
      /* Now we must predict the outcome of the next collision */
      /* We immediately average pre- and post-collision.  */
      /* TODO: need a reference for this.   */

      mode[4 + ii * LBQ ] = modes_from_pi_eq[0] + (0.5f + 0.5f*para.gamma_bulk[ii]) * (mode[4 + ii * LBQ] - modes_from_pi_eq[0]);
      mode[5 + ii * LBQ ] = modes_from_pi_eq[1] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[5 + ii * LBQ] - modes_from_pi_eq[1]);
      mode[6 + ii * LBQ ] = modes_from_pi_eq[2] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[6 + ii * LBQ] - modes_from_pi_eq[2]);
      mode[7 + ii * LBQ ] = modes_from_pi_eq[3] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[7 + ii * LBQ] - modes_from_pi_eq[3]);
      mode[8 + ii * LBQ ] = modes_from_pi_eq[4] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[8 + ii * LBQ] - modes_from_pi_eq[4]);
      mode[9 + ii * LBQ ] = modes_from_pi_eq[5] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[9 + ii * LBQ] - modes_from_pi_eq[5]);

      // Transform the stress tensor components according to the modes that
      // correspond to those used by U. Schiller. In terms of populations this
      // expression then corresponds exactly to those in Eqs. 116 - 121 in the
      // Duenweg and Ladd paper, when these are written out in populations.
      // But to ensure this, the expression in Schiller's modes has to be different!

      pi[0] += (   2.0f*(mode[0 + ii * LBQ] + mode[4 + ii * LBQ])
                + mode[6 + ii * LBQ] + 3.0f*mode[5 + ii * LBQ] )/6.0f;  // xx
      pi[1] += mode[7 + ii * LBQ];                                      // xy
      pi[2] += (   2.0f*(mode[0 + ii * LBQ] + mode[4 + ii * LBQ])
                + mode[6 + ii * LBQ] - 3.0f*mode[5 + ii * LBQ] )/6.0f;  // yy
      pi[3] += mode[8 + ii * LBQ];                                      // xz
      pi[4] += mode[9 + ii * LBQ];                                      // yz
      pi[5] += (   mode[0 + ii * LBQ] + mode[4 + ii * LBQ]
                - mode[6 + ii * LBQ] )/3.0f;                            // zz

    }
     
    for(int i = 0; i < 6; i++)
    {
      d_p_v[print_index].pi[i] = pi[i] / para.tau
                                       / para.tau
                                       / para.agrid;
    }
  }
  else
  {
    for(int ii = 0; ii < LB_COMPONENTS; ii++)
      d_p_v[print_index].rho[ii] = 0.0f;
     
    for(int i = 0; i < 3; i++)
      d_p_v[print_index].v[i] = 0.0f;

    for(int i = 0; i < 6; i++)
      d_p_v[print_index].pi[i] = 0.0f;
  }
}

/**function used to calculate hydrodynamic fields in MD units.
 * @param mode_single   Pointer to the local register values mode (Input)
 * @param d_v_single    Pointer to local device values (Input)
 * @param rho_out       Pointer to density (Output)
 * @param j_out         Pointer to momentum (Output)
 * @param pi_out        Pointer to pressure tensor (Output)
*/
__device__ void calc_values_from_m_in_LB_units(float *mode_single, LB_rho_v_gpu *d_v_single, float* rho_out, float* j_out, float* pi_out) {

  float modes_from_pi_eq[6];
  float j[6];
  float Rho; 

  // stress calculation

  for(int ii = 0; ii < LB_COMPONENTS; ii++)
  {
    // Set the rho ouput value

    Rho = d_v_single->rho[ii];
    rho_out[ii] = Rho;
    
    // note that d_v_single->v[] already includes the 1/2 f term, 
    // accounting for the pre- and post-collisional average

    j[0] = Rho * d_v_single->v[0];
    j[1] = Rho * d_v_single->v[1];
    j[2] = Rho * d_v_single->v[2];

    j_out[3*ii + 0] = j[0];
    j_out[3*ii + 1] = j[1];
    j_out[3*ii + 2] = j[2];    

    // equilibrium part of the stress modes, which comes from 
    // the equality between modes and stress tensor components

    modes_from_pi_eq[0] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] ) / Rho;
    modes_from_pi_eq[1] = ( j[0]*j[0] - j[1]*j[1] ) / Rho;
    modes_from_pi_eq[2] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] - 3.0f*j[2]*j[2] ) / Rho;
    modes_from_pi_eq[3] = j[0]*j[1] / Rho;
    modes_from_pi_eq[4] = j[0]*j[2] / Rho;
    modes_from_pi_eq[5] = j[1]*j[2] / Rho;
   
    // Now we must predict the outcome of the next collision
    // We immediately average pre- and post-collision.

    mode_single[4 + ii * LBQ ] = modes_from_pi_eq[0] + (0.5f + 0.5f* para.gamma_bulk[ii]) * (mode_single[4 + ii * LBQ] - modes_from_pi_eq[0]);
    mode_single[5 + ii * LBQ ] = modes_from_pi_eq[1] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[5 + ii * LBQ] - modes_from_pi_eq[1]);
    mode_single[6 + ii * LBQ ] = modes_from_pi_eq[2] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[6 + ii * LBQ] - modes_from_pi_eq[2]);
    mode_single[7 + ii * LBQ ] = modes_from_pi_eq[3] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[7 + ii * LBQ] - modes_from_pi_eq[3]);
    mode_single[8 + ii * LBQ ] = modes_from_pi_eq[4] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[8 + ii * LBQ] - modes_from_pi_eq[4]);
    mode_single[9 + ii * LBQ ] = modes_from_pi_eq[5] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[9 + ii * LBQ] - modes_from_pi_eq[5]);

    // Transform the stress tensor components according to the mode_singles.

    pi_out[6*ii + 0] = (   2.0f*(mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ])
                         + mode_single[6 + ii * LBQ] + 3.0f*mode_single[5 + ii * LBQ] )/6.0f;   // xx
    pi_out[6*ii + 1] = mode_single[7 + ii * LBQ];                                               // xy
    pi_out[6*ii + 2] = (   2.0f*(mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ])
                         + mode_single[6 + ii * LBQ] - 3.0f*mode_single[5 + ii * LBQ] )/6.0f;   // yy
    pi_out[6*ii + 3] = mode_single[8 + ii * LBQ];                                               // xz
    pi_out[6*ii + 4] = mode_single[9 + ii * LBQ];                                               // yz
    pi_out[6*ii + 5] = (   mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ]
                         - mode_single[6 + ii * LBQ] )/3.0f;                                    // zz
  }
}

/**function used to calc physical values of every node
 * @param n_a     Pointer to local node residing in array a for boundary flag(Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param d_v     Pointer to local device values (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param index   node index / thread index (Input)
*/

/* FIXME this function is basically un-used, think about removing/replacing it */
__device__ void calc_values(LB_nodes_gpu n_a, float *mode, LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f, unsigned int index){ 

  float Rho_tot=0.0f;
  float u_tot[3]={0.0f,0.0f,0.0f};

  if(n_a.boundary[index] != 1){
      #pragma unroll
      for(int ii=0;ii<LB_COMPONENTS;++ii) { 
          /** re-construct the real density
          * remember that the populations are stored as differences to their
          * equilibrium value */
          d_v[index].rho[ii]= mode[0 + ii * 4]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
          Rho_tot  += mode[0 + ii * 4]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
          u_tot[0] += mode[1 + ii * 4];
          u_tot[1] += mode[2 + ii * 4];
          u_tot[2] += mode[3 + ii * 4];
    
          /** if forces are present, the momentum density is redefined to
          * inlcude one half-step of the force action.  See the
          * Chapman-Enskog expansion in [Ladd & Verberg]. */
    
          u_tot[0] += 0.5f*node_f.force_density[(0+ii*3)*para.number_of_nodes + index];
          u_tot[1] += 0.5f*node_f.force_density[(1+ii*3)*para.number_of_nodes + index];
          u_tot[2] += 0.5f*node_f.force_density[(2+ii*3)*para.number_of_nodes + index];
      }
      u_tot[0]/=Rho_tot;
      u_tot[1]/=Rho_tot;
      u_tot[2]/=Rho_tot;
    
      d_v[index].v[0]=u_tot[0]; 
      d_v[index].v[1]=u_tot[1]; 
      d_v[index].v[2]=u_tot[2]; 
  } else { 
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii) { 
       d_v[index].rho[ii]   = 1.;
    }
    d_v[index].v[0] = 0.0f;
    d_v[index].v[1] = 0.0f; 
    d_v[index].v[2] = 0.0f; 
  }   
}


/** 
 * @param node_index  node index around (8) particle (Input)
 * @param *mode       Pointer to the local register values mode (Output)
 * @param n_a         Pointer to local node residing in array a(Input)
 * @param component_index   Shanchen component index        (Input)
*/
__device__ void calc_mode(float *mode, LB_nodes_gpu n_a, unsigned int node_index, int component_index){

  /** mass mode */
  mode[0] =   n_a.vd[( 0 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + node_index] 
            + n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index];

  /** momentum modes */
  mode[1] =   (n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index]);

  mode[2] =   (n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index]);

  mode[3] =   (n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index]);
}

/*********************************************************/
/** \name interpolation_three_point_coupling */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a             Pointer to local node residing in array a (Input)
 * @param *delta          Pointer for the weighting of particle position (Output)
 * @param *particle_data  Pointer to the particle position and velocity (Input)
 * @param *particle_force Pointer to the particle force (Input)
 * @param part_index      particle id / thread id (Input)
 * @param node_index      node index around (8) particle (Output)
 * @param *d_v            Pointer to local device values
 * @param *interpolated_u Pointer to the interpolated velocity (Output)
*/
__device__ __inline__ void interpolation_three_point_coupling(LB_nodes_gpu n_a, float* particle_position, unsigned int *node_index, LB_rho_v_gpu *d_v, float *delta, float *interpolated_u) {

  int my_center[3];
  float temp_delta[27];
  float mode[19*LB_COMPONENTS];

  /** see Duenweg and Ladd http://arxiv.org/abs/0803.2826 eqn. 301 */
  /** the i index is left node, nearest node, right node */
  for(int i=0; i<3; ++i){
    /** note the -0.5f is to account for the shift of the LB grid relative to the MD */
    float scaledpos = particle_position[i]/para.agrid-0.5f;
    /** the +0.5 is to turn the floorf into a round function */
    my_center[i] = (int)(floorf(scaledpos+0.5f));
    scaledpos = scaledpos-1.0f*my_center[i];
    temp_delta[0+3*i] = ( 5.0f - 3.0f*abs(scaledpos+1.0f) - sqrtf( -2.0f + 6.0f*abs(scaledpos+1.0f) - 3.0f*powf(scaledpos+1.0f,2) ) )/6.0f;
    temp_delta[1+3*i] = ( 1.0f + sqrtf( 1.0f - 3.0f*powf(scaledpos,2) ) )/3.0f;
    temp_delta[2+3*i] = ( 5.0f - 3.0f*abs(scaledpos-1.0f) - sqrtf( -2.0f + 6.0f*abs(scaledpos-1.0f) - 3.0f*powf(scaledpos-1.0f,2) ) )/6.0f;
  }

  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        delta[i+3*j+9*k+13] = temp_delta[i+1] * temp_delta[3+j+1] * temp_delta[6+k+1];
      }
    }
  }

  // modulo for negative numbers is strange at best, shift to make sure we are positive
  int x = my_center[0] + para.dim_x;
  int y = my_center[1] + para.dim_y;
  int z = my_center[2] + para.dim_z;
  /** Here we collect the nodes for the three point coupling scheme (27 nodes in 3d) with the analogous numbering scheme of the two point coupling scheme */
  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        node_index[i+3*j+9*k+13] = (x+i+para.dim_x)%para.dim_x + para.dim_x*((y+j+para.dim_y)%para.dim_y) + para.dim_x*para.dim_y*((z+k+para.dim_z)%para.dim_z);
      }
    }
  }

  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for(int i=0; i<27; ++i){
    float totmass=0.0f;
    calc_m_from_n(n_a,node_index[i],mode);
#pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;ii++){
      totmass+=mode[0]+para.rho[ii]*para.agrid*para.agrid*para.agrid;
    } 
    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary nodes
       to couple with velocity 0 to particles. This is necessary, since boundary nodes
       undergo the same LB dynamics as fluid nodes do. The flow within the boundaries
       does not interact with the physical fluid, since these populations are overwritten
       by the bounce back kernel. Particles close to walls can couple to this unphysical
       flow, though.
    */
    interpolated_u[0] += (mode[1]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] += (mode[2]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += (mode[3]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
  }
}


/*********************************************************/
/** \name calc_viscous_force_three_point_couple */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                Pointer to local node residing in array a (Input)
 * @param *delta             Pointer for the weighting of particle position (Output)
 * @param *delta_j           Pointer for the weighting of particle momentum (Output)
 * @param *particle_position Pointer to the particle position (Input)
 * @param *rn_part           Pointer to randomnumber array of the particle
 * @param node_index         node index around (8) particle (Output)
 * @param *d_v               Pointer to local device values
 * @param flag_cs            Determine if we are at the centre (0, typical) or at the source (1, swimmer only)
*/
__device__ void calc_viscous_force_three_point_couple(LB_nodes_gpu n_a, float *delta, CUDA_particle_data *particle_data, float *particle_force, unsigned int part_index, LB_randomnr_gpu *rn_part, float *delta_j, unsigned int *node_index, LB_rho_v_gpu *d_v, int flag_cs){

  float interpolated_u[3];
  float interpolated_rho[LB_COMPONENTS];
  float viscforce_density[3*LB_COMPONENTS];

  // Zero out workspace
#pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
#pragma unroll
    for(int jj=0; jj<3; ++jj){ 
      viscforce_density[jj+ii*3]=0.0f;
      delta_j[jj+ii*3]  =0.0f;
    }
  }
  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3*part_index+0] = flag_cs * particle_force[3*part_index+0];
  particle_force[3*part_index+1] = flag_cs * particle_force[3*part_index+1];
  particle_force[3*part_index+2] = flag_cs * particle_force[3*part_index+2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) * particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the particle
  position[0] += flag_cs * direction * particle_data[part_index].swim.quatu[0];
  position[1] += flag_cs * direction * particle_data[part_index].swim.quatu[1];
  position[2] += flag_cs * direction * particle_data[part_index].swim.quatu[2];
#endif

  // Do the velocity interpolation
  interpolation_three_point_coupling(n_a, position, node_index, d_v, delta, interpolated_u);

#ifdef ENGINE
  velocity[0] -= (particle_data[part_index].swim.v_swim)*particle_data[part_index].swim.quatu[0];
  velocity[1] -= (particle_data[part_index].swim.v_swim)*particle_data[part_index].swim.quatu[1];
  velocity[2] -= (particle_data[part_index].swim.v_swim)*particle_data[part_index].swim.quatu[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD units
  particle_data[part_index].swim.v_cs[0+3*flag_cs] = interpolated_u[0] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[1+3*flag_cs] = interpolated_u[1] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[2+3*flag_cs] = interpolated_u[2] * para.agrid / para.tau;
#endif

  /* for LB we do not reweight the friction force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    interpolated_rho[ii]=1.0;
  }

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot=0;

#pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    rhotot+=interpolated_rho[ii];
  }


  /* Viscous force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    viscforce_density[0+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[0] - interpolated_u[0]*para.agrid/para.tau)/rhotot;
    viscforce_density[1+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[1] - interpolated_u[1]*para.agrid/para.tau)/rhotot;
    viscforce_density[2+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[2] - interpolated_u[2]*para.agrid/para.tau)/rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
    viscforce_density[0+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[0]/rhotot;
    viscforce_density[1+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[1]/rhotot;
    viscforce_density[2+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[2]/rhotot;
#endif

    /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
#ifdef FLATNOISE
    random_01(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
    viscforce_density[1+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[1]-0.5f);
    random_01(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
#elif defined(GAUSSRANDOMCUT)
    gaussian_random_cut(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce_density[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random_cut(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#elif defined(GAUSSRANDOM)
    gaussian_random(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce_density[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#else
#error No noise type defined for the GPU LB
#endif    
    /** delta_j for transform momentum transfer to lattice units which is done in calc_node_force
      (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    particle_force[3*part_index+0] += (1-flag_cs) * viscforce_density[0+ii*3];
    particle_force[3*part_index+1] += (1-flag_cs) * viscforce_density[1+ii*3];
    particle_force[3*part_index+2] += (1-flag_cs) * viscforce_density[2+ii*3];

    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    delta_j[0+3*ii] -= (1-flag_cs)*viscforce_density[0+ii*3]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= (1-flag_cs)*viscforce_density[1+ii*3]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= (1-flag_cs)*viscforce_density[2+ii*3]*para.time_step*para.tau/para.agrid;
#ifdef ENGINE
    // add swimming force to source position
    delta_j[0+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[0]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[1]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[2]*para.time_step*para.tau/para.agrid;
#endif
  }
}

/**calcutlation of the node force caused by the particles, with atomicadd due to avoiding race conditions 
  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param *delta    Pointer for the weighting of particle position (Input)
 * @param *delta_j    Pointer for the weighting of particle momentum (Input)
 * @param node_index    node index around (8) particle (Input)
 * @param node_f        Pointer to the node force (Output).
*/
__device__ void calc_node_force_three_point_couple(float *delta, float *delta_j, unsigned int *node_index, LB_node_force_density_gpu node_f){
/* TODO: should the drag depend on the density?? */

  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        atomicadd(&(node_f.force_density[0*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[0]));
        atomicadd(&(node_f.force_density[1*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[1]));
        atomicadd(&(node_f.force_density[2*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[2]));
      }
    }
  }
}


/**calculate temperature of the fluid kernel
 * @param *cpu_jsquared   Pointer to result storage value (Output)
 * @param n_a             Pointer to local node residing in array a (Input)
*/
__global__ void temperature(LB_nodes_gpu n_a, float *cpu_jsquared, int *number_of_non_boundary_nodes ) {
  float mode[4];
  float jsquared = 0.0f;
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    if(!n_a.boundary[index])
    {
      for(int ii=0;ii<LB_COMPONENTS;++ii)
      {  
         calc_mode(mode, n_a, index,ii);
         jsquared = mode[1]*mode[1]+mode[2]*mode[2]+mode[3]*mode[3];
         atomicadd(cpu_jsquared, jsquared);
         atomicAdd(number_of_non_boundary_nodes, 1);
      }
    }
  }
}

/*********************************************************/
/** \name interpolation_two_point_coupling */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param *particle_position    Pointer to the particle position (Input)
 * @param node_index            node index around (8) particle (Output)
 * @param *mode                 Pointer to the 19 modes for current lattice point (Output)
 * @param *d_v                  Pointer to local device values
 * @param *delta                Pointer for the weighting of particle position (Output)
 * @param *interpolated_u       Pointer to the interpolated velocity (Output)
*/
__device__ __inline__ void interpolation_two_point_coupling(LB_nodes_gpu n_a, float *particle_position, unsigned int* node_index, float* mode, LB_rho_v_gpu *d_v, float* delta, float *interpolated_u) {
  int   left_node_index[3];
  float temp_delta[6];
  float temp_delta_half[6];

  // see ahlrichs + duenweg page 8227 equ (10) and (11)
#pragma unroll
  for(int i=0; i<3; ++i)
  {
    float scaledpos = particle_position[i]/para.agrid - 0.5f;
    left_node_index[i] = (int)(floorf(scaledpos));
    temp_delta[3+i] = scaledpos - left_node_index[i];
    temp_delta[i] = 1.0f - temp_delta[3+i];
    // further value used for interpolation of fluid velocity at part pos near boundaries
    temp_delta_half[3+i] = (scaledpos - left_node_index[i])*2.0f;
    temp_delta_half[i] = 2.0f - temp_delta_half[3+i];
  }

  delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
  delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
  delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
  delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
  delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
  delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
  delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
  delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

  // modulo for negative numbers is strange at best, shift to make sure we are positive
  int x = left_node_index[0] + para.dim_x;
  int y = left_node_index[1] + para.dim_y;
  int z = left_node_index[2] + para.dim_z;

  node_index[0] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[1] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[2] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[4] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[5] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[6] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);


  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for(int i=0; i<8; ++i)
  {
    float totmass=0.0f;

    calc_m_from_n(n_a,node_index[i],mode);

#pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;ii++)
    {
      totmass+=mode[0]+para.rho[ii]*para.agrid*para.agrid*para.agrid;
    } 

#ifdef SHANCHEN
    interpolated_u[0] += d_v[node_index[i]].v[0]/8.0f * (n_a.boundary[node_index[i]] == 0);  
    interpolated_u[1] += d_v[node_index[i]].v[1]/8.0f * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += d_v[node_index[i]].v[2]/8.0f * (n_a.boundary[node_index[i]] == 0);
#else
    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary nodes
       to couple with velocity 0 to particles. This is necessary, since boundary nodes
       undergo the same LB dynamics as fluid nodes do. The flow within the boundaries
       does not interact with the physical fluid, since these populations are overwritten
       by the bounce back kernel. Particles close to walls can couple to this unphysical
       flow, though.
    */
    interpolated_u[0] += (mode[1]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] += (mode[2]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += (mode[3]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
#endif
  }
}

/*********************************************************/
/** \name calc_viscous_force */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param partgrad1             particle gradient for the Shan-Chen
 * @param partgrad2             particle gradient for the Shan-Chen
 * @param partgrad3             particle gradient for the Shan-Chen
 * @param *delta                Pointer for the weighting of particle position (Output)
 * @param *delta_j              Pointer for the weighting of particle momentum (Output)
 * @param *particle_data        Pointer to the particle position and velocity (Input)
 * @param *particle_force       Pointer to the particle force (Input)
 * @param *fluid_composition    Pointer to the fluid composition (Input)
 * @param part_index            particle id / thread id (Input)
 * @param *rn_part              Pointer to randomnumber array of the particle
 * @param node_index            node index around (8) particle (Output)
 * @param *d_v                  Pointer to local device values
 * @param flag_cs               Determine if we are at the centre (0, typical) or at the source (1, swimmer only)
*/
__device__ void calc_viscous_force(LB_nodes_gpu n_a, float *delta, float * partgrad1, float * partgrad2, float * partgrad3, CUDA_particle_data *particle_data, float *particle_force, CUDA_fluid_composition * fluid_composition, unsigned int part_index, LB_randomnr_gpu *rn_part, float *delta_j, unsigned int *node_index, LB_rho_v_gpu *d_v, int flag_cs){

  float interpolated_u[3];
  float interpolated_rho[LB_COMPONENTS];
  float viscforce_density[3*LB_COMPONENTS];
  float scforce_density[3*LB_COMPONENTS];
  float mode[19*LB_COMPONENTS];
#ifdef SHANCHEN
  float gradrho1, gradrho2, gradrho3;
  float Rho;
#endif
  // Zero out workspace
  #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    #pragma unroll
    for(int jj=0; jj<3; ++jj)
    { 
      scforce_density[jj+ii*3]  =0.0f;
      viscforce_density[jj+ii*3]=0.0f;
      delta_j[jj+ii*3]  =0.0f;
    }
    
    #pragma unroll
    for(int jj=0; jj<8; ++jj)
    { 
      partgrad1[jj+ii*8]=0.0f;
      partgrad2[jj+ii*8]=0.0f;
      partgrad3[jj+ii*8]=0.0f;
    }
  }
  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3*part_index+0] = flag_cs * particle_force[3*part_index+0];
  particle_force[3*part_index+1] = flag_cs * particle_force[3*part_index+1];
  particle_force[3*part_index+2] = flag_cs * particle_force[3*part_index+2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) * particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the particle
  position[0] += flag_cs * direction * particle_data[part_index].swim.quatu[0];
  position[1] += flag_cs * direction * particle_data[part_index].swim.quatu[1];
  position[2] += flag_cs * direction * particle_data[part_index].swim.quatu[2];
#endif

  // Do the velocity interpolation
  interpolation_two_point_coupling(n_a, position, node_index, mode, d_v, delta, interpolated_u);

#ifdef ENGINE
  velocity[0] -= particle_data[part_index].swim.v_swim*particle_data[part_index].swim.quatu[0];
  velocity[1] -= particle_data[part_index].swim.v_swim*particle_data[part_index].swim.quatu[1];
  velocity[2] -= particle_data[part_index].swim.v_swim*particle_data[part_index].swim.quatu[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD units
  particle_data[part_index].swim.v_cs[0+3*flag_cs] = interpolated_u[0] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[1+3*flag_cs] = interpolated_u[1] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[2+3*flag_cs] = interpolated_u[2] * para.agrid / para.tau;
#endif

#ifdef SHANCHEN

 #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    float solvation2 = particle_data[part_index].solvation[2*ii + 1];
   
    interpolated_rho[ii]  = 0.0f;
    gradrho1 = gradrho2 = gradrho3 = 0.0f;
  
    // TODO: should one introduce a density-dependent friction ?
    calc_mode(mode, n_a, node_index[0],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[0] * Rho; 
    partgrad1[ii*8 + 0] += Rho * solvation2;
    partgrad2[ii*8 + 0] += Rho * solvation2;
    partgrad3[ii*8 + 0] += Rho * solvation2;
    gradrho1 -=(delta[0] + delta[1]) * Rho; 
    gradrho2 -=(delta[0] + delta[2]) * Rho; 
    gradrho3 -=(delta[0] + delta[4]) * Rho; 

    calc_mode(mode, n_a, node_index[1],ii); 
    Rho = mode[0] +  para.rho[ii]*para.agrid*para.agrid*para.agrid; 
    interpolated_rho[ii] += delta[1] * Rho; 
    partgrad1[ii*8 + 1] -= Rho * solvation2;
    partgrad2[ii*8 + 1] += Rho * solvation2;
    partgrad3[ii*8 + 1] += Rho * solvation2;
    gradrho1 +=(delta[1] + delta[0]) * Rho; 
    gradrho2 -=(delta[1] + delta[3]) * Rho; 
    gradrho3 -=(delta[1] + delta[5]) * Rho; 

    calc_mode(mode, n_a, node_index[2],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[2] * Rho; 
    partgrad1[ii*8 + 2] += Rho * solvation2;
    partgrad2[ii*8 + 2] -= Rho * solvation2;
    partgrad3[ii*8 + 2] += Rho * solvation2;
    gradrho1 -=(delta[2] + delta[3]) * Rho; 
    gradrho2 +=(delta[2] + delta[0]) * Rho; 
    gradrho3 -=(delta[2] + delta[6]) * Rho; 

    calc_mode(mode, n_a, node_index[3],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[3] * Rho; 
    partgrad1[ii*8 + 3] -= Rho * solvation2;
    partgrad2[ii*8 + 3] -= Rho * solvation2;
    partgrad3[ii*8 + 3] += Rho * solvation2;
    gradrho1 +=(delta[3] + delta[2]) * Rho; 
    gradrho2 +=(delta[3] + delta[1]) * Rho; 
    gradrho3 -=(delta[3] + delta[7]) * Rho; 

    calc_mode(mode, n_a, node_index[4],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[4] * Rho; 
    partgrad1[ii*8 + 4] += Rho * solvation2;
    partgrad2[ii*8 + 4] += Rho * solvation2;
    partgrad3[ii*8 + 4] -= Rho * solvation2;
    gradrho1 -=(delta[4] + delta[5]) * Rho; 
    gradrho2 -=(delta[4] + delta[6]) * Rho; 
    gradrho3 +=(delta[4] + delta[0]) * Rho; 

    calc_mode(mode, n_a, node_index[5],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[5] * Rho; 
    partgrad1[ii*8 + 5] -= Rho * solvation2;
    partgrad2[ii*8 + 5] += Rho * solvation2;
    partgrad3[ii*8 + 5] -= Rho * solvation2;
    gradrho1 +=(delta[5] + delta[4]) * Rho; 
    gradrho2 -=(delta[5] + delta[7]) * Rho; 
    gradrho3 +=(delta[5] + delta[1]) * Rho; 

    calc_mode(mode, n_a, node_index[6],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[6] * Rho; 
    partgrad1[ii*8 + 6] += Rho * solvation2;
    partgrad2[ii*8 + 6] -= Rho * solvation2;
    partgrad3[ii*8 + 6] -= Rho * solvation2;
    gradrho1 -=(delta[6] + delta[7]) * Rho; 
    gradrho2 +=(delta[6] + delta[4]) * Rho; 
    gradrho3 +=(delta[6] + delta[2]) * Rho; 

    calc_mode(mode, n_a, node_index[7],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[7] * Rho; 
    partgrad1[ii*8 + 7] -= Rho * solvation2;
    partgrad2[ii*8 + 7] -= Rho * solvation2;
    partgrad3[ii*8 + 7] -= Rho * solvation2;
    gradrho1 +=(delta[7] + delta[6]) * Rho; 
    gradrho2 +=(delta[7] + delta[5]) * Rho; 
    gradrho3 +=(delta[7] + delta[3]) * Rho; 

    /* normalize the gradient to md units TODO: is that correct?*/
    gradrho1 *= para.agrid; 
    gradrho2 *= para.agrid; 
    gradrho3 *= para.agrid; 

    // scforce_density is 0 at the interpolated point where the swimming force gets put back on the fluid
    // only add sc_force up for particle centre <=> (1-flag_cs) = 1
    scforce_density[0+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho1 ; 
    scforce_density[1+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho2 ;
    scforce_density[2+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho3 ;

    /* scforce_density is used also later...*/
    particle_force[3*part_index+0] += scforce_density[0+ii*3];
    particle_force[3*part_index+1] += scforce_density[1+ii*3];
    particle_force[3*part_index+2] += scforce_density[2+ii*3];
    // only set fluid_composition for particle centre <=> (1-flag_cs) = 1
    fluid_composition[part_index].weight[ii] = (1-flag_cs) * interpolated_rho[ii];
 }

#else // SHANCHEN is not defined

  /* for LB we do not reweight the friction force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){
    interpolated_rho[ii]=1.0;
  }

#endif // SHANCHEN

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot=0;

  #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){
    rhotot+=interpolated_rho[ii];
  }

  /* Viscous force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    viscforce_density[0+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[0] - interpolated_u[0]*para.agrid/para.tau)/rhotot;
    viscforce_density[1+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[1] - interpolated_u[1]*para.agrid/para.tau)/rhotot;
    viscforce_density[2+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[2] - interpolated_u[2]*para.agrid/para.tau)/rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
    viscforce_density[0+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[0]/rhotot;
    viscforce_density[1+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[1]/rhotot;
    viscforce_density[2+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[2]/rhotot;
#endif

    /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
#ifdef FLATNOISE
    random_01(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
    viscforce_density[1+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[1]-0.5f);
    random_01(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
#elif defined(GAUSSRANDOMCUT)
    gaussian_random_cut(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce_density[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random_cut(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#elif defined(GAUSSRANDOM)
    gaussian_random(rn_part);
    viscforce_density[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce_density[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random(rn_part);
    viscforce_density[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#else
#error No noise type defined for the GPU LB
#endif 

    /** delta_j for transform momentum transfer to lattice units which is done in calc_node_force
      (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */

    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    particle_force[3*part_index+0] += (1-flag_cs) * viscforce_density[0+ii*3];
    particle_force[3*part_index+1] += (1-flag_cs) * viscforce_density[1+ii*3];
    particle_force[3*part_index+2] += (1-flag_cs) * viscforce_density[2+ii*3];

    /* the average force from the particle to surrounding nodes is transmitted back to preserve momentum */
    for(int node=0 ; node < 8 ; node++ ) { 
      particle_force[3*part_index+0] -= (1-flag_cs) * partgrad1[node+ii*8]/8.0f;
      particle_force[3*part_index+1] -= (1-flag_cs) * partgrad2[node+ii*8]/8.0f;
      particle_force[3*part_index+2] -= (1-flag_cs) * partgrad3[node+ii*8]/8.0f;
    }

    /* note that scforce_density is zero if SHANCHEN is not #defined */
    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    delta_j[0+3*ii] -= (scforce_density[0+ii*3]+(1-flag_cs)*viscforce_density[0+ii*3])*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= (scforce_density[1+ii*3]+(1-flag_cs)*viscforce_density[1+ii*3])*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= (scforce_density[2+ii*3]+(1-flag_cs)*viscforce_density[2+ii*3])*para.time_step*para.tau/para.agrid;

#ifdef ENGINE
    // add swimming force to source position
    delta_j[0+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[0]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[1]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[2]*para.time_step*para.tau/para.agrid;
#endif

  }

#ifdef SHANCHEN
  for(int node=0 ; node < 8 ; node++ ) { 
    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ ) { 
      partgrad1[node+ii*8]*=(para.time_step*para.tau/para.agrid);
      partgrad2[node+ii*8]*=(para.time_step*para.tau/para.agrid);
      partgrad3[node+ii*8]*=(para.time_step*para.tau/para.agrid);
    }
  }
#endif

}



/**calculation of the node force caused by the particles, with atomicadd due to avoiding race conditions 
  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param *delta        Pointer for the weighting of particle position (Input)
 * @param partgrad1             particle gradient for the Shan-Chen
 * @param partgrad2             particle gradient for the Shan-Chen
 * @param partgrad3             particle gradient for the Shan-Chen
 * @param *delta_j      Pointer for the weighting of particle momentum (Input)
 * @param node_index    node index around (8) particle (Input)
 * @param node_f        Pointer to the node force (Output).
*/
__device__ void calc_node_force(float *delta, float *delta_j, float * partgrad1, float * partgrad2, float * partgrad3,  unsigned int *node_index, LB_node_force_density_gpu node_f){
/* TODO: should the drag depend on the density?? */
/* NOTE: partgrad is not zero only if SHANCHEN is defined. It is initialized in calc_node_force. Alternatively one could 
         specialize this function to the single component LB */ 
  for(int ii=0; ii < LB_COMPONENTS; ++ii)
  { 
    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[0+ii*3] + partgrad1[ii*8+0]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[1+ii*3] + partgrad2[ii*8+0]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[2+ii*3] + partgrad3[ii*8+0]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[0+ii*3] + partgrad1[ii*8+1]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[1+ii*3] + partgrad2[ii*8+1]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[2+ii*3] + partgrad3[ii*8+1]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[0+ii*3] + partgrad1[ii*8+2]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[1+ii*3] + partgrad2[ii*8+2]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[2+ii*3] + partgrad3[ii*8+2]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[0+ii*3] + partgrad1[ii*8+3]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[1+ii*3] + partgrad2[ii*8+3]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[2+ii*3] + partgrad3[ii*8+3]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[0+ii*3] + partgrad1[ii*8+4]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[1+ii*3] + partgrad2[ii*8+4]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[2+ii*3] + partgrad3[ii*8+4]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[0+ii*3] + partgrad1[ii*8+5]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[1+ii*3] + partgrad2[ii*8+5]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[2+ii*3] + partgrad3[ii*8+5]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[0+ii*3] + partgrad1[ii*8+6]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[1+ii*3] + partgrad2[ii*8+6]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[2+ii*3] + partgrad3[ii*8+6]));

    atomicadd(&(node_f.force_density[(0+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[0+ii*3] + partgrad1[ii*8+7]));
    atomicadd(&(node_f.force_density[(1+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[1+ii*3] + partgrad2[ii*8+7]));
    atomicadd(&(node_f.force_density[(2+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[2+ii*3] + partgrad3[ii*8+7]));
  }
}

/*********************************************************/
/** \name System setup and Kernel functions */
/*********************************************************/

/**kernel to calculate local populations from hydrodynamic fields.
 * The mapping is given in terms of the equilibrium distribution.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 * Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 * @param n_a   Pointer to the lattice site (Input).
 * @param *gpu_check additional check if gpu kernel are executed(Input).
 * @param *d_v    Pointer to local device values
 * @param *node_f          Pointer to node forces
*/
__global__ void calc_n_from_rho_j_pi(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f, int *gpu_check) {
   /* TODO: this can handle only a uniform density, something similar, but local, 
            has to be called every time the fields are set by the user ! */ 
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(index<para.number_of_nodes)
  {
    float mode[19*LB_COMPONENTS];

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      /** default values for fields in lattice units */
      gpu_check[0] = 1;
     
      float Rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;
      float v[3] = { 0.0f, 0.0f, 0.0f };
      float pi[6] = { Rho*c_sound_sq, 0.0f, Rho*c_sound_sq, 0.0f, 0.0f, Rho*c_sound_sq };
     
      float rhoc_sq = Rho*c_sound_sq;
      float avg_rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;
      float local_rho, local_j[3], *local_pi, trace;
     
      local_rho  = Rho;
     
      local_j[0] = Rho * v[0];
      local_j[1] = Rho * v[1];
      local_j[2] = Rho * v[2];
     
      local_pi = pi;
     
      /** reduce the pressure tensor to the part needed here. 
          NOTE: this not true anymore for SHANCHEN 
          if the densities are not uniform. FIXME*/

      local_pi[0] -= rhoc_sq;
      local_pi[2] -= rhoc_sq;
      local_pi[5] -= rhoc_sq;
     
      trace = local_pi[0] + local_pi[2] + local_pi[5];
     
      float rho_times_coeff;
      float tmp1,tmp2;
     
      /** update the q=0 sublattice */
      n_a.vd[(0 + ii*LBQ ) * para.number_of_nodes + index] = 1.0f/3.0f * (local_rho-avg_rho) - 1.0f/2.0f*trace;
     
      /** update the q=1 sublattice */
      rho_times_coeff = 1.0f/18.0f * (local_rho-avg_rho);
     
      n_a.vd[(1 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[0] + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(2 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[0] + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(3 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[1] + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(4 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[1] + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(5 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[2] + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
      n_a.vd[(6 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[2] + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
     
      /** update the q=2 sublattice */
      rho_times_coeff = 1.0f/36.0f * (local_rho-avg_rho);
     
      tmp1 = local_pi[0] + local_pi[2];
      tmp2 = 2.0f*local_pi[1];
      n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[1]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[1]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[1]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[1]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      tmp1 = local_pi[0] + local_pi[5];
      tmp2 = 2.0f*local_pi[3];
     
      n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      tmp1 = local_pi[2] + local_pi[5];
      tmp2 = 2.0f*local_pi[4];
     
      n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[1]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[1]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[1]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[1]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      /**set different seed for randomgen on every node */
      n_a.seed[index] = para.your_seed + index;
    }

    calc_m_from_n(n_a,index,mode);
    update_rho_v(mode,index,node_f,d_v);
  }
}

/** kernel to calculate local populations from hydrodynamic fields
 * from given flow field velocities. The mapping is given in terms of
 * the equilibrium distribution.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 * Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 * @param n_a              the current nodes array (double buffering!)
 * @param single_nodeindex the node to set the velocity for
 * @param velocity         the velocity to set
 * @param *d_v             Pointer to local device values (Input)
 * @param *node_f          Pointer to node forces (Input)
 */ 
__global__ void set_u_from_rho_v_pi( LB_nodes_gpu n_a, int single_nodeindex, float *velocity, LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f ) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
  {
    float local_rho;
    float local_j[3];
    float local_pi[6];
    float trace, avg_rho;
    float rho_times_coeff;
    float tmp1, tmp2; 

    float mode_for_pi[19*LB_COMPONENTS];
    float rho_from_m[1*LB_COMPONENTS];
    float j_from_m[3*LB_COMPONENTS];
    float pi_from_m[6*LB_COMPONENTS];

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Reset the d_v

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);

    // Calculate the density, velocity, and pressure tensor
    // in LB unit for this node

    calc_values_from_m_in_LB_units( mode_for_pi, &d_v[single_nodeindex], rho_from_m, j_from_m, pi_from_m);

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      // Take LB component density and calculate the equilibrium part

      local_rho = rho_from_m[ii];
      avg_rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;

      // Take LB component velocity and make it a momentum

      local_j[0] = local_rho * velocity[0];
      local_j[1] = local_rho * velocity[1];
      local_j[2] = local_rho * velocity[2];

      // Take LB component pressure tensor and put in equilibrium

      local_pi[0] = pi_from_m[6*ii + 0];
      local_pi[1] = pi_from_m[6*ii + 1];
      local_pi[2] = pi_from_m[6*ii + 2];
      local_pi[3] = pi_from_m[6*ii + 3];
      local_pi[4] = pi_from_m[6*ii + 4];
      local_pi[5] = pi_from_m[6*ii + 5];

      trace = local_pi[0] + local_pi[2] + local_pi[5];

      // update the q=0 sublattice

      n_a.vd[(0 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/3.0f * (local_rho - avg_rho) - 1.0f/2.0f*trace;

      // update the q=1 sublattice

      rho_times_coeff = 1.0f/18.0f * (local_rho - avg_rho);

      n_a.vd[(1 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[0]
                                                                        + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(2 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[0]
                                                                        + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(3 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[1]
                                                                        + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(4 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[1]
                                                                        + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(5 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[2]
                                                                        + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
      n_a.vd[(6 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[2]
                                                                        + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;

      // update the q=2 sublattice

      rho_times_coeff = 1.0f/36.0f * (local_rho - avg_rho);

      tmp1 = local_pi[0] + local_pi[2];
      tmp2 = 2.0f*local_pi[1];

      n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[1])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[1])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[1])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[1])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;

      tmp1 = local_pi[0] + local_pi[5];
      tmp2 = 2.0f*local_pi[3];

      n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;

      tmp1 = local_pi[2] + local_pi[5];
      tmp2 = 2.0f*local_pi[4];

      n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[1]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[1]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[1]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[1]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
    }

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Update the density and velocity field for this mode

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);
  }
}



/**calculate mass of the whole fluid kernel
 * @param *sum    Pointer to result storage value (Output)
 * @param n_a     Pointer to local node residing in array a (Input)
*/
__global__ void calc_mass(LB_nodes_gpu n_a, float *sum) {
  float mode[4];

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      calc_mode(mode, n_a, index,ii);
      float Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      atomicadd(&(sum[0]), Rho);
    }
  }
}

/** (re-)initialization of the node force density / set up of external force density in lb units
 * @param node_f  Pointer to local node force density (Input)
*/
__global__ void reinit_node_force(LB_node_force_density_gpu node_f){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    {
      node_f.force_density[(0+ii*3)*para.number_of_nodes + index] = 0.0f;
      node_f.force_density[(1+ii*3)*para.number_of_nodes + index] = 0.0f;
      node_f.force_density[(2+ii*3)*para.number_of_nodes + index] = 0.0f;
    }
  }
}


/**set extern force on single nodes kernel
 * @param n_extern_node_force_densities   number of nodes (Input)
 * @param *extern_node_force_densities    Pointer to extern node force array (Input)
 * @param node_f                node force struct (Output)
*/
__global__ void init_extern_node_force_densities(int n_extern_node_force_densities, LB_extern_nodeforcedensity_gpu *extern_node_force_densities, LB_node_force_density_gpu node_f){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  float factor=powf(para.agrid,2)*para.tau*para.tau;
  if(index<n_extern_node_force_densities)
  {
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    {
      node_f.force_density[(0+ii*3)*para.number_of_nodes + extern_node_force_densities[index].index] = extern_node_force_densities[index].force_density[0] * factor;
      node_f.force_density[(1+ii*3)*para.number_of_nodes + extern_node_force_densities[index].index] = extern_node_force_densities[index].force_density[1] * factor;
      node_f.force_density[(2+ii*3)*para.number_of_nodes + extern_node_force_densities[index].index] = extern_node_force_densities[index].force_density[2] * factor;
    }
  }
}

#ifdef SHANCHEN

/** 
 * @param single_nodeindex  Single node index        (Input)
 * @param component_index   Shanchen component index        (Input)
 * @param n_a               Pointer to local node residing in array a(Input)
*/
__device__ __inline__ float calc_massmode(LB_nodes_gpu n_a, int single_nodeindex, int component_index){

  /** mass mode */
  float mode;
  mode =   n_a.vd[( 0 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] 
         + n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex];

 mode += para.rho[component_index]*para.agrid*para.agrid*para.agrid;

 return mode;
}

__device__ __inline__ void calc_shanchen_contribution(LB_nodes_gpu n_a,int component_index, int x, int y, int z, float *p){ 

  float tmp_p[3]={0.0f,0.0f,0.0f};
  float pseudo;
  int index;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/18.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/18.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/18.0f;

  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/18.0f;

  index  = x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[2]+=pseudo/18.0f;

  index  = x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[2]-=pseudo/18.0f;

  index  = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[1]+=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[1]-=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[1]-=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[1]+=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;
  
  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;
 
  p[0]=tmp_p[0];
  p[1]=tmp_p[1];
  p[2]=tmp_p[2];
}

/** function to calc shanchen forces 
 * @param n_a     Pointer to local node residing in array a(Input)
 * @param node_f  Pointer to local node force (Input)
*/
__global__ void lb_shanchen_GPU(LB_nodes_gpu n_a,LB_node_force_density_gpu node_f){
#ifndef D3Q19
#error Lattices other than D3Q19 not supported
#endif
#if ( LB_COMPONENTS == 1  ) 
  #warning shanchen forces not implemented 
#else  

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int xyz[3];
  float pseudo;

  if(index<para.number_of_nodes)
  if( n_a.boundary[index] == 0 )
  {

    /* ShanChen forces are not reset at the end of the integration cycle, 
       in order to compute properly the hydrodynamic fields, so we have
       to reset them here. For the standard LB this is not needed */
     reset_LB_force_densities(index, node_f) ;
     /*Let's first identify the neighboring nodes */
     index_to_xyz(index, xyz);
     int x = xyz[0];
     int y = xyz[1];
     int z = xyz[2];
     
     #pragma unroll
     for(int ii=0;ii<LB_COMPONENTS;ii++)
     { 
       float p[3]={0.0f,0.0f,0.0f};
       pseudo =  calc_massmode(n_a,index,ii);

       #pragma unroll
       for(int jj=0;jj<LB_COMPONENTS;jj++)
       { 
             float tmpp[3]={0.0f,0.0f,0.0f};
             calc_shanchen_contribution(n_a, jj, x,y,z, tmpp);

             // FIXME  coupling HAS to be rescaled with agrid....
             p[0] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[0];
             p[1] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[1];
             p[2] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[2];
       }

       node_f.force_density[(0+ii*3)*para.number_of_nodes + index]+=p[0];
       node_f.force_density[(1+ii*3)*para.number_of_nodes + index]+=p[1];
       node_f.force_density[(2+ii*3)*para.number_of_nodes + index]+=p[2];
/* copy to be used when resetting force densities */
       node_f.scforce_density[(0+ii*3)*para.number_of_nodes + index]=p[0];
       node_f.scforce_density[(1+ii*3)*para.number_of_nodes + index]=p[1];
       node_f.scforce_density[(2+ii*3)*para.number_of_nodes + index]=p[2];
    }
  }
#endif 
  return; 
}

#endif //SHANCHEN

/** kernel to set the local density
 *
 * @param n_a               the current nodes array (double buffering!)
 * @param single_nodeindex  the node to set the velocity for
 * @param rho               the density to set
 * @param d_v                Pointer to the local modes
*/
__global__ void set_rho(LB_nodes_gpu n_a,  LB_rho_v_gpu *d_v, int single_nodeindex,float *rho) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  /*Note: this sets the velocities to zero */
  if(index == 0)
  {
    float local_rho;

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
       /** default values for fields in lattice units */
       local_rho = (rho[ii]-para.rho[ii])*para.agrid*para.agrid*para.agrid;
       d_v[single_nodeindex].rho[ii]=rho[ii];

       n_a.vd[(0  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/ 3.0f * local_rho ;
       n_a.vd[(1  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(2  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(3  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(4  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(5  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(6  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(7  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(8  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(9  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
    }
  }
}

/**set the boundary flag for all boundary nodes
 * @param boundary_node_list    The indices of the boundary nodes
 * @param boundary_index_list   The flag representing the corresponding boundary
 * @param number_of_boundnodes  The number of boundary nodes
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param n_b                   Pointer to local node residing in array b (Input)
*/
__global__ void init_boundaries(int *boundary_node_list, int *boundary_index_list, int number_of_boundnodes, LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<number_of_boundnodes)
  {
    n_a.boundary[boundary_node_list[index]] = boundary_index_list[index];
    n_b.boundary[boundary_node_list[index]] = boundary_index_list[index];
  }
}

/**reset the boundary flag of every node
 * @param n_a   Pointer to local node residing in array a (Input)
 * @param n_b   Pointer to local node residing in array b (Input)
*/
__global__ void reset_boundaries(LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  size_t index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
    n_a.boundary[index] = n_b.boundary[index] = 0;
}

/** integrationstep of the lb-fluid-solver
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param n_b     Pointer to local node residing in array b (Input)
 * @param *d_v    Pointer to local device values (Input)
 * @param node_f  Pointer to local node force density (Input)
 * @param ek_parameters_gpu  Pointer to the parameters for the electrokinetics (Input)
*/


__global__ void integrate(LB_nodes_gpu n_a, LB_nodes_gpu n_b, LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f, EK_parameters* ek_parameters_gpu) {
  /**every node is connected to a thread via the index*/
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  /**the 19 moments (modes) are only temporary register values */
  float mode[19*LB_COMPONENTS];
  LB_randomnr_gpu rng;

  if( index < para.number_of_nodes )
  {
    /** storing the seed into a register value*/
    rng.seed = n_a.seed[index];
    /**calc_m_from_n*/
    calc_m_from_n(n_a, index, mode);
    /**lb_relax_modes*/
    relax_modes(mode, index, node_f,d_v);
    /**lb_thermalize_modes */
    if (para.fluct)
    {
      thermalize_modes(mode, index, &rng);
    }
#if  defined(EXTERNAL_FORCES)  ||   defined (SHANCHEN)  
    /**if external force is used apply node force */
    apply_forces(index, mode, node_f,d_v);
#endif
    /**lb_calc_n_from_modes_push*/
    normalize_modes(mode);
    /**calc of velocity densities and streaming with pbc*/
    calc_n_from_modes_push(n_b, mode, index);
    /** rewriting the seed back to the global memory*/
    n_b.seed[index] = rng.seed;
  }  
}

/** part interaction kernel
 * @param n_a                Pointer to local node residing in array a (Input)
 * @param *particle_data     Pointer to the particle position and velocity (Input)
 * @param *particle_force    Pointer to the particle force (Input)
 * @param *part              Pointer to the rn array of the particles (Input)
 * @param node_f             Pointer to local node force (Input)
 * @param *fluid_composition Pointer to the local fluid composition for the Shanchen
 * @param *d_v               Pointer to local device values
*/
__global__ void calc_fluid_particle_ia(LB_nodes_gpu n_a, CUDA_particle_data *particle_data, float *particle_force, CUDA_fluid_composition * fluid_composition, LB_node_force_density_gpu node_f, CUDA_particle_seed *part, LB_rho_v_gpu *d_v, bool couple_virtual){

  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[8];
  float delta[8];
  float delta_j[3*LB_COMPONENTS]; 
  float partgrad1[8*LB_COMPONENTS]; 
  float partgrad2[8*LB_COMPONENTS]; 
  float partgrad3[8*LB_COMPONENTS]; 
  LB_randomnr_gpu rng_part;
  if(part_index<para.number_of_particles)
  {
#if defined(VIRTUAL_SITES)
    if ( !particle_data[part_index].is_virtual || couple_virtual )
#endif
    {
      rng_part.seed = part[part_index].seed;

      /**force acting on the particle. delta_j will be used later to compute the force that acts back onto the fluid. */
      calc_viscous_force(n_a, delta, partgrad1, partgrad2, partgrad3, particle_data, particle_force, fluid_composition,part_index, &rng_part, delta_j, node_index, d_v, 0);
      calc_node_force(delta, delta_j, partgrad1, partgrad2, partgrad3, node_index, node_f); 

#ifdef ENGINE
      if ( particle_data[part_index].swim.swimming ) {
        calc_viscous_force(n_a, delta, partgrad1, partgrad2, partgrad3, particle_data, particle_force, fluid_composition,part_index, &rng_part, delta_j, node_index, d_v, 1);
        calc_node_force(delta, delta_j, partgrad1, partgrad2, partgrad3, node_index, node_f);
      }
#endif

      /**force which acts back to the fluid node */
      part[part_index].seed = rng_part.seed;
    }
  }
}

/** part interaction kernel
 * @param n_a       Pointer to local node residing in array a (Input)
 * @param *particle_data    Pointer to the particle position and velocity (Input)
 * @param *particle_force   Pointer to the particle force (Input)
 * @param *part       Pointer to the rn array of the particles (Input)
 * @param node_f      Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__global__ void calc_fluid_particle_ia_three_point_couple(LB_nodes_gpu n_a, CUDA_particle_data *particle_data, float *particle_force, LB_node_force_density_gpu node_f, CUDA_particle_seed *part, LB_rho_v_gpu *d_v){

  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[27];
  float delta[27];
  float delta_j[3*LB_COMPONENTS]; 
  LB_randomnr_gpu rng_part;
  if(part_index<para.number_of_particles){

    rng_part.seed = part[part_index].seed;
    /**force acting on the particle. delta_j will be used later to compute the force that acts back onto the fluid. */
    calc_viscous_force_three_point_couple(n_a, delta, particle_data, particle_force, part_index, &rng_part, delta_j, node_index,d_v,0);
    calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);

#ifdef ENGINE
    if ( particle_data[part_index].swim.swimming ) {
      calc_viscous_force_three_point_couple(n_a, delta, particle_data, particle_force, part_index, &rng_part, delta_j, node_index,d_v,1);
      calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);
    }
#endif

    /**force which acts back to the fluid node */
    part[part_index].seed = rng_part.seed;    
  }
}


#ifdef LB_BOUNDARIES_GPU
/**Bounce back boundary kernel
 * @param n_a         Pointer to local node residing in array a (Input)
 * @param n_b         Pointer to local node residing in array b (Input)
 * @param lb_boundary_velocity    The constant velocity at the boundary, set by the user (Input)
 * @param lb_boundary_force       The force on the boundary nodes (Output)
*/
__global__ void apply_boundaries(LB_nodes_gpu n_curr, float* lb_boundary_velocity, float* lb_boundary_force){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
    bounce_back_boundaries(n_curr, index, lb_boundary_velocity, lb_boundary_force);
}
#ifdef SHANCHEN
__global__ void lb_shanchen_set_boundaries(LB_nodes_gpu n_curr){
/* This implements neutral boundary conditions for the shanchen fluid (i.e., 90 deg contact angle) */

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int xyz[3];
  if(index<para.number_of_nodes){
   if(n_curr.boundary[index] != 0 ) { 
    index_to_xyz(index, xyz);
    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];
    unsigned int to_index_x,to_index_y,to_index_z,to_index;
    int c[3],count=0;

    for(int ii=0;ii<LB_COMPONENTS;ii++)
       for(int comp=0;comp<19;comp++)
          n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index]  =  0.0 ; 
    for(c[0]=-1;c[0]<=1;c[0]++){
       for(c[1]=-1;c[1]<=1;c[1]++){
          for(c[2]=-1;c[2]<=1;c[2]++){
             to_index_x = (x+c[0]+para.dim_x)%para.dim_x; 
             to_index_y = (y+c[1]+para.dim_y)%para.dim_y; 
             to_index_z = (z+c[2]+para.dim_z)%para.dim_z; 
             to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z;  
	     if(n_curr.boundary[to_index] == 0 ) { 
                  for(int ii=0;ii<LB_COMPONENTS;ii++){
                     for(int comp=0;comp<19;comp++){ /* We copy all velocities: at the end we will need 
                                                        only the density mode, but this introduces no overhead anyway */
	                 n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index] += 
	                    n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + to_index] ;
                         count++;
                     }
                  }
             }	      	
          }
       }
    }
    if(count>0)
      for(int ii=0;ii<LB_COMPONENTS;ii++)
        for(int comp=0;comp<19;comp++)
           n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index]  /= count ; 
  }
 }
}
#endif /* SHANCHEN */


#endif

/** get physical values of the nodes (density, velocity, ...)
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param *p_v    Pointer to local print values (Output)
 * @param *d_v    Pointer to local device values (Input)
 * @param node_f  The forces on the LB nodes
*/
__global__ void get_mesoscopic_values_in_MD_units(LB_nodes_gpu n_a, LB_rho_v_pi_gpu *p_v,LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index < para.number_of_nodes)
  {
    float mode[19*LB_COMPONENTS];
    calc_m_from_n(n_a, index, mode);
    calc_values_in_MD_units(n_a, mode, p_v, d_v, node_f, index, index);
  }
}

/** get boundary flags
 *  @param n_a                Pointer to local node residing in array a (Input)
 *  @param device_bound_array Pointer to local device values (Input)
 */
__global__ void lb_get_boundaries(LB_nodes_gpu n_a, unsigned int *device_bound_array){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
   device_bound_array[index] = n_a.boundary[index];
}



/**print single node values kernel
 * @param single_nodeindex  index of the node (Input)
 * @param *d_p_v            Pointer to result storage array (Input)
 * @param n_a               Pointer to local node residing in array a (Input)
 * @param *d_v    Pointer to local device values
 * @param node_f  Pointer to local node force
*/
__global__ void lb_print_node(int single_nodeindex, LB_rho_v_pi_gpu *d_p_v, LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_density_gpu node_f){

  float mode[19*LB_COMPONENTS];
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
  {
    calc_m_from_n(n_a, single_nodeindex, mode);
     
    /* the following actually copies rho and v from d_v, and calculates pi */
    calc_values_in_MD_units(n_a, mode, d_p_v, d_v, node_f, single_nodeindex, 0);
  }
}

__global__ void momentum(LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_density_gpu node_f, float *sum) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    float j[3]={0.0f,0.0f,0.0f};
    float mode[4];

    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ )
    { 
      calc_mode(mode, n_a, index,ii);

      j[0] += mode[1]+node_f.force_density[(0+ii*3)*para.number_of_nodes + index];
      j[1] += mode[2]+node_f.force_density[(1+ii*3)*para.number_of_nodes + index];
      j[2] += mode[3]+node_f.force_density[(2+ii*3)*para.number_of_nodes + index];
    }

#ifdef LB_BOUNDARIES_GPU
    if(n_a.boundary[index])
      j[0]=j[1]=j[2]=0.0f;
#endif

    atomicadd(&(sum[0]), j[0]); 
    atomicadd(&(sum[1]), j[1]); 
    atomicadd(&(sum[2]), j[2]); 
  }
}
__global__ void remove_momentum(LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_density_gpu node_f, float *sum) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(index<para.number_of_nodes){
    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ ) { 
        node_f.force_density[(0+ii*3)*para.number_of_nodes + index]-=sum[0]/para.number_of_nodes;
        node_f.force_density[(1+ii*3)*para.number_of_nodes + index]-=sum[1]/para.number_of_nodes;
        node_f.force_density[(2+ii*3)*para.number_of_nodes + index]-=sum[2]/para.number_of_nodes;
    }
  }
}

/**print single node boundary flag
 * @param single_nodeindex  index of the node (Input)
 * @param *device_flag      Pointer to result storage array (Input)
 * @param n_a               Pointer to local node residing in array a (Input)
*/
__global__ void lb_get_boundary_flag(int single_nodeindex, unsigned int *device_flag, LB_nodes_gpu n_a){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
    device_flag[0] = n_a.boundary[single_nodeindex];
}

/**********************************************************************/
/* Host functions to setup and call kernels*/
/**********************************************************************/

void lb_get_para_pointer(LB_parameters_gpu** pointeradress) {
  if(cudaGetSymbolAddress((void**) pointeradress, para) != cudaSuccess)
  {
    fprintf(stderr, "Trouble getting address of LB parameters.\n"); //TODO give proper error message
    errexit();
  }
}

void lb_get_lbpar_pointer(LB_parameters_gpu** pointeradress) {
  *pointeradress = &lbpar_gpu;
}


void lb_get_boundary_force_pointer(float** pointeradress) {
#ifdef LB_BOUNDARIES_GPU
  *pointeradress = lb_boundary_force;
#endif
}

void lb_get_device_values_pointer(LB_rho_v_gpu** pointeradress) {
  *pointeradress = device_rho_v;
}

/**initialization for the lb gpu fluid called from host
 * @param *lbpar_gpu  Pointer to parameters to setup the lb field
*/
void lb_init_GPU(LB_parameters_gpu *lbpar_gpu){
#define free_realloc_and_clear(var,size)\
  { if( (var) != nullptr ) cuda_safe_mem(cudaFree((var))); cuda_safe_mem(cudaMalloc((void**)&var, size)); cudaMemset(var,0,size); } 

  size_of_rho_v     = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_gpu);
  size_of_rho_v_pi  = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_pi_gpu);


  /** Allocate structs in device memory*/
  /* see the notes to the stucture device_rho_v_pi above...*/
  if(extended_values_flag==0) 
  {
    free_realloc_and_clear(device_rho_v, size_of_rho_v);
  }
  else 
  {
    free_realloc_and_clear(device_rho_v_pi, size_of_rho_v_pi);
  }

  /* TODO: this is a almost a copy copy of  device_rho_v think about eliminating it, and maybe pi can be added to device_rho_v in this case*/
  free_realloc_and_clear(print_rho_v_pi  , size_of_rho_v_pi);
  free_realloc_and_clear(nodes_a.vd      , lbpar_gpu->number_of_nodes * 19 * LB_COMPONENTS * sizeof(float));
  free_realloc_and_clear(nodes_b.vd      , lbpar_gpu->number_of_nodes * 19 * LB_COMPONENTS * sizeof(float));   
  free_realloc_and_clear(node_f.force_density    , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(lbForceFloat));
#if defined(VIRTUAL_SITES_INERTIALESS_TRACERS) || defined(EK_DEBUG)
  free_realloc_and_clear(node_f.force_density_buf    , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(lbForceFloat));
#endif
#ifdef SHANCHEN
  free_realloc_and_clear(node_f.scforce_density  , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(float));
#endif

  free_realloc_and_clear(nodes_a.seed    , lbpar_gpu->number_of_nodes * sizeof( unsigned int));
  free_realloc_and_clear(nodes_a.boundary, lbpar_gpu->number_of_nodes * sizeof( unsigned int));
  free_realloc_and_clear(nodes_b.seed    , lbpar_gpu->number_of_nodes * sizeof( unsigned int));
  free_realloc_and_clear(nodes_b.boundary, lbpar_gpu->number_of_nodes * sizeof( unsigned int));

  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));

  /**check flag if lb gpu init works*/
  free_realloc_and_clear(gpu_check, sizeof(int));

  if(h_gpu_check!=nullptr)
    free(h_gpu_check);  

  h_gpu_check = (int*)Utils::malloc(sizeof(int));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, (nodes_a, nodes_b));

  #ifdef SHANCHEN
  // TODO FIXME: 
  /* We must add shan-chen forces, which are zero only if the densities are uniform*/
  #endif

  /** calc of velocitydensities from given parameters and initialize the Node_Force array with zero */
  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, (node_f));
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, (nodes_a, device_rho_v, node_f, gpu_check));
 
  intflag = 1;
  current_nodes = &nodes_a;
  h_gpu_check[0] = 0;
  cuda_safe_mem(cudaMemcpy(h_gpu_check, gpu_check, sizeof(int), cudaMemcpyDeviceToHost));
//fprintf(stderr, "initialization of lb gpu code %i\n", lbpar_gpu->number_of_nodes);
  cudaThreadSynchronize();

#if __CUDA_ARCH__ >= 200
  if(!h_gpu_check[0])
  {
    fprintf(stderr, "initialization of lb gpu code failed! \n");
    errexit();
  }
#endif
}

/** reinitialization for the lb gpu fluid called from host
 * @param *lbpar_gpu  Pointer to parameters to setup the lb field
*/
void lb_reinit_GPU(LB_parameters_gpu *lbpar_gpu){

  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
  
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  /** calc of velocity densities from given parameters and initialize the Node_Force array with zero */
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, (nodes_a, device_rho_v, node_f, gpu_check));
}

void lb_realloc_particles_GPU_leftovers(LB_parameters_gpu *lbpar_gpu){

  //copy parameters, especially number of parts to gpu mem
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
}

#ifdef LB_BOUNDARIES_GPU
/** setup and call boundaries from the host
 * @param host_n_lb_boundaries  number of LB boundaries
 * @param number_of_boundnodes  number of boundnodes
 * @param host_boundary_node_list     The indices of the boundary nodes
 * @param host_boundary_index_list    The flag representing the corresponding boundary
 * @param host_lb_boundary_velocity   The constant velocity at the boundary, set by the user (Input)
*/
void lb_init_boundaries_GPU(int host_n_lb_boundaries, int number_of_boundnodes, int *host_boundary_node_list, int* host_boundary_index_list, float* host_lb_boundary_velocity){
  if (this_node != 0) return;
  
  size_of_boundindex = number_of_boundnodes*sizeof(int);
  cuda_safe_mem(cudaMalloc((void**)&boundary_node_list, size_of_boundindex));
  cuda_safe_mem(cudaMalloc((void**)&boundary_index_list, size_of_boundindex));
  cuda_safe_mem(cudaMemcpy(boundary_index_list, host_boundary_index_list, size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(boundary_node_list, host_boundary_node_list, size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMalloc((void**)&lb_boundary_force   , 3*host_n_lb_boundaries*sizeof(float)));
  cuda_safe_mem(cudaMalloc((void**)&lb_boundary_velocity, 3*host_n_lb_boundaries*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(lb_boundary_velocity, host_lb_boundary_velocity, 3*LBBoundaries::lbboundaries.size()*sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, (nodes_a, nodes_b));

  if (LBBoundaries::lbboundaries.size() == 0 && !pdb_boundary_lattice)
  {
    cudaThreadSynchronize();
    return;
  }

  if(number_of_boundnodes == 0)
  {
    fprintf(stderr, "WARNING: boundary cmd executed but no boundary node found!\n");
  }
  else
  {
    int threads_per_block_bound = 64;
    int blocks_per_grid_bound_y = 4;
    int blocks_per_grid_bound_x = (number_of_boundnodes + threads_per_block_bound * blocks_per_grid_bound_y - 1) /(threads_per_block_bound * blocks_per_grid_bound_y);
    dim3 dim_grid_bound = make_uint3(blocks_per_grid_bound_x, blocks_per_grid_bound_y, 1);

    KERNELCALL(init_boundaries, dim_grid_bound, threads_per_block_bound, (boundary_node_list, boundary_index_list, number_of_boundnodes, nodes_a, nodes_b));
  }

  cudaThreadSynchronize();
}
#endif
/**setup and call extern single node force initialization from the host
 * @param *lbpar_gpu    Pointer to host parameter struct
*/
void lb_reinit_extern_nodeforce_GPU(LB_parameters_gpu *lbpar_gpu){

  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu))); 

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, (node_f));

}
/**setup and call extern single node force initialization from the host
 * @param n_extern_node_force_densities       number of nodes on which the external force has to be applied
 * @param *host_extern_node_force_densities   Pointer to the host extern node forces
 * @param *lbpar_gpu                Pointer to host parameter struct
*/
void lb_init_extern_nodeforcedensities_GPU(int n_extern_node_force_densities, LB_extern_nodeforcedensity_gpu *host_extern_node_force_densities, LB_parameters_gpu *lbpar_gpu){

  size_of_extern_node_force_densities = n_extern_node_force_densities*sizeof(LB_extern_nodeforcedensity_gpu);
  cuda_safe_mem(cudaMalloc((void**)&extern_node_force_densities, size_of_extern_node_force_densities));
  cuda_safe_mem(cudaMemcpy(extern_node_force_densities, host_extern_node_force_densities, size_of_extern_node_force_densities, cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu))); 

  int threads_per_block_exf = 64;
  int blocks_per_grid_exf_y = 4;
  int blocks_per_grid_exf_x = (n_extern_node_force_densities + threads_per_block_exf * blocks_per_grid_exf_y - 1) / 
                              (threads_per_block_exf * blocks_per_grid_exf_y);
  dim3 dim_grid_exf = make_uint3(blocks_per_grid_exf_x, blocks_per_grid_exf_y, 1);

  KERNELCALL(init_extern_node_force_densities, dim_grid_exf, threads_per_block_exf, (n_extern_node_force_densities, extern_node_force_densities, node_f));
  cudaFree(extern_node_force_densities);
}

/**setup and call particle kernel from the host
*/
void lb_calc_particle_lattice_ia_gpu(bool couple_virtual){
  if (lbpar_gpu.number_of_particles) 
  {
    /** call of the particle kernel */
    /** values for the particle kernel */
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x = (lbpar_gpu.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1) / 
                                      (threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

    if ( lbpar_gpu.lb_couple_switch & LB_COUPLE_TWO_POINT )
    {
      KERNELCALL( calc_fluid_particle_ia, dim_grid_particles, threads_per_block_particles, 
                  ( *current_nodes, gpu_get_particle_pointer(), 
                    gpu_get_particle_force_pointer(), gpu_get_fluid_composition_pointer(),
                    node_f, gpu_get_particle_seed_pointer(), device_rho_v, couple_virtual )
                );
    }
    else { /** only other option is the three point coupling scheme */
#ifdef SHANCHEN
#if __CUDA_ARCH__ >= 200
      fprintf (stderr, "The three point particle coupling is not currently compatible with the Shan-Chen implementation of the LB\n");
      errexit(); 
#endif
#endif
      KERNELCALL( calc_fluid_particle_ia_three_point_couple, dim_grid_particles, threads_per_block_particles,
                   ( *current_nodes, gpu_get_particle_pointer(),
                     gpu_get_particle_force_pointer(), node_f,
                     gpu_get_particle_seed_pointer(), device_rho_v )
                );
    }
  }
}

/** setup and call kernel for getting macroscopic fluid values of all nodes
 * @param *host_values struct to save the gpu values
*/
void lb_get_values_GPU(LB_rho_v_pi_gpu *host_values){

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) / 
                          (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL( get_mesoscopic_values_in_MD_units, dim_grid, threads_per_block,
              ( *current_nodes, print_rho_v_pi, device_rho_v, node_f ) );
  cuda_safe_mem( cudaMemcpy( host_values, print_rho_v_pi, size_of_rho_v_pi, cudaMemcpyDeviceToHost ) );

}

/** get all the boundary flags for all nodes
 *  @param host_bound_array here go the values of the boundary flag
 */
void lb_get_boundary_flags_GPU(unsigned int* host_bound_array){
   
  unsigned int* device_bound_array;
  cuda_safe_mem(cudaMalloc((void**)&device_bound_array, lbpar_gpu.number_of_nodes*sizeof(unsigned int)));
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) / (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(lb_get_boundaries, dim_grid, threads_per_block, (*current_nodes, device_bound_array));

  cuda_safe_mem(cudaMemcpy(host_bound_array, device_bound_array, lbpar_gpu.number_of_nodes*sizeof(unsigned int), cudaMemcpyDeviceToHost));

  cudaFree(device_bound_array);

}

/** setup and call kernel for getting macroscopic fluid values of a single node*/
void lb_print_node_GPU(int single_nodeindex, LB_rho_v_pi_gpu *host_print_values){ 
      
  LB_rho_v_pi_gpu *device_print_values;
  cuda_safe_mem(cudaMalloc((void**)&device_print_values, sizeof(LB_rho_v_pi_gpu)));
  int threads_per_block_print = 1;
  int blocks_per_grid_print_y = 1;
  int blocks_per_grid_print_x = 1;
  dim3 dim_grid_print = make_uint3(blocks_per_grid_print_x, blocks_per_grid_print_y, 1);

  KERNELCALL(lb_print_node, dim_grid_print, threads_per_block_print, (single_nodeindex, device_print_values, *current_nodes, device_rho_v, node_f));

  cuda_safe_mem(cudaMemcpy(host_print_values, device_print_values, sizeof(LB_rho_v_pi_gpu), cudaMemcpyDeviceToHost));
  cudaFree(device_print_values);

}

/** setup and call kernel to calculate the total momentum of the hole fluid
 * @param *mass value of the mass calcutated on the GPU
*/
void lb_calc_fluid_mass_GPU(double* mass){

  float* tot_mass;
  float cpu_mass =  0.0f ;
  cuda_safe_mem(cudaMalloc((void**)&tot_mass, sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_mass, &cpu_mass, sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(calc_mass, dim_grid, threads_per_block,(*current_nodes, tot_mass));

  cuda_safe_mem(cudaMemcpy(&cpu_mass, tot_mass, sizeof(float), cudaMemcpyDeviceToHost));
  
  cudaFree(tot_mass);
  mass[0] = (double)(cpu_mass);
}

/** setup and call kernel to calculate the total momentum of the whole fluid
 *  @param host_mom value of the momentum calcutated on the GPU
 */
void lb_calc_fluid_momentum_GPU(double* host_mom){

  float* tot_momentum;
  float host_momentum[3] = { 0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void**)&tot_momentum, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3*sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3*sizeof(float), cudaMemcpyDeviceToHost));
  
  cudaFree(tot_momentum);
  host_mom[0] = (double)(host_momentum[0]* lbpar_gpu.agrid/lbpar_gpu.tau);
  host_mom[1] = (double)(host_momentum[1]* lbpar_gpu.agrid/lbpar_gpu.tau);
  host_mom[2] = (double)(host_momentum[2]* lbpar_gpu.agrid/lbpar_gpu.tau);
}

/** setup and call kernel to remove the net momentum of the whole fluid
 */
void lb_remove_fluid_momentum_GPU(void){
  float* tot_momentum;
  float host_momentum[3] = { 0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void**)&tot_momentum, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3*sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3*sizeof(float), cudaMemcpyDeviceToHost));

  KERNELCALL(remove_momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cudaFree(tot_momentum);
}


/** setup and call kernel to calculate the temperature of the hole fluid
 *  @param host_temp value of the temperatur calcutated on the GPU
*/
void lb_calc_fluid_temperature_GPU(double* host_temp){

  int host_number_of_non_boundary_nodes = 0;
  int *device_number_of_non_boundary_nodes;
  cuda_safe_mem(cudaMalloc((void**)&device_number_of_non_boundary_nodes, sizeof(int)));
  cuda_safe_mem(cudaMemcpy(device_number_of_non_boundary_nodes, &host_number_of_non_boundary_nodes, sizeof(int), cudaMemcpyHostToDevice));

  float host_jsquared = 0.0f;
  float* device_jsquared;
  cuda_safe_mem(cudaMalloc((void**)&device_jsquared, sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_jsquared, &host_jsquared, sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(temperature, dim_grid, threads_per_block,(*current_nodes, device_jsquared, device_number_of_non_boundary_nodes));

  cuda_safe_mem(cudaMemcpy(&host_number_of_non_boundary_nodes, device_number_of_non_boundary_nodes, sizeof(int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(&host_jsquared, device_jsquared, sizeof(float), cudaMemcpyDeviceToHost));

  // TODO: check that temperature calculation is properly implemented for shanchen
  *host_temp=0;

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      *host_temp += (double)(host_jsquared*1./(3.0f*lbpar_gpu.rho[ii]*host_number_of_non_boundary_nodes*lbpar_gpu.tau*lbpar_gpu.tau*lbpar_gpu.agrid));
  }
}


#ifdef SHANCHEN
void lb_calc_shanchen_GPU(){
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

#ifdef LB_BOUNDARIES_GPU
  if (LBBoundaries::lbboundaries.size() != 0)
  {
    KERNELCALL(lb_shanchen_set_boundaries, dim_grid, threads_per_block,(*current_nodes));
    cudaThreadSynchronize();
  }
#endif
  KERNELCALL(lb_shanchen_GPU, dim_grid, threads_per_block,(*current_nodes, node_f));
}

#endif // SHANCHEN

/** setup and call kernel for getting macroscopic fluid values of all nodes
 * @param *host_checkpoint_vd struct to save the gpu populations
 * @param *host_checkpoint_seed struct to save the nodes' seeds for the lb on the gpu
 * @param *host_checkpoint_boundary struct to save the boundary nodes
 * @param *host_checkpoint_force struct to save the forces on the nodes
 */
void lb_save_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, unsigned int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force){

  cuda_safe_mem(cudaMemcpy(host_checkpoint_vd, current_nodes->vd, lbpar_gpu.number_of_nodes * 19 * sizeof(float), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_seed, current_nodes->seed, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_boundary, current_nodes->boundary, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_force, node_f.force_density, lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat), cudaMemcpyDeviceToHost));

}

/** setup and call kernel for setting macroscopic fluid values of all nodes
 * @param *host_checkpoint_vd struct to save the gpu populations
 * @param *host_checkpoint_seed struct to save the nodes' seeds for the lb on the gpu
 * @param *host_checkpoint_boundary struct to save the boundary nodes
 * @param *host_checkpoint_force struct to save the forces on the nodes
*/
void lb_load_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, unsigned int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force){

  current_nodes = &nodes_a;
  intflag = 1;

  cuda_safe_mem(cudaMemcpy(current_nodes->vd, host_checkpoint_vd, lbpar_gpu.number_of_nodes * 19 * sizeof(float), cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaMemcpy(current_nodes->seed, host_checkpoint_seed, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(current_nodes->boundary, host_checkpoint_boundary, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(node_f.force_density, host_checkpoint_force, lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat), cudaMemcpyHostToDevice));
}

/** setup and call kernel to get the boundary flag of a single node
 *  @param single_nodeindex number of the node to get the flag for
 *  @param host_flag her goes the value of the boundary flag
 */
void lb_get_boundary_flag_GPU(int single_nodeindex, unsigned int* host_flag){
   
  unsigned int* device_flag;
  cuda_safe_mem(cudaMalloc((void**)&device_flag, sizeof(unsigned int)));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(lb_get_boundary_flag, dim_grid_flag, threads_per_block_flag, (single_nodeindex, device_flag, *current_nodes));

  cuda_safe_mem(cudaMemcpy(host_flag, device_flag, sizeof(unsigned int), cudaMemcpyDeviceToHost));

  cudaFree(device_flag);
}

/** set the density at a single node
 *  @param single_nodeindex the node to set the velocity for 
 *  @param *host_rho the density to set
 */
void lb_set_node_rho_GPU(int single_nodeindex, float* host_rho){
   
  float* device_rho;
  cuda_safe_mem(cudaMalloc((void**)&device_rho, LB_COMPONENTS*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_rho, host_rho, LB_COMPONENTS*sizeof(float), cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);
  KERNELCALL(set_rho, dim_grid_flag, threads_per_block_flag, (*current_nodes, device_rho_v, single_nodeindex, device_rho)); 
  cudaFree(device_rho);
}

/** set the net velocity at a single node
 *  @param single_nodeindex the node to set the velocity for 
 *  @param host_velocity the velocity to set
 */
void lb_set_node_velocity_GPU(int single_nodeindex, float* host_velocity){
   
  float* device_velocity;
  cuda_safe_mem(cudaMalloc((void**)&device_velocity, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_velocity, host_velocity, 3*sizeof(float), cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(set_u_from_rho_v_pi, dim_grid_flag, threads_per_block_flag, (*current_nodes, single_nodeindex, device_velocity, device_rho_v, node_f));

  cudaFree(device_velocity);
}

/** reinit of params 
 * @param *lbpar_gpu struct containing the paramters of the fluid
*/
void reinit_parameters_GPU(LB_parameters_gpu *lbpar_gpu){
  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
}

/**integration kernel for the lb gpu fluid update called from host */
void lb_integrate_GPU() {

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

#ifdef LB_BOUNDARIES_GPU
  if (LBBoundaries::lbboundaries.size() > 0)
  {
    cuda_safe_mem(cudaMemset( lb_boundary_force, 0, 3*LBBoundaries::lbboundaries.size()*sizeof(float)));
  }
#endif

  /**call of fluid step*/
  /* NOTE: if pi is needed at every integration step, one should call an extended version 
           of the integrate kernel, or pass also device_rho_v_pi and make sure that either 
           it or device_rho_v are nullptr depending on extended_values_flag */ 
  if (intflag == 1)
  {
    KERNELCALL(integrate, dim_grid, threads_per_block, (nodes_a, nodes_b, device_rho_v, node_f, lb_ek_parameters_gpu));
    current_nodes = &nodes_b;
    intflag = 0;
  }
  else
  {
    KERNELCALL(integrate, dim_grid, threads_per_block, (nodes_b, nodes_a, device_rho_v, node_f, lb_ek_parameters_gpu));
    current_nodes = &nodes_a;
    intflag = 1;
  }

#ifdef LB_BOUNDARIES_GPU
  if (LBBoundaries::lbboundaries.size() > 0)
    {
      KERNELCALL(apply_boundaries, dim_grid, threads_per_block, (*current_nodes, lb_boundary_velocity, lb_boundary_force));
    }
#endif
}

void lb_gpu_get_boundary_forces(double* forces) {
#ifdef LB_BOUNDARIES_GPU
  float* temp = (float*) Utils::malloc(3*LBBoundaries::lbboundaries.size()*sizeof(float));
  cuda_safe_mem(cudaMemcpy(temp, lb_boundary_force, 3*LBBoundaries::lbboundaries.size()*sizeof(float), cudaMemcpyDeviceToHost));

  for (int i =0; i<3*LBBoundaries::lbboundaries.size(); i++)
  {
    forces[i]=(double)temp[i];
  }
  free(temp);
#endif
}

struct lb_lbfluid_mass_of_particle
{
  __device__ float operator()(CUDA_particle_data particle) const
  {
#ifdef MASS
    return particle.mass;
#else
    return 1.;
#endif
  };
};

void lb_lbfluid_remove_total_momentum()
{
  // calculate momentum of fluid and particles
  float total_momentum[3] = { 0.0f, 0.0f, 0.0f };
  lb_lbfluid_calc_linear_momentum(total_momentum, /*include_particles*/ 1, /*include_lbfluid*/ 1);

  thrust::device_ptr<CUDA_particle_data> ptr(gpu_get_particle_pointer());
  float particles_mass = thrust::transform_reduce(
    ptr,
    ptr + lbpar_gpu.number_of_particles,
    lb_lbfluid_mass_of_particle(),
    0.0f,
    thrust::plus<float>());

  // lb_calc_fluid_mass_GPU has to be called with double but we don't
  // want narrowing warnings, that's why we narrow it down by hand.
  double lb_calc_fluid_mass_res;
  lb_calc_fluid_mass_GPU( &lb_calc_fluid_mass_res );
  float fluid_mass = lb_calc_fluid_mass_res;

  /* Momentum fraction of the particles */
  auto const part_frac = particles_mass / (fluid_mass + particles_mass);
  /* Mometum per particle */
  float momentum_particles[3] = {
    -total_momentum[0]*part_frac,
    -total_momentum[1]*part_frac,
    -total_momentum[2]*part_frac
  };

  auto const fluid_frac = fluid_mass / (fluid_mass + particles_mass);
  float momentum_fluid[3] = {
    -total_momentum[0]*fluid_frac,
    -total_momentum[1]*fluid_frac,
    -total_momentum[2]*fluid_frac
  };

  lb_lbfluid_particles_add_momentum( momentum_particles );
  lb_lbfluid_fluid_add_momentum( momentum_fluid );
}

__global__ void lb_lbfluid_fluid_add_momentum_kernel(
  float momentum[3],
  LB_nodes_gpu n_a,
  LB_node_force_density_gpu node_f,
  LB_rho_v_gpu *d_v)
{
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int number_of_nodes = para.number_of_nodes;
#ifdef LB_BOUNDARIES_GPU
  number_of_nodes -= para.number_of_boundnodes;
#endif
  if( index < para.number_of_nodes )
  {
    if( n_a.boundary[index] == 0 )
    {
      float force_factor=powf(para.agrid,2)*para.tau*para.tau;
      for(int i = 0 ; i < LB_COMPONENTS ; ++i )
      {
        // add force density onto each node (momentum / time_step / Volume)
        node_f.force_density[(0+i*3)*para.number_of_nodes + index] += momentum[0] / para.tau / (number_of_nodes * powf(para.agrid,3)) * force_factor;
        node_f.force_density[(1+i*3)*para.number_of_nodes + index] += momentum[1] / para.tau / (number_of_nodes * powf(para.agrid,3)) * force_factor;
        node_f.force_density[(2+i*3)*para.number_of_nodes + index] += momentum[2] / para.tau / (number_of_nodes * powf(para.agrid,3)) * force_factor;
      }
    }
  }
}


void lb_lbfluid_fluid_add_momentum( float momentum_host[3] )
{
  float* momentum_device;
  cuda_safe_mem(cudaMalloc((void**)&momentum_device,3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(momentum_device, momentum_host, 3*sizeof(float), cudaMemcpyHostToDevice));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1)/(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL( lb_lbfluid_fluid_add_momentum_kernel, dim_grid, threads_per_block, (momentum_device, *current_nodes, node_f, device_rho_v));
}


/**set the populations of a specific node on the GPU
 * @param n_a            Pointer to local node residing in array a (Input)
 * @param population     Pointer to new population (Input)
 * @param x              x-coordinate of node (Input)
 * @param y              y-coordinate of node (Input)
 * @param z              z-coordinate of node (Input)
 * @param c              LB component (for SHANCHEN) (Input)
*/
__global__ void lb_lbfluid_set_population_kernel(LB_nodes_gpu n_a, float population[LBQ], int x, int y, int z, int c)
{
  int xyz[3] = { x, y, z };
  int index = xyz_to_index( xyz );

  for (int i = 0; i < LBQ; ++i)
  {
    n_a.vd[( i + c*LBQ ) * para.number_of_nodes + index] = population[i];
  }
}


/**interface to set the populations of a specific node for the GPU
 * @param xyz            coordinates of node (Input)
 * @param population     Pointer to population (Input)
 * @param c              LB component (for SHANCHEN) (Input)
*/
void lb_lbfluid_set_population( int xyz[3], float population_host[LBQ], int c )
{
  float* population_device;
  cuda_safe_mem(cudaMalloc((void**)&population_device,LBQ*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(population_device, population_host, LBQ*sizeof(float), cudaMemcpyHostToDevice));

  dim3 dim_grid = make_uint3(1, 1, 1);
  KERNELCALL( lb_lbfluid_set_population_kernel, dim_grid, 1,
              (*current_nodes, population_device, xyz[0], xyz[1], xyz[2], c));

  cuda_safe_mem(cudaFree(population_device));
}


/**get the populations of a specific node on the GPU
 * @param n_a            Pointer to local node residing in array a (Input)
 * @param population     Pointer to population (Output)
 * @param x              x-coordinate of node (Input)
 * @param y              y-coordinate of node (Input)
 * @param z              z-coordinate of node (Input)
 * @param c              LB component (for SHANCHEN) (Input)
*/
__global__ void lb_lbfluid_get_population_kernel(LB_nodes_gpu n_a, float population[LBQ], int x, int y, int z, int c)
{
  int xyz[3] = { x, y, z };
  int index = xyz_to_index( xyz );

  for (int i = 0; i < LBQ; ++i)
  {
    population[i] = n_a.vd[( i + c*LBQ ) * para.number_of_nodes + index];
  }
}


/**interface to get the populations of a specific node for the GPU
 * @param xyz            coordinates of node (Input)
 * @param population     Pointer to population (Output)
 * @param c              LB component (for SHANCHEN) (Input)
*/
void lb_lbfluid_get_population( int xyz[3], float population_host[LBQ], int c )
{
  float* population_device;
  cuda_safe_mem(cudaMalloc((void**)&population_device,LBQ*sizeof(float)));

  dim3 dim_grid = make_uint3(1, 1, 1);
  KERNELCALL( lb_lbfluid_get_population_kernel, dim_grid, 1,
              (*current_nodes, population_device, xyz[0], xyz[1], xyz[2], c));

  cuda_safe_mem(cudaMemcpy(population_host, population_device, LBQ*sizeof(float), cudaMemcpyDeviceToHost));

  cuda_safe_mem(cudaFree(population_device));
}

struct two_point_interpolation { 
    LB_nodes_gpu current_nodes_gpu;
    LB_rho_v_gpu *d_v_gpu;
    two_point_interpolation(LB_nodes_gpu _current_nodes_gpu, LB_rho_v_gpu *_d_v_gpu) : current_nodes_gpu(_current_nodes_gpu), d_v_gpu(_d_v_gpu) {};
	__device__ float3 operator()(const float3 &position) const {
        unsigned int node_index[8];
        float delta[8];
        float u[3];
        float mode[19*LB_COMPONENTS];
        float _position[3] = {position.x, position.y, position.z};
        interpolation_two_point_coupling(current_nodes_gpu, _position, node_index, mode, d_v_gpu, delta, u);
        return make_float3(u[0], u[1], u[2]);
	} 
};

void lb_lbfluid_get_interpolated_velocity_at_positions(double const *positions, double *velocities, int length) {
    thrust::host_vector<float3> positions_host(length);
    for (int p=0; p < 3 * length; p+=3) {
        // Cast double coming from python to float.
        positions_host[p/3].x = static_cast<float>(positions[p]);
        positions_host[p/3].y = static_cast<float>(positions[p+1]);
        positions_host[p/3].z = static_cast<float>(positions[p+2]);
    }
    thrust::device_vector<float3> positions_device = positions_host;
    thrust::device_vector<float3> velocities_device(length);
    thrust::transform(positions_device.begin(), positions_device.end(), velocities_device.begin(), two_point_interpolation(*current_nodes, device_rho_v));
    thrust::host_vector<float3> velocities_host = velocities_device;
    int index = 0;
    for (auto v : velocities_host) {
        velocities[index] = static_cast<double>(v.x) * lbpar_gpu.agrid/lbpar_gpu.tau;
        velocities[index+1] = static_cast<double>(v.y) * lbpar_gpu.agrid/lbpar_gpu.tau;
        velocities[index+2] = static_cast<double>(v.z) * lbpar_gpu.agrid/lbpar_gpu.tau;
        index += 3;
    }
}

#endif /* LB_GPU */
