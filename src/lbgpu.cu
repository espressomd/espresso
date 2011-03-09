/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
 */

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>

extern "C" {
#include "lbgpu.h"
}

#ifdef LB_GPU

/**defining structures residing in global memory */
/** struct for phys. values */
static LB_values_gpu *device_values = NULL;
/** structs for velocity densities */
static LB_nodes_gpu nodes_a;
static LB_nodes_gpu nodes_b;
/** struct for particle force */
static LB_particle_force_gpu *particle_force = NULL;
/** struct for particle position and veloctiy */
static LB_particle_gpu *particle_data = NULL;
/** struct for node force */
static LB_node_force_gpu node_f;
/** struct for storing particle rn seed */
static LB_particle_seed_gpu *part = NULL;

static LB_extern_nodeforce_gpu *extern_nodeforces = NULL;
/** pointer for bound index array*/
static int *boundindex;
/** pointers for additional cuda check flag*/
static int *gpu_check = NULL;
static int *h_gpu_check = NULL;

/** values for the kernel call */
static int threads_per_block;
static int blocks_per_grid;

/** values for the particle kernel */
static int threads_per_block_particles;
static int blocks_per_grid_particles;

/** values for the boundary init kernel */
static int threads_per_block_bound;
static int blocks_per_grid_bound;

static int threads_per_block_exf;
static int blocks_per_grid_exf;

static int threads_per_block_print;
static int blocks_per_grid_print;
static unsigned int intflag = 1;
/**defining size values for allocating global memory */
static size_t size_of_values;
static size_t size_of_forces;
static size_t size_of_positions;
static size_t size_of_seed;
static size_t size_of_boundindex;
static size_t size_of_extern_nodeforces;

/**parameters residing in constant memory */
static __constant__ LB_parameters_gpu para;
static __constant__ int number_of_bnodes;
static __constant__ float c_sound_sq = 1.f/3.f;

/**cudasteams for parallel computing on cpu and gpu */
cudaStream_t stream[1];

cudaError_t err;
/*-------------------------------------------------------*/
/*********************************************************/
/**device funktions call by kernel funktions */
/*********************************************************/
/*-------------------------------------------------------*/


/*-------------------------------------------------------*/

/** atomic add function for sveral cuda architectures */
/*@{
 * @param 
}*/
/*-------------------------------------------------------*/
__device__ inline void atomicadd(float* address, float value){
#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110
#warning Using slower atomicAdd emulation
// float-atomic-add from 
// [url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]http://forums.nvidia.com/index.php?showtop...st&p=991561[/url]
  float old = value;
  while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
#else
#error I need at least compute capability 1.1
#endif
}
/*-------------------------------------------------------*/
/**randomgenerator which generates numbers [0,1] */
/*@{
 * @param *rn	Pointer to randomnumber array of the local node or particle 
}*/
/*-------------------------------------------------------*/
__device__ void random_01(LB_randomnr_gpu *rn){

  const float mxi = 1.f/(float)(1ul<<31);
  unsigned int curr = rn->seed;

  curr = 1103515245 * curr + 12345;
  rn->randomnr[0] = (float)(curr & ((1ul<<31)-1))*mxi;
  curr = 1103515245 * curr + 12345;
  rn->randomnr[1] = (float)(curr & ((1ul<<31)-1))*mxi;
  rn->seed = curr;

}
/*-------------------------------------------------------*/
/** gaussian random nummber generator for thermalisation */
/*@{
 * @param *rn	Pointer to randomnumber array of the local node node or particle 
}*/
/*-------------------------------------------------------*/
__device__ void gaussian_random(LB_randomnr_gpu *rn){

  float x1, x2;
  float r2, fac;
  /* On every second call two gaussian random numbers are calculated
   via the Box-Muller transformation.*/
  /* draw two uniform random numbers in the unit circle */
  do {
    random_01(rn);
    x1 = 2.f*rn->randomnr[0]-1.f;
    x2 = 2.f*rn->randomnr[1]-1.f;
    r2 = x1*x1 + x2*x2;
  } while (r2 >= 1.f || r2 == 0.f);

  /* perform Box-Muller transformation */
  fac = sqrtf(-2.f*__logf(r2)/r2);
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  
}
/*-------------------------------------------------------*/
/**tranformation from 1d array-index to xyz */
/*@{
 * @param index		node index / thread index (Input)
 * @param xyz		Pointer to calculated xyz array (Output)
 */
/*-------------------------------------------------------*/
__device__ void index_to_xyz(unsigned int index, unsigned int *xyz){

  xyz[0] = index%para.dim_x;
  index /= para.dim_x;
  xyz[1] = index%para.dim_y;
  index /= para.dim_y;
  xyz[2] = index;
}
/*-------------------------------------------------------*/
/**calculation of the modes from the velocitydensities (space-transform.)*/
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Output)
}*/
/*-------------------------------------------------------*/
__device__ void calc_m_from_n(LB_nodes_gpu n_a, unsigned int index, float *mode){

  /* mass mode */
  mode[0] = n_a.vd[0][index] + n_a.vd[1][index] + n_a.vd[2][index]
          + n_a.vd[3][index] + n_a.vd[4][index] + n_a.vd[5][index]
          + n_a.vd[6][index] + n_a.vd[7][index] + n_a.vd[8][index]
          + n_a.vd[9][index] + n_a.vd[10][index] + n_a.vd[11][index] + n_a.vd[12][index]
          + n_a.vd[13][index] + n_a.vd[14][index] + n_a.vd[15][index] + n_a.vd[16][index]
          + n_a.vd[17][index] + n_a.vd[18][index];

  /* momentum modes */
  mode[1] = (n_a.vd[1][index] - n_a.vd[2][index]) + (n_a.vd[7][index] - n_a.vd[8][index])
          + (n_a.vd[9][index] - n_a.vd[10][index]) + (n_a.vd[11][index] - n_a.vd[12][index])
          + (n_a.vd[13][index] - n_a.vd[14][index]);
  mode[2] = (n_a.vd[3][index] - n_a.vd[4][index]) + (n_a.vd[7][index] - n_a.vd[8][index])
          - (n_a.vd[9][index] - n_a.vd[10][index]) + (n_a.vd[15][index] - n_a.vd[16][index])
          + (n_a.vd[17][index] - n_a.vd[18][index]);
  mode[3] = (n_a.vd[5][index] - n_a.vd[6][index]) + (n_a.vd[11][index] - n_a.vd[12][index])
          - (n_a.vd[13][index] - n_a.vd[14][index]) + (n_a.vd[15][index] - n_a.vd[16][index])
          - (n_a.vd[17][index] - n_a.vd[18][index]);

  /* stress modes */
  mode[4] = -(n_a.vd[0][index]) + n_a.vd[7][index] + n_a.vd[8][index] + n_a.vd[9][index] + n_a.vd[10][index]
          + n_a.vd[11][index] + n_a.vd[12][index] + n_a.vd[13][index] + n_a.vd[14][index]
          + n_a.vd[15][index] + n_a.vd[16][index] + n_a.vd[17][index] + n_a.vd[18][index];
  mode[5] = n_a.vd[1][index] + n_a.vd[2][index] - (n_a.vd[3][index] + n_a.vd[4][index])
          + (n_a.vd[11][index] + n_a.vd[12][index]) + (n_a.vd[13][index] + n_a.vd[14][index])
          - (n_a.vd[15][index] + n_a.vd[16][index]) - (n_a.vd[17][index] + n_a.vd[18][index]);
  mode[6] = (n_a.vd[1][index] + n_a.vd[2][index]) + (n_a.vd[3][index] + n_a.vd[4][index])
          - (n_a.vd[11][index] + n_a.vd[12][index]) - (n_a.vd[13][index] + n_a.vd[14][index])
          - (n_a.vd[15][index] + n_a.vd[16][index]) - (n_a.vd[17][index] + n_a.vd[18][index])
          - 2.f*(n_a.vd[5][index] + n_a.vd[6][index] - (n_a.vd[7][index] + n_a.vd[8][index])
          - (n_a.vd[9][index] +n_a.vd[10][index]));
  mode[7] = n_a.vd[7][index] + n_a.vd[8][index] - (n_a.vd[9][index] + n_a.vd[10][index]);
  mode[8] = n_a.vd[11][index] + n_a.vd[12][index] - (n_a.vd[13][index] + n_a.vd[14][index]);
  mode[9] = n_a.vd[15][index] + n_a.vd[16][index] - (n_a.vd[17][index] + n_a.vd[18][index]);

  /* kinetic modes */
  mode[10] = -2.f*(n_a.vd[1][index] - n_a.vd[2][index]) + (n_a.vd[7][index] - n_a.vd[8][index])
           + (n_a.vd[9][index] - n_a.vd[10][index]) + (n_a.vd[11][index] - n_a.vd[12][index])
           + (n_a.vd[13][index] - n_a.vd[14][index]);
  mode[11] = -2.f*(n_a.vd[3][index] - n_a.vd[4][index]) + (n_a.vd[7][index] - n_a.vd[8][index])
           - (n_a.vd[9][index] - n_a.vd[10][index]) + (n_a.vd[15][index] - n_a.vd[16][index])
           + (n_a.vd[17][index] - n_a.vd[18][index]);
  mode[12] = -2.f*(n_a.vd[5][index] - n_a.vd[6][index]) + (n_a.vd[11][index] - n_a.vd[12][index])
           - (n_a.vd[13][index] - n_a.vd[14][index]) + (n_a.vd[15][index] - n_a.vd[16][index])
           - (n_a.vd[17][index] - n_a.vd[18][index]);
  mode[13] = (n_a.vd[7][index] - n_a.vd[8][index]) + (n_a.vd[9][index] - n_a.vd[10][index])
           - (n_a.vd[11][index] - n_a.vd[12][index]) - (n_a.vd[13][index] - n_a.vd[14][index]);
  mode[14] = (n_a.vd[7][index] - n_a.vd[8][index]) - (n_a.vd[9][index] - n_a.vd[10][index])
           - (n_a.vd[15][index] - n_a.vd[16][index]) - (n_a.vd[17][index] - n_a.vd[18][index]);
  mode[15] = (n_a.vd[11][index] - n_a.vd[12][index]) - (n_a.vd[13][index] - n_a.vd[14][index])
           - (n_a.vd[15][index] - n_a.vd[16][index]) + (n_a.vd[17][index] - n_a.vd[18][index]);
  mode[16] = n_a.vd[0][index] + n_a.vd[7][index] + n_a.vd[8][index] + n_a.vd[9][index] + n_a.vd[10][index]
           + n_a.vd[11][index] + n_a.vd[12][index] + n_a.vd[13][index] + n_a.vd[14][index]
           + n_a.vd[15][index] + n_a.vd[16][index] + n_a.vd[17][index] + n_a.vd[18][index]
           - 2.f*((n_a.vd[1][index] + n_a.vd[2][index]) + (n_a.vd[3][index] + n_a.vd[4][index])
           + (n_a.vd[5][index] + n_a.vd[6][index]));
  mode[17] = -(n_a.vd[1][index] + n_a.vd[2][index]) + (n_a.vd[3][index] + n_a.vd[4][index])
           + (n_a.vd[11][index] + n_a.vd[12][index]) + (n_a.vd[13][index] + n_a.vd[14][index])
           - (n_a.vd[15][index] + n_a.vd[16][index]) - (n_a.vd[17][index] + n_a.vd[18][index]);
  mode[18] = -(n_a.vd[1][index] + n_a.vd[2][index]) - (n_a.vd[3][index] + n_a.vd[4][index])
           - (n_a.vd[11][index] + n_a.vd[12][index]) - (n_a.vd[13][index] + n_a.vd[14][index])
           - (n_a.vd[15][index] + n_a.vd[16][index]) - (n_a.vd[17][index] + n_a.vd[18][index])
           + 2.f*((n_a.vd[5][index] + n_a.vd[6][index]) + (n_a.vd[7][index] + n_a.vd[8][index])
           + (n_a.vd[9][index] + n_a.vd[10][index]));

}
/*-------------------------------------------------------*/
/**lb_relax_modes, means collision*/
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Input/Output)
 * @param node_f	Pointer to local node force (Input)
}*/
/*-------------------------------------------------------*/
__device__ void relax_modes(float *mode, unsigned int index, LB_node_force_gpu node_f){

  float Rho = mode[0] + para.rho*para.agrid*para.agrid*para.agrid;
  float j[3], pi_eq[6];

  /* re-construct the real density
  * remember that the populations are stored as differences to their
  * equilibrium value */

  j[0] = mode[1];
  j[1] = mode[2];
  j[2] = mode[3];

  /* if forces are present, the momentum density is redefined to
  * inlcude one half-step of the force action.  See the
  * Chapman-Enskog expansion in [Ladd & Verberg]. */

  j[0] += 0.5f*node_f.force[0][index];
  j[1] += 0.5f*node_f.force[1][index];
  j[2] += 0.5f*node_f.force[2][index];

  /* equilibrium part of the stress modes (eq13 schiller)*/
  pi_eq[0] = ((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2]))/Rho;
  pi_eq[1] = ((j[0]*j[0])-(j[1]*j[1]))/Rho;
  pi_eq[2] = (((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2])) - 3.0f*(j[2]*j[2]))/Rho;
  pi_eq[3] = j[0]*j[1]/Rho;
  pi_eq[4] = j[0]*j[2]/Rho;
  pi_eq[5] = j[1]*j[2]/Rho;

  /* relax the stress modes (eq14 schiller)*/
  mode[4] = pi_eq[0] + para.gamma_bulk*(mode[4] - pi_eq[0]);
  mode[5] = pi_eq[1] + para.gamma_shear*(mode[5] - pi_eq[1]);
  mode[6] = pi_eq[2] + para.gamma_shear*(mode[6] - pi_eq[2]);
  mode[7] = pi_eq[3] + para.gamma_shear*(mode[7] - pi_eq[3]);
  mode[8] = pi_eq[4] + para.gamma_shear*(mode[8] - pi_eq[4]);
  mode[9] = pi_eq[5] + para.gamma_shear*(mode[9] - pi_eq[5]);

  /* relax the ghost modes (project them out) */
  /* ghost modes have no equilibrium part due to orthogonality */
  mode[10] = para.gamma_odd*mode[10];
  mode[11] = para.gamma_odd*mode[11];
  mode[12] = para.gamma_odd*mode[12];
  mode[13] = para.gamma_odd*mode[13];
  mode[14] = para.gamma_odd*mode[14];
  mode[15] = para.gamma_odd*mode[15];
  mode[16] = para.gamma_even*mode[16];
  mode[17] = para.gamma_even*mode[17];
  mode[18] = para.gamma_even*mode[18];

}
/*-------------------------------------------------------*/
/**thermalization of the modes with gaussian random numbers*/
/*@{
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Input/Output)
 * @param *rn		Pointer to randomnumber array of the local node
 */
/*-------------------------------------------------------*/
__device__ void thermalize_modes(float *mode, unsigned int index, LB_randomnr_gpu *rn){

  float rootrho = sqrt(mode[0]+para.rho*para.agrid*para.agrid*para.agrid);

  /* stress modes */
  gaussian_random(rn);
  mode[4] += rootrho*(para.mu*(2.f/3.f)*(1.f-(para.gamma_bulk*para.gamma_bulk))) * rn->randomnr[1];
  mode[5] += rootrho*(para.mu*(4.f/9.f)*(1.f-(para.gamma_shear*para.gamma_shear))) * rn->randomnr[0];

  gaussian_random(rn);
  mode[6] += rootrho*(para.mu*(4.f/3.f)*(1.f-(para.gamma_shear*para.gamma_shear))) * rn->randomnr[1];
  mode[7] += rootrho*(para.mu*(1.f/9.f)*(1.f-(para.gamma_shear*para.gamma_shear))) * rn->randomnr[0];

  gaussian_random(rn);
  mode[8] += rootrho*(para.mu*(1.f/9.f)*(1.f-(para.gamma_shear*para.gamma_shear))) * rn->randomnr[1];
  mode[9] += rootrho*(para.mu*(1.f/9.f)*(1.f-(para.gamma_shear*para.gamma_shear))) * rn->randomnr[0];
 
  /* ghost modes */
  gaussian_random(rn);
  mode[10] += rootrho*(para.mu*(2.f/3.f)) * rn->randomnr[1];
  mode[11] += rootrho*(para.mu*(2.f/3.f)) * rn->randomnr[0];

  gaussian_random(rn);
  mode[12] += rootrho*(para.mu*(2.f/3.f)) * rn->randomnr[1];
  mode[13] += rootrho*(para.mu*(2.f/9.f)) * rn->randomnr[0];

  gaussian_random(rn);
  mode[14] += rootrho*(para.mu*(2.f/9.f)) * rn->randomnr[1];
  mode[15] += rootrho*(para.mu*(2.f/9.f)) * rn->randomnr[0];

  gaussian_random(rn);
  mode[16] += rootrho*(para.mu*(2.f)) * rn->randomnr[1];
  mode[17] += rootrho*(para.mu*(4.f/9.f)) * rn->randomnr[0];

  gaussian_random(rn);
  mode[18] += rootrho*(para.mu*(4.f/3.f)) * rn->randomnr[1];

}
/*-------------------------------------------------------*/
/**normalization of the modes need befor backtransformation into velocity space*/
/*@{
 * @param mode		Pointer to the local register values mode (Input/Output)
}*/
/*-------------------------------------------------------*/
__device__ void normalize_modes(float* mode){

  /* normalization factors enter in the back transformation */
  /* the following values are the (weighted) lengths of the vectors */
  mode[0] *= 1.f;
  mode[1] *= 3.f;
  mode[2] *= 3.f;
  mode[3] *= 3.f;
  mode[4] *= 3.f/2.f;
  mode[5] *= 9.f/4.f;
  mode[6] *= 3.f/4.f;
  mode[7] *= 9.f;
  mode[8] *= 9.f;
  mode[9] *= 9.f;
  mode[10] *= 3.f/2.f;
  mode[11] *= 3.f/2.f;
  mode[12] *= 3.f/2.f;
  mode[13] *= 9.f/2.f;
  mode[14] *= 9.f/2.f;
  mode[15] *= 9.f/2.f;
  mode[16] *= 1.f/2.f;
  mode[17] *= 9.f/4.f;
  mode[18] *= 3.f/4.f;

}
/*-------------------------------------------------------*/
/**backtransformation from modespace to desityspace and streaming with the push method using pbc*/
/*@{
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Input)
 * @param *n_b		Pointer to local node residing in array b (Output)
}*/
/*-------------------------------------------------------*/
__device__ void calc_n_from_modes_push(LB_nodes_gpu n_b, float *mode, unsigned int index){

  unsigned int xyz[3];
  index_to_xyz(index, xyz);
  unsigned int x = xyz[0];
  unsigned int y = xyz[1];
  unsigned int z = xyz[2];

  n_b.vd[0][x + para.dim_x*y + para.dim_x*para.dim_y*z] = 1.f/3.f * (mode[0] - mode[4] + mode[16]);
  n_b.vd[1][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = 1.f/18.f * (mode[0] + mode[1] + mode[5] + mode[6] - mode[17] - mode[18] - 2.f*(mode[10] + mode[16]));
  n_b.vd[2][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = 1.f/18.f * (mode[0] - mode[1] + mode[5] + mode[6] - mode[17] - mode[18] + 2.f*(mode[10] - mode[16]));
  n_b.vd[3][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/18.f * (mode[0] + mode[2] - mode[5] + mode[6] + mode[17] - mode[18] - 2.f*(mode[11] + mode[16]));
  n_b.vd[4][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/18.f * (mode[0] - mode[2] - mode[5] + mode[6] + mode[17] - mode[18] + 2.f*(mode[11] - mode[16]));
  n_b.vd[5][x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = 1.f/18.f * (mode[0] + mode[3] - 2.f*(mode[6] + mode[12] + mode[16] - mode[18]));
  n_b.vd[6][x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = 1.f/18.f * (mode[0] - mode[3] - 2.f*(mode[6] - mode[12] + mode[16] - mode[18]));
  n_b.vd[7][(x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/36.f * (mode[0] + mode[1] + mode[2] + mode[4] + 2.f*mode[6] + mode[7] + mode[10] + mode[11] + mode[13] + mode[14] + mode[16] + 2.f*mode[18]);
  n_b.vd[8][(para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/36.f * (mode[0] - mode[1] - mode[2] + mode[4] + 2.f*mode[6] + mode[7] - mode[10] - mode[11] - mode[13] - mode[14] + mode[16] + 2.f*mode[18]);
  n_b.vd[9][(x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/36.f * (mode[0] + mode[1] - mode[2] + mode[4] + 2.f*mode[6] - mode[7] + mode[10] - mode[11] + mode[13] - mode[14] + mode[16] + 2.f*mode[18]);
  n_b.vd[10][(para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = 1.f/36.f * (mode[0] - mode[1] + mode[2] + mode[4] + 2.f*mode[6] - mode[7] - mode[10] + mode[11] - mode[13] + mode[14] + mode[16] + 2.f*mode[18]);
  n_b.vd[11][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = 1.f/36.f * (mode[0] + mode[1] + mode[3] + mode[4] + mode[5] - mode[6] + mode[8] + mode[10] + mode[12] - mode[13] + mode[15] + mode[16] + mode[17] - mode[18]);
  n_b.vd[12][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = 1.f/36.f * (mode[0] - mode[1] - mode[3] + mode[4] + mode[5] - mode[6] + mode[8] - mode[10] - mode[12] + mode[13] - mode[15] + mode[16] + mode[17] - mode[18]);
  n_b.vd[13][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = 1.f/36.f * (mode[0] + mode[1] - mode[3] + mode[4] + mode[5] - mode[6] - mode[8] + mode[10] - mode[12] - mode[13] - mode[15] + mode[16] + mode[17] - mode[18]);
  n_b.vd[14][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = 1.f/36.f * (mode[0] - mode[1] + mode[3] + mode[4] + mode[5] - mode[6] - mode[8] - mode[10] + mode[12] + mode[13] + mode[15] + mode[16] + mode[17] - mode[18]);
  n_b.vd[15][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = 1.f/36.f * (mode[0] + mode[2] + mode[3] + mode[4] - mode[5] - mode[6] + mode[9] + mode[11] + mode[12] - mode[14] - mode[15] + mode[16] - mode[17] - mode[18]);
  n_b.vd[16][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = 1.f/36.f * (mode[0] - mode[2] - mode[3] + mode[4] - mode[5] - mode[6] + mode[9] - mode[11] - mode[12] + mode[14] + mode[15] + mode[16] - mode[17] - mode[18]);
  n_b.vd[17][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = 1.f/36.f * (mode[0] + mode[2] - mode[3] + mode[4] - mode[5] - mode[6] - mode[9] + mode[11] - mode[12] - mode[14] + mode[15] + mode[16] - mode[17] - mode[18]);
  n_b.vd[18][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = 1.f/36.f * (mode[0] - mode[2] + mode[3] + mode[4] - mode[5] - mode[6] - mode[9] - mode[11] + mode[12] + mode[14] - mode[15] + mode[16] - mode[17] - mode[18]);

}
/*-------------------------------------------------------*/
/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
/*@{
 * @param index			node index / thread index (Input)
 * @param n_b			Pointer to local node residing in array b (Input)
 * @param n_a			Pointer to local node residing in array a (Output) (temp stored in buffer a)
}*/
/*-------------------------------------------------------*/
__device__ void bounce_back_read(LB_nodes_gpu n_b, LB_nodes_gpu n_a, unsigned int index){
    
  unsigned int xyz[3];

  if(n_b.boundary[index] == 1){
    index_to_xyz(index, xyz);
    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];

    /* stream vd from boundary node back to origin node */
    n_a.vd[1][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = n_b.vd[2][index];
    n_a.vd[2][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = n_b.vd[1][index];
    n_a.vd[3][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[4][index];
    n_a.vd[4][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[3][index];
    n_a.vd[5][x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_b.vd[6][index];
    n_a.vd[6][x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_b.vd[5][index];
    n_a.vd[7][(x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[8][index];
    n_a.vd[8][(para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[7][index];
    n_a.vd[9][(x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[10][index];
    n_a.vd[10][(para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_b.vd[9][index];
    n_a.vd[11][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_b.vd[12][index];
    n_a.vd[12][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_b.vd[11][index]; 
    n_a.vd[13][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_b.vd[14][index]; 
    n_a.vd[14][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_b.vd[13][index]; 
    n_a.vd[15][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_b.vd[16][index];
    n_a.vd[16][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_b.vd[15][index];
    n_a.vd[17][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_b.vd[18][index]; 
    n_a.vd[18][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_b.vd[17][index];
  }
}
/*-------------------------------------------------------*/
/*@{
 * @param index			node index / thread index (Input)
 * @param n_b			Pointer to local node residing in array b (Input)
 * @param n_a			Pointer to local node residing in array a (Output) (temp stored in buffer a)
}*/
/*-------------------------------------------------------*/
__device__ void bounce_back_write(LB_nodes_gpu n_b, LB_nodes_gpu n_a, unsigned int index){

  unsigned int xyz[3];

  if(n_b.boundary[index] == 1){
    index_to_xyz(index, xyz);
    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];

    n_b.vd[1][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = n_a.vd[1][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z];
    n_b.vd[2][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z] = n_a.vd[2][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z];
    n_b.vd[3][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[3][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[4][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[4][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[5][x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_a.vd[5][x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)];
    n_b.vd[6][x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_a.vd[6][x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)];
    n_b.vd[7][(x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[7][(x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[8][(para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[8][(para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[9][(x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[9][(x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[10][(para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z] = n_a.vd[10][(para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z];
    n_b.vd[11][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_a.vd[11][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)];
    n_b.vd[12][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_a.vd[12][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)];
    n_b.vd[13][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_a.vd[13][(x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)];
    n_b.vd[14][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_a.vd[14][(para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z)];
    n_b.vd[15][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_a.vd[15][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)];
    n_b.vd[16][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_a.vd[16][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)];
    n_b.vd[17][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] = n_a.vd[17][x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)];
    n_b.vd[18][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)] = n_a.vd[18][x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z)];
  }
}

/*-------------------------------------------------------*/
/**used for reset the field to create qickanddirty non periodic bc !!only for advanced users!!*/
/*@{
 * @param index			node index / thread index (Input)
 * @param n_a			Pointer to local node residing in array a (Input)
 * @param *n_b			Pointer to local node residing in array b (Output)
}*/
/*-------------------------------------------------------*/
__device__ void reset_pop(LB_nodes_gpu n_b, LB_nodes_gpu n_a, unsigned int index){	

  //float avg_rho = para.rho*para.agrid*para.agrid*para.agrid;
  /* delete populations */
#if 1
  if(n_b.boundary[index] == 2){
    n_b.vd[0][index] = n_b.vd[0][index+1];
    n_b.vd[1][index] = n_b.vd[1][index+1];
    n_b.vd[2][index] = n_b.vd[2][index+1];
    n_b.vd[3][index] = n_b.vd[3][index+1];
    n_b.vd[4][index] = n_b.vd[4][index+1];
    n_b.vd[5][index] = n_b.vd[5][index+1];
    n_b.vd[6][index] = n_b.vd[6][index+1];
    n_b.vd[7][index]  = n_b.vd[7][index+1];
    n_b.vd[8][index]  = n_b.vd[8][index+1];
    n_b.vd[9][index]  = n_b.vd[9][index+1];
    n_b.vd[10][index] = n_b.vd[10][index+1];
    n_b.vd[11][index] = n_b.vd[11][index+1];
    n_b.vd[12][index] = n_b.vd[12][index+1];
    n_b.vd[13][index] = n_b.vd[13][index+1];
    n_b.vd[14][index] = n_b.vd[14][index+1];
    n_b.vd[15][index] = n_b.vd[15][index+1];
    n_b.vd[16][index] = n_b.vd[16][index+1];
    n_b.vd[17][index] = n_b.vd[17][index+1];
    n_b.vd[18][index] = n_b.vd[18][index+1];
  }
#endif
}
/*-------------------------------------------------------*/
/** add of (external) forces within the modespace, needed for particle-interaction */
/*@{
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Input/Output)
 * @param node_f	Pointer to local node force (Input)
}*/
/*-------------------------------------------------------*/
__device__ void apply_forces(unsigned int index, float *mode, LB_node_force_gpu node_f) {

  float Rho, u[3], C[6];
  Rho = mode[0] + para.rho*para.agrid*para.agrid*para.agrid;

  /* hydrodynamic momentum density is redefined when forces present */
  u[0] = (mode[1] + 0.5f*node_f.force[0][index])/Rho;
  u[1] = (mode[2] + 0.5f*node_f.force[1][index])/Rho;
  u[2] = (mode[3] + 0.5f*node_f.force[2][index])/Rho;

  C[0] = (1.f + para.gamma_bulk)*u[0]*node_f.force[0][index] + 1.f/3.f*(para.gamma_bulk-para.gamma_shear)*(u[0]*node_f.force[0][index] + u[1]*node_f.force[1][index] + u[2]*node_f.force[2][index]);
  C[2] = (1.f + para.gamma_bulk)*u[1]*node_f.force[1][index] + 1.f/3.f*(para.gamma_bulk-para.gamma_shear)*(u[0]*node_f.force[0][index] + u[1]*node_f.force[1][index] + u[2]*node_f.force[2][index]);
  C[5] = (1.f + para.gamma_bulk)*u[2]*node_f.force[2][index] + 1.f/3.f*(para.gamma_bulk-para.gamma_shear)*(u[0]*node_f.force[0][index] + u[1]*node_f.force[1][index] + u[2]*node_f.force[2][index]);
  C[1] = 1.f/2.f*(1.f+para.gamma_shear)*(u[0]*node_f.force[1][index]+u[1]*node_f.force[0][index]);
  C[3] = 1.f/2.f*(1.f+para.gamma_shear)*(u[0]*node_f.force[2][index]+u[2]*node_f.force[0][index]);
  C[4] = 1.f/2.f*(1.f+para.gamma_shear)*(u[1]*node_f.force[2][index]+u[2]*node_f.force[1][index]);

  /* update momentum modes */
  mode[1] += node_f.force[0][index];
  mode[2] += node_f.force[1][index];
  mode[3] += node_f.force[2][index];
  	
  /* update stress modes */
  mode[4] += C[0] + C[2] + C[5];
  mode[5] += C[0] - C[2];
  mode[6] += C[0] + C[2] - 2.f*C[5];
  mode[7] += C[1];
  mode[8] += C[3];
  mode[9] += C[4];

#ifdef EXTERNAL_FORCES
  if(para.external_force){
    node_f.force[0][index] = para.ext_force[0]*para.agrid*para.agrid*para.tau*para.tau;
    node_f.force[1][index] = para.ext_force[1]*para.agrid*para.agrid*para.tau*para.tau;
    node_f.force[2][index] = para.ext_force[2]*para.agrid*para.agrid*para.tau*para.tau;
  }
#else
  /* reset force */
  node_f.force[0][index] = 0.f;
  node_f.force[1][index] = 0.f;
  node_f.force[2][index] = 0.f;
#endif
}
/*-------------------------------------------------------*/
/**function used to calc physical values of every node*/
/*@{
 * @param index		node index / thread index (Input)
 * @param mode		Pointer to the local register values mode (Input)
 * @param n_a		Pointer to local node residing in array a for boundary flag(Input)
 * @param *d_v		Pointer to local device values (Input/Output)
}*/
/*-------------------------------------------------------*/
__device__ void calc_values(LB_nodes_gpu n_a, float *mode, LB_values_gpu *d_v, unsigned int index, unsigned int singlenode){

  float Rho = mode[0] + para.rho*para.agrid*para.agrid*para.agrid;
	
  /**implemented due to the problem of division via zero*/
  if(n_a.boundary[index] == 1){
    Rho = 1.0f;
    mode[1] = 0.f;
    mode[2] = 0.f;
    mode[3] = 0.f;
  }

  if(singlenode == 1){
    d_v[0].rho = Rho;
    d_v[0].v[0] = mode[1]/Rho;
    d_v[0].v[1] = mode[2]/Rho;
    d_v[0].v[2] = mode[3]/Rho;
  }
  else{
    d_v[index].rho = Rho;
    d_v[index].v[0] = mode[1]/Rho;
    d_v[index].v[1] = mode[2]/Rho;
    d_v[index].v[2] = mode[3]/Rho;
  }
#if 0
  if(singlenode == 1){
    /* equilibrium part of the stress modes */
    /**to print out the stress tensor entries, ensure that in lbgpu.h struct the values are available*/
    d_v[0].pi[0] = ((mode[1]*mode[1]) + (mode[2]*mode[2]) + (mode[3]*mode[3]))/para.rho;
    d_v[0].pi[1] = ((mode[1]*mode[1]) - (mode[2]*mode[2]))/para.rho;
    d_v[0].pi[2] = ((mode[1]*mode[1]) + (mode[2]*mode[2])  + (mode[3]*mode[3])) - 3.0f*(mode[3]*mode[3]))/para.rho;
    d_v[0].pi[3] = mode[1]*mode[2]/para.rho;
    d_v[0].pi[4] = mode[1]*mode[3]/para.rho;
    d_v[0].pi[5] = mode[2]*mode[3]/para.rho;
  else{
    d_v[index].pi[0] = ((mode[1]*mode[1]) + (mode[2]*mode[2]) + (mode[3]*mode[3]))/para.rho;
    d_v[index].pi[1] = ((mode[1]*mode[1]) - (mode[2]*mode[2]))/para.rho;
    d_v[index].pi[2] = ((mode[1]*mode[1]) + (mode[2]*mode[2])  + (mode[3]*mode[3])) - 3.0f*(mode[3]*mode[3]))/para.rho;
    d_v[index].pi[3] = mode[1]*mode[2]/para.rho;
    d_v[index].pi[4] = mode[1]*mode[3]/para.rho;
    d_v[index].pi[5] = mode[2]*mode[3]/para.rho;
  }
#endif
}
/*-------------------------------------------------------*/
/*@{
 * @param node_index	node index around (8) particle (Input)
 * @param *mode			Pointer to the local register values mode (Output)
 * @param n_a			Pointer to local node residing in array a(Input)
}*/
/*-------------------------------------------------------*/
__device__ void calc_mode(float *mode, LB_nodes_gpu n_a, unsigned int node_index){
	
  /* mass mode */
  mode[0] = n_a.vd[0][node_index] + n_a.vd[1][node_index] + n_a.vd[2][node_index]
          + n_a.vd[3][node_index] + n_a.vd[4][node_index] + n_a.vd[5][node_index]
          + n_a.vd[6][node_index] + n_a.vd[7][node_index] + n_a.vd[8][node_index]
          + n_a.vd[9][node_index] + n_a.vd[10][node_index] + n_a.vd[11][node_index] + n_a.vd[12][node_index]
          + n_a.vd[13][node_index] + n_a.vd[14][node_index] + n_a.vd[15][node_index] + n_a.vd[16][node_index]
          + n_a.vd[17][node_index] + n_a.vd[18][node_index];

  /* momentum modes */
  mode[1] = (n_a.vd[1][node_index] - n_a.vd[2][node_index]) + (n_a.vd[7][node_index] - n_a.vd[8][node_index])
          + (n_a.vd[9][node_index] - n_a.vd[10][node_index]) + (n_a.vd[11][node_index] - n_a.vd[12][node_index])
          + (n_a.vd[13][node_index] - n_a.vd[14][node_index]);
  mode[2] = (n_a.vd[3][node_index] - n_a.vd[4][node_index]) + (n_a.vd[7][node_index] - n_a.vd[8][node_index])
          - (n_a.vd[9][node_index] - n_a.vd[10][node_index]) + (n_a.vd[15][node_index] - n_a.vd[16][node_index])
          + (n_a.vd[17][node_index] - n_a.vd[18][node_index]);
  mode[3] = (n_a.vd[5][node_index] - n_a.vd[6][node_index]) + (n_a.vd[11][node_index] - n_a.vd[12][node_index])
          - (n_a.vd[13][node_index] - n_a.vd[14][node_index]) + (n_a.vd[15][node_index] - n_a.vd[16][node_index])
          - (n_a.vd[17][node_index] - n_a.vd[18][node_index]);
}
/*********************************************************/
/** \name Coupling part */
/*********************************************************/
/*-------------------------------------------------------*/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))*/
/*@{
 * @param n_a				Pointer to local node residing in array a (Input)
 * @param *delta			Pointer for the weighting of particle position (Output)
 * @param *delta_j			Pointer for the weighting of particle momentum (Output)
 * @param *particle_data	Pointer to the particle position and velocity (Input)
 * @param *particle_force	Pointer to the particle force (Input)
 * @param part_index		particle id / thread id (Input)
 * @param *rn				Pointer to randomnumber array of the particle
 * @param node_index		node index around (8) particle (Output)
}*/
/*-------------------------------------------------------*/
__device__ void calc_viscous_force(LB_nodes_gpu n_a, float *delta, LB_particle_gpu *particle_data, LB_particle_force_gpu *particle_force, unsigned int part_index, LB_randomnr_gpu *rn_part, float *delta_j, unsigned int *node_index){
	
  float mode[4];
  unsigned int my_left[3];
  float interpolated_u1, interpolated_u2, interpolated_u3;
  float Rho;
  interpolated_u1 = interpolated_u2 = interpolated_u3 = 0.f;

  float temp_delta[6];
  float temp_delta_half[6];

/** see ahlrichs + duennweg page 8227 equ (10) and (11) */
  #pragma unroll
  for(int i=0; i<3; ++i){
    float scaledpos = particle_data[part_index].p[i]/para.agrid;
    my_left[i] = (unsigned int)(floorf(scaledpos));
    temp_delta[3+i] = scaledpos - my_left[i];
    temp_delta[i] = 1.f - temp_delta[3+i];
    /**further value used for interpolation of fluid velocity at part pos near boundaries */
    temp_delta_half[3+i] = (scaledpos - my_left[i])*2.f;
    temp_delta_half[i] = 2.f - temp_delta_half[3+i];
  }

  delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
  delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
  delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
  delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
  delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
  delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
  delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
  delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

  unsigned int x = my_left[0];
  unsigned int y = my_left[1];
  unsigned int z = my_left[2];

  node_index[0] = x                + para.dim_x*y                  + para.dim_x*para.dim_y*z;
  node_index[1] = (x+1)%para.dim_x + para.dim_x*y                  + para.dim_x*para.dim_y*z;
  node_index[2] = x                + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  node_index[4] = x                + para.dim_x*y                  + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[5] = (x+1)%para.dim_x + para.dim_x*y                  + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[6] = x                + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
	
#if 0
	/** calc of the interpolated verlocity at the position of the particle !!!still under investigation and development!!!*/
  if(n_a.boundary[node_index[0]] == 1){
    delta[1] = temp_delta_half[3] * temp_delta[1] * temp_delta[2];
    delta[2] = temp_delta[0] * temp_delta_half[4] * temp_delta[2];
    delta[4] = temp_delta[0] * temp_delta[1] * temp_delta_half[5];
  }
  if(n_a.boundary[node_index[1]] == 1){		
    delta[0] = temp_delta_half[0] * temp_delta[1] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta_half[4] * temp_delta[2];
    delta[5] = temp_delta[3] * temp_delta[1] * temp_delta_half[5];
  }
  if(n_a.boundary[node_index[2]] == 1){		
    delta[0] = temp_delta_half[0] * temp_delta[1] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta_half[4] * temp_delta[2];
    delta[6] = temp_delta[0] * temp_delta[4] * temp_delta_half[5];
  }
  if(n_a.boundary[node_index[3]] == 1){		
    delta[1] = temp_delta[3] * temp_delta_half[1] * temp_delta[2];
    delta[2] = temp_delta_half[0] * temp_delta[4] * temp_delta[2];
    delta[7] = temp_delta[3] * temp_delta[4] * temp_delta_half[5];
  }
  if(n_a.boundary[node_index[4]] == 1){		
    delta[0] = temp_delta[0] * temp_delta[1] * temp_delta_half[2];
    delta[5] = temp_delta_half[3] * temp_delta[1] * temp_delta[5];
    delta[6] = temp_delta[0] * temp_delta_half[4] * temp_delta[5];
  }
  if(n_a.boundary[node_index[5]] == 1){		
    delta[1] = temp_delta[3] * temp_delta[1] * temp_delta_half[2];
    delta[4] = temp_delta_half[0] * temp_delta[1] * temp_delta[5];
    delta[7] = temp_delta[3] * temp_delta_half[4] * temp_delta[5];
  }
  if(n_a.boundary[node_index[6]] == 1){		
    delta[2] = temp_delta[0] * temp_delta[4] * temp_delta_half[2];
    delta[4] = temp_delta[0] * temp_delta_half[1] * temp_delta[5];
    delta[7] = temp_delta_half[3] * temp_delta[4] * temp_delta[5];
  }
  if(n_a.boundary[node_index[7]] == 1){		
    delta[3] = temp_delta[3] * temp_delta[4] * temp_delta_half[2];
    delta[5] = temp_delta[3] * temp_delta_half[1] * temp_delta[5];
    delta[6] = temp_delta_half[0] * temp_delta[4] * temp_delta[5];
  }

  if(n_a.boundary[node_index[0]] == 1)delta[0] = 0.f;

  if(n_a.boundary[node_index[1]] == 1)delta[1] = 0.f;

  if(n_a.boundary[node_index[2]] == 1)delta[2] = 0.f;

  if(n_a.boundary[node_index[3]] == 1)delta[3] = 0.f;

  if(n_a.boundary[node_index[4]] == 1)delta[4] = 0.f;

  if(n_a.boundary[node_index[5]] == 1)delta[5] = 0.f;

  if(n_a.boundary[node_index[6]] == 1)delta[6] = 0.f;

  if(n_a.boundary[node_index[7]] == 1)delta[7] = 0.f;
#endif

 #pragma unroll
  for(int i=0; i<8; ++i){
    calc_mode(mode, n_a, node_index[i]);
    Rho = mode[0] + para.rho*para.agrid*para.agrid*para.agrid;	
    interpolated_u1 += delta[i]*mode[1]/(Rho);
    interpolated_u2 += delta[i]*mode[2]/(Rho);
    interpolated_u3 += delta[i]*mode[3]/(Rho);
  }

	/* calculate viscous force
	* take care to rescale velocities with time_step and transform to MD units
 	* (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
#ifdef LB_ELECTROHYDRODYNAMICS
  particle_force[part_index].f[0] = - para.friction * (particle_data[part_index].v[0]/para.time_step - interpolated_u1*para.agrid/para.tau - particle_data[part_index].mu_E[0]);
  particle_force[part_index].f[1] = - para.friction * (particle_data[part_index].v[1]/para.time_step - interpolated_u2*para.agrid/para.tau - particle_data[part_index].mu_E[1]);
  particle_force[part_index].f[2] = - para.friction * (particle_data[part_index].v[2]/para.time_step - interpolated_u3*para.agrid/para.tau - particle_data[part_index].mu_E[2]);
#else
  particle_force[part_index].f[0] = - para.friction * (particle_data[part_index].v[0]/para.time_step - interpolated_u1*para.agrid/para.tau);
  particle_force[part_index].f[1] = - para.friction * (particle_data[part_index].v[1]/para.time_step - interpolated_u2*para.agrid/para.tau);
  particle_force[part_index].f[2] = - para.friction * (particle_data[part_index].v[2]/para.time_step - interpolated_u3*para.agrid/para.tau);
#endif
	/** add stochastik force of zero mean (Ahlrichs, Duennweg equ. 15)*/
#ifdef GAUSSRANDOM
  gaussian_random(rn_part);
  particle_force[part_index].f[0] += para.lb_coupl_pref2*rn_part->randomnr[0];
  particle_force[part_index].f[1] += para.lb_coupl_pref2*rn_part->randomnr[1];
  gaussian_random(rn_part);
  particle_force[part_index].f[2] += para.lb_coupl_pref2*rn_part->randomnr[0];
#else
  random_01(rn_part);
  particle_force[part_index].f[0] += para.lb_coupl_pref*(rn_part->randomnr[0]-0.5f);
  particle_force[part_index].f[1] += para.lb_coupl_pref*(rn_part->randomnr[1]-0.5f);
  random_01(rn_part);
  particle_force[part_index].f[2] += para.lb_coupl_pref*(rn_part->randomnr[0]-0.5f);
#endif
	  
  /* delta_j for transform momentum transfer to lattice units which is done in calc_node_force
  (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  delta_j[0] = - particle_force[part_index].f[0]*para.time_step*para.tau/para.agrid;
  delta_j[1] = - particle_force[part_index].f[1]*para.time_step*para.tau/para.agrid;
  delta_j[2] = - particle_force[part_index].f[2]*para.time_step*para.tau/para.agrid;  	
															  																	  
}
/*-------------------------------------------------------*/
/**calcutlation of the node force caused by the particles, with atomicadd due to avoiding race conditions 
	(Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))*/
/*@{

 * @param *delta			Pointer for the weighting of particle position (Input)
 * @param *delta_j			Pointer for the weighting of particle momentum (Input)
 * @param node_index		node index around (8) particle (Input)
 * @param node_f    		Pointer to the node force (Output).
}*/
/*-------------------------------------------------------*/
__device__ void calc_node_force(float *delta, float *delta_j, unsigned int *node_index, LB_node_force_gpu node_f){

  atomicadd(&(node_f.force[0][node_index[0]]), (delta[0]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[0]]), (delta[0]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[0]]), (delta[0]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[1]]), (delta[1]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[1]]), (delta[1]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[1]]), (delta[1]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[2]]), (delta[2]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[2]]), (delta[2]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[2]]), (delta[2]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[3]]), (delta[3]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[3]]), (delta[3]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[3]]), (delta[3]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[4]]), (delta[4]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[4]]), (delta[4]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[4]]), (delta[4]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[5]]), (delta[5]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[5]]), (delta[5]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[5]]), (delta[5]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[6]]), (delta[6]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[6]]), (delta[6]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[6]]), (delta[6]*delta_j[2]));

  atomicadd(&(node_f.force[0][node_index[7]]), (delta[7]*delta_j[0]));
  atomicadd(&(node_f.force[1][node_index[7]]), (delta[7]*delta_j[1]));
  atomicadd(&(node_f.force[2][node_index[7]]), (delta[7]*delta_j[2]));

}
/*-------------------------------------------------------*/
/**additional check if the particles are within the box */
/*@{
 * @param *particle_data	Pointer to the particle data (Input).
 * @param part_index		index of the particle == thread index (Input).		
}*/
/*-------------------------------------------------------*/
__device__ void check_part_posis(LB_particle_gpu *particle_data, unsigned int part_index){

  if(particle_data[part_index].p[0]/para.agrid < 0.f || particle_data[part_index].p[0]/para.agrid > para.dim_x){
    printf("particle out of box! (dim_x) \t %u \t %f \n", part_index, particle_data[part_index].p[0]); 
  }
  if(particle_data[part_index].p[1]/para.agrid < 0.f || particle_data[part_index].p[1]/para.agrid > para.dim_y){
    printf("particle out of box! (dim_y) \t %u \t %f \n", part_index, particle_data[part_index].p[1]); 
  }
  if(particle_data[part_index].p[2]/para.agrid < 0.f || particle_data[part_index].p[2]/para.agrid > para.dim_z){
    printf("particle out of box! (dim_z) \t %u \t %f \n", part_index, particle_data[part_index].p[2]); 
  }
}
/*-------------------------------------------------------*/
/**kernel to calculate local populations from hydrodynamic fields given by the tcl values.
 *
 * The mapping is given in terms of the equilibrium distribution.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 * Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 * @param n_a		 Pointer to the lattice site (Input).
 * @param node_f    Pointer to the node force (Input).
 * @param *gpu_check additional check if gpu kernel are executed(Input).
 */
/*-------------------------------------------------------*/
__global__ void calc_n_equilibrium(LB_nodes_gpu n_a, LB_node_force_gpu node_f, int *gpu_check) {

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){

    /*temp gesetzt aus lb_reinit_fluid() wären Anfangs-Werte die aus tcl übergeben werden*/
    /* default values for fields in lattice units */
    gpu_check[0] = 1;

    float Rho = para.rho*para.agrid*para.agrid*para.agrid;
    float v[3] = { 0.0f, 0.0f, 0.0f };
    float pi[6] = { Rho*c_sound_sq, 0.0f, Rho*c_sound_sq, 0.0f, 0.0f, Rho*c_sound_sq };

    float rhoc_sq = Rho*c_sound_sq;
    float avg_rho = para.rho*para.agrid*para.agrid*para.agrid;
    float local_rho, local_j[3], *local_pi, trace;

    local_rho  = Rho;

    local_j[0] = Rho * v[0];
    local_j[1] = Rho * v[1];
    local_j[2] = Rho * v[2];

    local_pi = pi;

    /* reduce the pressure tensor to the part needed here */
    local_pi[0] -= rhoc_sq;
    local_pi[2] -= rhoc_sq;
    local_pi[5] -= rhoc_sq;

    trace = local_pi[0] + local_pi[2] + local_pi[5];

    float rho_times_coeff;
    float tmp1,tmp2;

    /* update the q=0 sublattice */
    n_a.vd[0][index] = 1.f/3.f * (local_rho-avg_rho) - 1.f/2.f*trace;

    /* update the q=1 sublattice */
    rho_times_coeff = 1.f/18.f * (local_rho-avg_rho);

    n_a.vd[1][index] = rho_times_coeff + 1.f/6.f*local_j[0] + 1.f/4.f*local_pi[0] - 1.f/12.f*trace;
    n_a.vd[2][index] = rho_times_coeff - 1.f/6.f*local_j[0] + 1.f/4.f*local_pi[0] - 1.f/12.f*trace;
    n_a.vd[3][index] = rho_times_coeff + 1.f/6.f*local_j[1] + 1.f/4.f*local_pi[2] - 1.f/12.f*trace;
    n_a.vd[4][index] = rho_times_coeff - 1.f/6.f*local_j[1] + 1.f/4.f*local_pi[2] - 1.f/12.f*trace;
    n_a.vd[5][index] = rho_times_coeff + 1.f/6.f*local_j[2] + 1.f/4.f*local_pi[5] - 1.f/12.f*trace;
    n_a.vd[6][index] = rho_times_coeff - 1.f/6.f*local_j[2] + 1.f/4.f*local_pi[5] - 1.f/12.f*trace;

    /* update the q=2 sublattice */
    rho_times_coeff = 1.f/36.f * (local_rho-avg_rho);

    tmp1 = local_pi[0] + local_pi[2];
    tmp2 = 2.0f*local_pi[1];
    n_a.vd[7][index]  = rho_times_coeff + 1.f/12.f*(local_j[0]+local_j[1]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[8][index]  = rho_times_coeff - 1.f/12.f*(local_j[0]+local_j[1]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[9][index]  = rho_times_coeff + 1.f/12.f*(local_j[0]-local_j[1]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;
    n_a.vd[10][index] = rho_times_coeff - 1.f/12.f*(local_j[0]-local_j[1]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;

    tmp1 = local_pi[0] + local_pi[5];
    tmp2 = 2.0f*local_pi[3];

    n_a.vd[11][index] = rho_times_coeff + 1.f/12.f*(local_j[0]+local_j[2]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[12][index] = rho_times_coeff - 1.f/12.f*(local_j[0]+local_j[2]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[13][index] = rho_times_coeff + 1.f/12.f*(local_j[0]-local_j[2]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;
    n_a.vd[14][index] = rho_times_coeff - 1.f/12.f*(local_j[0]-local_j[2]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;

    tmp1 = local_pi[2] + local_pi[5];
    tmp2 = 2.0f*local_pi[4];

    n_a.vd[15][index] = rho_times_coeff + 1.f/12.f*(local_j[1]+local_j[2]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[16][index] = rho_times_coeff - 1.f/12.f*(local_j[1]+local_j[2]) + 1.f/8.f*(tmp1+tmp2) - 1.f/24.f*trace;
    n_a.vd[17][index] = rho_times_coeff + 1.f/12.f*(local_j[1]-local_j[2]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;
    n_a.vd[18][index] = rho_times_coeff - 1.f/12.f*(local_j[1]-local_j[2]) + 1.f/8.f*(tmp1-tmp2) - 1.f/24.f*trace;

    /*set different seed for randomgen on every node */
    n_a.seed[index] = para.your_seed + index;
  }
}
/*-------------------------------------------------------*/
/** kernel for the initalisation of the particle force array*/
/*@{
 * @param *particle_force	Pointer to local particle force (Output)
 * @param *part				Pointer to the particle rn seed storearray (Output)
}*/
/*-------------------------------------------------------*/
__global__ void init_particle_force(LB_particle_force_gpu *particle_force, LB_particle_seed_gpu *part){
	
  unsigned int part_index = blockDim.x * blockIdx.x + threadIdx.x;
	
  if(part_index<para.number_of_particles){
    particle_force[part_index].f[0] = 0.0f;
    particle_force[part_index].f[1] = 0.0f;
    particle_force[part_index].f[2] = 0.0f;
	
    part[part_index].seed = para.your_seed + part_index;
  }
			
}
/*-------------------------------------------------------*/
/** kernel for the initalisation of the partikel force array */
/*@{
 * @param *particle_force	pointer to local particle force (Input)
}*/
/*-------------------------------------------------------*/
__global__ void reset_particle_force(LB_particle_force_gpu *particle_force){
	
  unsigned int part_index = blockDim.x * blockIdx.x + threadIdx.x;
	
  if(part_index<para.number_of_particles){
    particle_force[part_index].f[0] = 0.0f;
    particle_force[part_index].f[1] = 0.0f;
    particle_force[part_index].f[2] = 0.0f;
  }			
}
/*-------------------------------------------------------*/
/** (re-)initialization of the node force / set up of external force in lb units */
/*@{
 * @param node_f		Pointer to local node force (Input)
}*/
/*-------------------------------------------------------*/
__global__ void reinit_node_force(LB_node_force_gpu node_f){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    if(para.external_force){
      node_f.force[0][index] = para.ext_force[0]*para.agrid*para.agrid*para.tau*para.tau;
      node_f.force[1][index] = para.ext_force[1]*para.agrid*para.agrid*para.tau*para.tau;
      node_f.force[2][index] = para.ext_force[2]*para.agrid*para.agrid*para.tau*para.tau;
    }
    else{
      node_f.force[0][index] = 0.0f;
      node_f.force[1][index] = 0.0f;
      node_f.force[2][index] = 0.0f;
    }
  }
}
#if 1
/*-------------------------------------------------------*/
/**hard coded boundary kernel for custom made boundaries */
/**just for advanced LB users to setup special boundaries or mark some nodes with
	the boundary flag e.g. to reset this nodes*/
/*@{
 * @param n_a		Pointer to local node residing in array a (Input/Output)
 * @param n_b		Pointer to local node residing in array b (Input/Output)
}*/
/*-------------------------------------------------------*/
__global__ void init_boundaries_hardcoded(LB_nodes_gpu n_a, LB_nodes_gpu n_b){
	
  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){	
#if 1
    unsigned int xyz[3];
    index_to_xyz(index, xyz);
    unsigned int x = xyz[0];
#endif
#if 0
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];

    n_a.boundary[index] = n_b.boundary[index] = 0;
	
    /* bottomplate || topplate*/	
    if(index < para.dim_x*para.dim_y || index >= (para.dim_z-1)*para.dim_x*para.dim_y){
      n_a.boundary[index] = n_b.boundary[index] = 1;	
    }
#endif
#if 0
    if(x == (para.dim_x/8) && y > ((35*para.dim_y)/100)  && y < ((65*para.dim_y)/100) && z > (25*para.dim_z/100) && z < (75*para.dim_z/100)){
      n_a.boundary[index] = n_b.boundary[index] = 1;
    }
#endif	
#if 1
    if(x == 0){
      n_a.boundary[index] = n_b.boundary[index] = 2;
    }
#endif	
  }
}
#endif
/*-------------------------------------------------------*/
/**set the boundary flag for all boundary nodes */
/*@{
 * @param *boundindex	Pointer to the 1d index of the boundnode (Input)
 * @param n_a			Pointer to local node residing in array a (Input)
 * @param n_b			Pointer to local node residing in array b (Input)
}*/
/*-------------------------------------------------------*/
__global__ void init_boundaries(int *boundindex, int number_of_boundnodes, LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<number_of_boundnodes){
    n_a.boundary[boundindex[index]] = n_b.boundary[boundindex[index]] = 1;
  }	
}
/*-------------------------------------------------------*/
/**reset the boundary flag of every node */
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param n_b		Pointer to local node residing in array b (Input)	
}*/
/*-------------------------------------------------------*/
__global__ void reset_boundaries(LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    n_a.boundary[index] = n_b.boundary[index] = 0;
  }
}
/*-------------------------------------------------------*/
/** integrationstep of the lb-fluid-solver */
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param n_b		Pointer to local node residing in array b (Input)
 * @param *d_v		Pointer to local device values (Input)
 * @param node_f	Pointer to local node force (Input)
}*/
/*-------------------------------------------------------*/
__global__ void integrate(LB_nodes_gpu n_a, LB_nodes_gpu n_b, LB_values_gpu *d_v, LB_node_force_gpu node_f){
    
  /**every node is connected to a thread via the index*/
  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
  /**the 19 moments (modes) are only temporary register values */
  float mode[19];
  LB_randomnr_gpu rng;

  if(index<para.number_of_nodes){
    /** storing the seed into a register value*/
    rng.seed = n_a.seed[index];
    /**calc_m_from_n*/
    calc_m_from_n(n_a, index, mode);
    /**lb_relax_modes*/
    relax_modes(mode, index, node_f);
    /**lb_thermalize_modes */
    if (para.fluct) thermalize_modes(mode, index, &rng);
#ifdef EXTERNAL_FORCES
    /**if external force is used apply node force */
    if (para.external_force) apply_forces(index, mode, node_f);
#else
    /**if partcles are used apply node forces*/
    if (para.number_of_particles) apply_forces(index, mode, node_f); 
#endif
    /**lb_calc_n_from_modes_push*/
    normalize_modes(mode);
    /**calc of velocity densities and streaming with pbc*/
    calc_n_from_modes_push(n_b, mode, index);
    /** rewriting the seed back to the global memory*/
    n_b.seed[index] = rng.seed;
  }  
}
/*-------------------------------------------------------*/
/** part interaction kernel */
/*@{
 * @param n_a				Pointer to local node residing in array a (Input)
 * @param *particle_data	Pointer to the particle position and velocity (Input)
 * @param *particle_force	Pointer to the particle force (Input)
 * @param *part				Pointer to the rn array of the particles (Input)
 * @param node_f			Pointer to local node force (Input)
}*/
/*-------------------------------------------------------*/
__global__ void calc_fluid_particle_ia(LB_nodes_gpu n_a, LB_particle_gpu *particle_data, LB_particle_force_gpu *particle_force, LB_node_force_gpu node_f, LB_particle_seed_gpu *part){
	
  unsigned int part_index = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[8];
  float delta[8];
  float delta_j[3];
  LB_randomnr_gpu rng_part;
	
  if(part_index<para.number_of_particles){
#if 0		
    /** check if particles are in the box*/
    check_part_posis(particle_data, part_index);		
#endif
    rng_part.seed = part[part_index].seed;
    /**calc of the force which act on the particle */
    calc_viscous_force(n_a, delta, particle_data, particle_force, part_index, &rng_part, delta_j, node_index);
    /**calc of the force which acts back to the fluid node */
    calc_node_force(delta, delta_j, node_index, node_f);
    part[part_index].seed = rng_part.seed;		
  }
}
/*-------------------------------------------------------*/
/**Bounce back boundary read kernel*/
/*@{
 * @param n_a					Pointer to local node residing in array a (Input)
 * @param n_b					Pointer to local node residing in array b (Input)
}*/
/*-------------------------------------------------------*/
__global__ void bb_read(LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    bounce_back_read(n_b, n_a, index);
  }
}
/*-------------------------------------------------------*/
/**Bounce back boundary write kernel*/
/*@{
 * @param n_a					Pointer to local node residing in array a (Input)
 * @param n_b					Pointer to local node residing in array b (Input)
}*/
/*-------------------------------------------------------*/
__global__ void bb_write(LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    bounce_back_write(n_b, n_a, index);
  }
}
/*-------------------------------------------------------*/
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param n_b		Pointer to local node residing in array b (Input)
}*/
/*-------------------------------------------------------*/
__global__ void reset_population(LB_nodes_gpu n_b, LB_nodes_gpu n_a){
 	
  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    reset_pop(n_b, n_a, index);
  }
}
/*-------------------------------------------------------*/
/*@{
 * @param n_a		Pointer to local node residing in array a (Input)
 * @param *d_v		Pointer to local device values (Input)
}*/
/*-------------------------------------------------------*/
__global__ void values(LB_nodes_gpu n_a, LB_values_gpu *d_v){

  float mode[19];
  unsigned int singlenode = 0;
  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes){
    calc_mode(mode, n_a, index);
    calc_values(n_a, mode, d_v, index, singlenode);
  }
}

__global__ void init_extern_nodeforces(int n_extern_nodeforces, LB_extern_nodeforce_gpu *extern_nodeforces, LB_node_force_gpu node_f){

  unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

  if(index<n_extern_nodeforces){
    node_f.force[0][extern_nodeforces[index].index] = extern_nodeforces[index].force[0]*para.agrid*para.agrid*para.tau*para.tau;
    node_f.force[1][extern_nodeforces[index].index] = extern_nodeforces[index].force[1]*para.agrid*para.agrid*para.tau*para.tau;
    node_f.force[2][extern_nodeforces[index].index] = extern_nodeforces[index].force[2]*para.agrid*para.agrid*para.tau*para.tau;
  }
}
__global__ void lb_print_node(int single_nodeindex, LB_values_gpu *d_p_v, LB_nodes_gpu n_a){
	
  float mode[19];
  unsigned int singlenode = 1;

  if(blockDim.x * blockIdx.x + threadIdx.x == 0){
    calc_mode(mode, n_a, single_nodeindex);
    calc_values(n_a, mode, d_p_v, single_nodeindex, singlenode);
  }	
}

void cuda_safe_mem(cudaError_t err){
    if( cudaSuccess != err) {                                             
      fprintf(stderr, "Could not allocate gpu memory.\n");
      exit(EXIT_FAILURE);
    }
}
void cuda_safe_kernel(cudaError_t err){
    if( cudaSuccess != err) {                                             
      fprintf(stderr, "cuda kernel failed! (maybe your dimensions are to large).\n");
      exit(EXIT_FAILURE);
    }
}
/**********************************************************************/
/* Host funktions to setup and call kernels*/
/**********************************************************************/
/**-------------------------------------------------------*/
/*@{
 * @param *lbpar_gpu	Pointer to parameters to setup the lb field
}*/
/**-------------------------------------------------------*/
void lb_init_GPU(LB_parameters_gpu *lbpar_gpu){

  // Allocate lattice-struct in device memory
  size_of_values = lbpar_gpu->number_of_nodes * sizeof(LB_values_gpu);
  size_of_forces = lbpar_gpu->number_of_particles * sizeof(LB_particle_force_gpu);
  size_of_positions = lbpar_gpu->number_of_particles * sizeof(LB_particle_gpu);
  size_of_seed = lbpar_gpu->number_of_particles * sizeof(LB_particle_seed_gpu);

  cuda_safe_mem(cudaMalloc((void**)&device_values, size_of_values));

  for(int i=0; i<19; i++){
    cuda_safe_mem(cudaMalloc((void**)&nodes_a.vd[i], lbpar_gpu->number_of_nodes * sizeof(float)));
    cuda_safe_mem(cudaMalloc((void**)&nodes_b.vd[i], lbpar_gpu->number_of_nodes * sizeof(float)));                                           
  }
  cuda_safe_mem(cudaMalloc((void**)&nodes_a.seed, lbpar_gpu->number_of_nodes * sizeof(unsigned int)));
  cuda_safe_mem(cudaMalloc((void**)&nodes_a.boundary, lbpar_gpu->number_of_nodes * sizeof(unsigned int)));
  cuda_safe_mem(cudaMalloc((void**)&nodes_b.seed, lbpar_gpu->number_of_nodes * sizeof(unsigned int)));
  cuda_safe_mem(cudaMalloc((void**)&nodes_b.boundary, lbpar_gpu->number_of_nodes * sizeof(unsigned int)));

  for(int i=0; i<3; i++){
    cuda_safe_mem(cudaMalloc((void**)&node_f.force[i], lbpar_gpu->number_of_nodes * sizeof(float)));
  }
  
  cuda_safe_mem(cudaMalloc((void**)&particle_force, size_of_forces));
  cuda_safe_mem(cudaMalloc((void**)&particle_data, size_of_positions));	
  cuda_safe_mem(cudaMalloc((void**)&part, size_of_seed));
	
  /**write parameters in const memory*/
  cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu));
  /**check flag if lb gpu init works*/
  cuda_safe_mem(cudaMalloc((void**)&gpu_check, sizeof(int)));
  h_gpu_check = (int*)malloc(sizeof(int));

  /** values for the kernel call */
  threads_per_block = 128;
  blocks_per_grid = (lbpar_gpu->number_of_nodes + threads_per_block - 1) /(threads_per_block);

  /** values for the particle kernel */
  threads_per_block_particles = 64;
  blocks_per_grid_particles = (lbpar_gpu->number_of_particles + threads_per_block_particles - 1)/(threads_per_block_particles);

  reset_boundaries<<<blocks_per_grid, threads_per_block>>>(nodes_a, nodes_b);
  cuda_safe_kernel(cudaGetLastError());	
  /** calc of veloctiydensities from given parameters and initialize the Node_Force array with zero */
  calc_n_equilibrium<<<blocks_per_grid, threads_per_block>>>(nodes_a, node_f, gpu_check);
  cuda_safe_kernel(cudaGetLastError());	
  /** init part forces with zero*/
  if(lbpar_gpu->number_of_particles) init_particle_force<<<blocks_per_grid_particles, threads_per_block_particles>>>(particle_force, part);
  cuda_safe_kernel(cudaGetLastError());	
  reinit_node_force<<<blocks_per_grid, threads_per_block>>>(node_f);
  cuda_safe_kernel(cudaGetLastError());	

  cudaStreamCreate(&stream[0]);

  h_gpu_check[0] = 0;
  cudaMemcpy(h_gpu_check, gpu_check, sizeof(int), cudaMemcpyDeviceToHost);

  cudaThreadSynchronize();

  if(!h_gpu_check[0]){
    fprintf(stderr, "initialization of lb gpu code failed! \n");
    errexit();	
  }	
}
/**-------------------------------------------------------------------------*/
/**setup and call particle reallocation from the host */
/*@{
 * @param *lbpar_gpu	Pointer to parameters to setup the lb field
}*/
/**-------------------------------------------------------------------------*/
void lb_realloc_particle_GPU(LB_parameters_gpu *lbpar_gpu){

  cudaFree(particle_force);
  cudaFree(particle_data);
  cudaFree(part);

  cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu));

  size_of_forces = lbpar_gpu->number_of_particles * sizeof(LB_particle_force_gpu);
  size_of_positions = lbpar_gpu->number_of_particles * sizeof(LB_particle_gpu);
  size_of_seed = lbpar_gpu->number_of_particles * sizeof(LB_particle_seed_gpu);

  cuda_safe_mem(cudaMalloc((void**)&particle_force, size_of_forces));
  cuda_safe_mem(cudaMalloc((void**)&particle_data, size_of_positions));
  cuda_safe_mem(cudaMalloc((void**)&part, size_of_seed));

  /** values for the particle kernel */
  threads_per_block_particles = 64;
  blocks_per_grid_particles = (lbpar_gpu->number_of_particles + threads_per_block_particles - 1)/(threads_per_block_particles);

  if(lbpar_gpu->number_of_particles) init_particle_force<<<blocks_per_grid_particles, threads_per_block_particles>>>(particle_force, part);
	
  if(lbpar_gpu->number_of_particles) reinit_node_force<<<blocks_per_grid, threads_per_block>>>(node_f);
	
}

/**-------------------------------------------------------------------------*/
/**setup and call boundaries from the host */
/*@{
 * @param *host_boundindex		Pointer to the host bound index
 * @param number_of_boundnodes	number of boundnodes
}*/
/**-------------------------------------------------------------------------*/
void lb_init_boundaries_GPU(int number_of_boundnodes, int *host_boundindex){

  size_of_boundindex = number_of_boundnodes*sizeof(int);
  cudaMemcpyToSymbol(number_of_bnodes, &number_of_boundnodes, sizeof(int));
  cudaMalloc((void**)&boundindex, size_of_boundindex);
  cudaMemcpy(boundindex, host_boundindex, size_of_boundindex, cudaMemcpyHostToDevice);

  reset_boundaries<<<blocks_per_grid, threads_per_block>>>(nodes_a, nodes_b);
  cuda_safe_kernel(cudaGetLastError());	

  threads_per_block_bound = 128;
  blocks_per_grid_bound = (number_of_boundnodes + threads_per_block_bound -1)/(threads_per_block_bound);

#if 0
  init_boundaries_hardcoded<<<blocks_per_grid_bound, threads_per_block_bound>>>(nodes_a, nodes_b);
#endif
  init_boundaries<<<blocks_per_grid_bound, threads_per_block_bound>>>(boundindex, number_of_boundnodes, nodes_a, nodes_b);
  cuda_safe_kernel(cudaGetLastError());	
  calc_n_equilibrium<<<blocks_per_grid, threads_per_block>>>(nodes_a, node_f, gpu_check);
  cuda_safe_kernel(cudaGetLastError());	
  cudaThreadSynchronize();
}
void lb_init_extern_nodeforces_GPU(int n_extern_nodeforces, LB_extern_nodeforce_gpu *host_extern_nodeforces, LB_parameters_gpu *lbpar_gpu){

  size_of_extern_nodeforces = n_extern_nodeforces*sizeof(LB_extern_nodeforce_gpu);
  cudaMalloc((void**)&extern_nodeforces, size_of_extern_nodeforces);
  cudaMemcpy(extern_nodeforces, host_extern_nodeforces, size_of_extern_nodeforces, cudaMemcpyHostToDevice);

  if(para.external_force == 0)cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)); 

  threads_per_block_exf = 128;
  blocks_per_grid_exf = (n_extern_nodeforces + threads_per_block_exf -1)/(threads_per_block_exf);
	
  init_extern_nodeforces<<<blocks_per_grid_exf, threads_per_block_exf>>>(n_extern_nodeforces, extern_nodeforces, node_f);
  cuda_safe_kernel(cudaGetLastError());	
}

/**-------------------------------------------------------------------------*/
/**setup and call particle kernel from the host */
/*@{
 * @param **host_data		Pointer to the host particle positions and velocities
}*/
/**-------------------------------------------------------------------------*/
void lb_particle_GPU(LB_particle_gpu *host_data){
  	
  /** get espresso md particle values*/
  cudaMemcpyAsync(particle_data, host_data, size_of_positions, cudaMemcpyHostToDevice, stream[0]);

  /** call of the particle kernel */
  calc_fluid_particle_ia<<<blocks_per_grid_particles, threads_per_block_particles, 0, stream[0]>>>(nodes_a, particle_data, particle_force, node_f, part);
  cuda_safe_kernel(cudaGetLastError());	
}
/** setup and call kernel to copy particle forces to host */
void lb_copy_forces_GPU(LB_particle_force_gpu *host_forces){

  /** Copy result from device memory to host memory*/
  cudaMemcpy(host_forces, particle_force, size_of_forces, cudaMemcpyDeviceToHost);

  /** reset part forces with zero*/
  reset_particle_force<<<blocks_per_grid_particles, threads_per_block_particles, 0,  stream[0]>>>(particle_force);
  cuda_safe_kernel(cudaGetLastError());	
  cudaThreadSynchronize();
}

/** setup and call kernel for getting macroscopic fluid values of all nodes*/
void lb_get_values_GPU(LB_values_gpu *host_values){

  values<<<blocks_per_grid, threads_per_block>>>(nodes_a, device_values);
  cuda_safe_kernel(cudaGetLastError());	

  cudaMemcpy(host_values, device_values, size_of_values, cudaMemcpyDeviceToHost);

  cudaThreadSynchronize();

}
/** setup and call kernel for getting macroscopic fluid values of a single node*/
void lb_print_node_GPU(int single_nodeindex, LB_values_gpu *host_print_values){ 
      
  LB_values_gpu *device_print_values;
  cudaMalloc((void**)&device_print_values, sizeof(LB_values_gpu));	
  threads_per_block_print = 1;
  blocks_per_grid_print = 1;
  lb_print_node<<<blocks_per_grid_print, threads_per_block_print>>>(single_nodeindex, device_print_values, nodes_a);
  cudaMemcpy(host_print_values, device_print_values, sizeof(LB_values_gpu), cudaMemcpyDeviceToHost);
  cuda_safe_kernel(cudaGetLastError());

  cudaThreadSynchronize();
}

/**-------------------------------------------------------------------------*/
			/**setup and call integrate kernel from the host */
/**-------------------------------------------------------------------------*/
void lb_integrate_GPU(){
		
  /**call of fluid step*/
  if (intflag == 1){
    (integrate<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_a, nodes_b, device_values, node_f));
    cuda_safe_kernel(cudaGetLastError());
#if 0		
    reset_population<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_b, nodes_a);
    cuda_safe_kernel(cudaGetLastError());
#endif
#ifdef LB_BOUNDARIES_GPU		
    if (lb_boundaries_bb_gpu == 1) bb_read<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_a, nodes_b);
      cuda_safe_kernel(cudaGetLastError());			
    if (lb_boundaries_bb_gpu == 1) bb_write<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_a, nodes_b);
      cuda_safe_kernel(cudaGetLastError());
#endif
    intflag = 0;
  }
  else{
    integrate<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_b, nodes_a, device_values, node_f);
    cuda_safe_kernel(cudaGetLastError());
#if 0		
    reset_population<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_a, nodes_b);
    cuda_safe_kernel(cudaGetLastError());
#endif
#ifdef LB_BOUNDARIES_GPU		
    if (lb_boundaries_bb_gpu == 1) bb_read<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_b, nodes_a);
      cuda_safe_kernel(cudaGetLastError());
			
    if (lb_boundaries_bb_gpu == 1) bb_write<<<blocks_per_grid, threads_per_block, 0,  stream[0]>>>(nodes_b, nodes_a);
      cuda_safe_kernel(cudaGetLastError());
#endif
    intflag = 1;
  }             
}

/**-------------------------------------------------------------------------*/
			/** free gpu memory kernel called from the host */
/**-------------------------------------------------------------------------*/
void lb_free_GPU(){
  // Free device memory
  cudaFree(device_values);
  cudaFree(&para);
  cudaFree(&nodes_a);
  cudaFree(&nodes_b);
  cudaFree(particle_force);
  cudaFree(particle_data);
  cudaFree(&node_f);
  cudaFree(part);
  cudaFree(&number_of_bnodes);
  cudaStreamDestroy(stream[0]);
}
#endif /* LB_GPU */
