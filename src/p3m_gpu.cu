/* 
   Copyright (C) 2010,2011,2012,2013 The ESPResSo project

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

/** \file lbgpu.cu
 *
 * Cuda (.cu) file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.h.
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include <cufft.h>
#include "cuda_common.h"
#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "p3m_gpu.h"
#include "utils.h"
#ifdef __cplusplus
}
#endif

#ifdef ELECTROSTATICS

struct dummytypename {
  CUFFT_TYPE_COMPLEX *charge_mesh;
  CUFFT_TYPE_COMPLEX *force_mesh;
  REAL_TYPE *G_hat, *G_hat_host;
  cufftHandle fft_plan;
  int cao, mesh;
  REAL_TYPE alpha;
  int npart;
  REAL_TYPE box;
} p3m_gpu_data;


#define SQR(A) ((A)*(A))

void static Aliasing_sums_ik ( int cao, REAL_TYPE box, REAL_TYPE alpha, int mesh, int NX, int NY, int NZ,
                        REAL_TYPE *Zaehler, REAL_TYPE *Nenner ) {
    REAL_TYPE S1,S2,S3;
    REAL_TYPE fak1,fak2,zwi;
    int    MX,MY,MZ;
    REAL_TYPE NMX,NMY,NMZ;
    REAL_TYPE NM2;
    REAL_TYPE expo, TE;
    REAL_TYPE Leni = 1.0/box;

    fak1 = 1.0/ ( REAL_TYPE ) mesh;
    fak2 = SQR ( PI/ ( alpha ) );

    Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

    for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
      NMX = ( ( NX > mesh/2 ) ? NX - mesh : NX ) + mesh*MX;
      S1 = pow ( sinc(fak1*NMX ), 2*cao );
      for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
	NMY = ( ( NY > mesh/2 ) ? NY - mesh : NY ) + mesh*MY;
	S2   = S1*pow ( sinc (fak1*NMY ), 2*cao );
	for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
	  NMZ = ( ( NZ > mesh/2 ) ? NZ - mesh : NZ ) + mesh*MZ;
	  S3   = S2*pow ( sinc( fak1*NMZ ), 2*cao );

	  NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
	  *Nenner += S3;

	  expo = fak2*NM2;
	  TE = exp ( -expo );
	  zwi  = S3 * TE/NM2;
	  Zaehler[0] += NMX*zwi*Leni;
	  Zaehler[1] += NMY*zwi*Leni;
	  Zaehler[2] += NMZ*zwi*Leni;
	}
      }
    }
}

/* Calculate influence function */
void static calculate_influence_function ( int cao, int mesh, REAL_TYPE box, REAL_TYPE alpha, REAL_TYPE *G_hat ) {

  int    NX,NY,NZ;
  REAL_TYPE Dnx,Dny,Dnz;
  REAL_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner=0.0;
  REAL_TYPE zwi;
  int ind = 0;
  REAL_TYPE Leni = 1.0/box;

  for ( NX=0; NX<mesh; NX++ ) {
    for ( NY=0; NY<mesh; NY++ ) {
      for ( NZ=0; NZ<mesh; NZ++ ) {
	ind = NX*mesh*mesh + NY * mesh + NZ;
	  
	if ( ( NX==0 ) && ( NY==0 ) && ( NZ==0 ) )
	  G_hat[ind]=0.0;
	else if ( ( NX% ( mesh/2 ) == 0 ) && ( NY% ( mesh/2 ) == 0 ) && ( NZ% ( mesh/2 ) == 0 ) )
	  G_hat[ind]=0.0;
	else {
	  Aliasing_sums_ik ( cao, box, alpha, mesh, NX, NY, NZ, Zaehler, &Nenner );
		  
	  Dnx = ( NX > mesh/2 ) ? NX - mesh : NX;
	  Dny = ( NY > mesh/2 ) ? NY - mesh : NY;
	  Dnz = ( NZ > mesh/2 ) ? NZ - mesh : NZ;
	    
	  zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
	  zwi /= ( ( SQR ( Dnx*Leni ) + SQR ( Dny*Leni ) + SQR ( Dnz*Leni ) ) * SQR ( Nenner ) );
	  G_hat[ind] = 2.0 * zwi / PI;
	}
      }
    }
  }
}

//NOTE :if one wants to use the function below it requires cuda compute capability 1.3
#ifdef _P3M_GPU_REAL_DOUBLE
__device__ double atomicAdd (double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

/** atomic add function for several cuda architectures 
*/

#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
//atomicAdd supports floats already, do nothing
#elif __CUDA_ARCH__ >= 110
#warning Using slower atomicAdd emulation
__device__ inline void atomicAdd(float* address, float value){
  // float-atomic-add from 
// [url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]
  float old = value;
  while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
}
#else
#error I need at least compute capability 1.1
#endif



__device__ unsigned int getThreadIndexP3M() { //rename is dumb but can't import same fnc from cuda_common

  return blockIdx.y * gridDim.x * blockDim.x +
         blockDim.x * blockIdx.x +
         threadIdx.x;
}


// __global__ void add_p3m_farfield_force_gpu( LB_parameters_gpu* lb_parameters_gpu,
//                                             CUDA_particle_data* lb_particle_gpu,
//                                             CUDA_particle_force* lb_particle_force_gpu
//                                           ) {

//   unsigned int index = getThreadIndex();

//   if( index < lb_parameters_gpu->number_of_particles ) {
    
//     lb_particle_force_gpu[ index ].f[0] = 1.0f;
//     lb_particle_force_gpu[ index ].f[1] = 2.0f;
//     lb_particle_force_gpu[ index ].f[2] = 3.0f;
//   }
// }


template<int dim>
__global__ void apply_diff_op( CUFFT_TYPE_COMPLEX *mesh, const int mesh_size, CUFFT_TYPE_COMPLEX *force_mesh,  const REAL_TYPE box ) {
  int linear_index = mesh_size*mesh_size*blockIdx.x + mesh_size * blockIdx.y + threadIdx.x;
  int n;

  switch( dim ) {
  case 0:
    n = blockIdx.x;
    break;
  case 1:
    n = blockIdx.y;
    break;
  case 2:
    n = threadIdx.x;
    break;
  }

  n = ( n == mesh_size/2 ) ? 0.0 : n;
  n = ( n > mesh_size/2) ? n - mesh_size : n;
 
  force_mesh[linear_index].x =  -2.0 * PI * n * mesh[linear_index].y / box;
  force_mesh[linear_index].y =   2.0 * PI * n * mesh[linear_index].x / box;
}


__device__ inline int wrap_index(const int ind, const int mesh) {
  if(ind < 0)
    return ind + mesh;
  else if(ind >= mesh)
    return ind - mesh;
  else 
    return ind;	   
}

__device__ REAL_TYPE caf(int i, REAL_TYPE x, int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
      return 0.0;
    }
  }
  case 5 : {
    switch (i) {
    case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
    case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
    case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
    case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
    case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    default:
      return 0.0;
    }
  }
  case 6 : {
    switch (i) {
    case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
    case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
    case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
    case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
    case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
    case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    default:
      return 0.0;
    }
  }
  case 7 : {
    switch (i) {
    case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
    case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
    case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
    case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
    case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
    case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
    case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    default:
      return 0.0;
    }
  }
  }}
  return 0.0;
}

__global__ void apply_influence_function( CUFFT_TYPE_COMPLEX *mesh, int mesh_size, REAL_TYPE *G_hat ) {
  int linear_index = mesh_size*mesh_size*blockIdx.x + mesh_size * blockIdx.y + threadIdx.x;
  mesh[linear_index].x *= G_hat[linear_index];
  mesh[linear_index].y *= G_hat[linear_index];
}

__global__ void assign_charges(const CUDA_particle_data * const pdata,
CUFFT_TYPE_COMPLEX *mesh, const int m_size, const int cao, const REAL_TYPE pos_shift, const
REAL_TYPE hi) {
      /** id of the particle **/
      int id = blockIdx.x;
      /** position relative to the closest gird point **/
      REAL_TYPE m_pos[3];
      /** index of the nearest mesh point **/
      int nmp_x, nmp_y, nmp_z;      
      
      CUDA_particle_data p = pdata[id];

      m_pos[0] = p.p[0] * hi - pos_shift;
      m_pos[1] = p.p[1] * hi - pos_shift;
      m_pos[2] = p.p[2] * hi - pos_shift;

      nmp_x = (int) floor(m_pos[0] + 0.5);
      nmp_y = (int) floor(m_pos[1] + 0.5);
      nmp_z = (int) floor(m_pos[2] + 0.5);

      m_pos[0] -= nmp_x;
      m_pos[1] -= nmp_y;
      m_pos[2] -= nmp_z;

      nmp_x = wrap_index(nmp_x + threadIdx.x, m_size);
      nmp_y = wrap_index(nmp_y + threadIdx.y, m_size);
      nmp_z = wrap_index(nmp_z + threadIdx.z, m_size);

      atomicAdd( &(mesh[m_size*m_size*nmp_x +  m_size*nmp_y + nmp_z].x), caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*p.q);
}

__global__ void assign_forces(const CUDA_particle_data * const pdata, CUFFT_TYPE_COMPLEX *mesh, const int m_size, const int cao, const REAL_TYPE pos_shift, const
			      REAL_TYPE hi, CUDA_particle_force * lb_particle_force_gpu, REAL_TYPE prefactor, int dim) {
      /** id of the particle **/
      int id = blockIdx.x;
      /** position relative to the closest gird point **/
      REAL_TYPE m_pos[3];
      /** index of the nearest mesh point **/
      int nmp_x, nmp_y, nmp_z;      

      CUDA_particle_data p = pdata[id];

      m_pos[0] = p.p[0] * hi - pos_shift;
      m_pos[1] = p.p[1] * hi - pos_shift;
      m_pos[2] = p.p[2] * hi - pos_shift;

      nmp_x = (int) floor(m_pos[0] + 0.5);
      nmp_y = (int) floor(m_pos[1] + 0.5);
      nmp_z = (int) floor(m_pos[2] + 0.5);

      m_pos[0] -= nmp_x;
      m_pos[1] -= nmp_y;
      m_pos[2] -= nmp_z;

      nmp_x = wrap_index(nmp_x + threadIdx.x, m_size);
      nmp_y = wrap_index(nmp_y + threadIdx.y, m_size);
      nmp_z = wrap_index(nmp_z + threadIdx.z, m_size);

      atomicAdd( &(lb_particle_force_gpu[id].f[dim]), (float)(-prefactor*mesh[m_size*m_size*nmp_x +  m_size*nmp_y + nmp_z].x*caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*p.q));
      
}

extern "C" {

  /* Init the internal datastructures of the P3M GPU.
   * Mainly allocation on the device and influence function calculation.
   * Be advised: this needs mesh^3*5*sizeof(REAL_TYPE) of device memory. 
   */

  void p3m_gpu_init(int cao, int mesh, REAL_TYPE alpha, REAL_TYPE box) {
    gpu_init_particle_comm();
    
 
    if ( this_node == 0 ) {
      p3m_gpu_data.npart = gpu_get_global_particle_vars_pointer_host()->number_of_particles;
      p3m_gpu_data.alpha = alpha;
      p3m_gpu_data.cao = cao;
      p3m_gpu_data.mesh = mesh;
      p3m_gpu_data.box = box;
      int mesh3 = mesh*mesh*mesh;

      cudaMalloc((void **)&(p3m_gpu_data.charge_mesh), mesh3*sizeof(CUFFT_TYPE_COMPLEX));
      cudaMalloc((void **)&(p3m_gpu_data.force_mesh), mesh3*sizeof(CUFFT_TYPE_COMPLEX));
      cudaMalloc((void **)&(p3m_gpu_data.G_hat), mesh3*sizeof(REAL_TYPE));

      p3m_gpu_data.G_hat_host = (REAL_TYPE *)malloc(mesh3*sizeof(REAL_TYPE));

      // Calculate influence function of host.
      calculate_influence_function( cao, mesh, box, alpha, p3m_gpu_data.G_hat_host);

      // Copy influence function to device.
      cudaMemcpy( p3m_gpu_data.G_hat, p3m_gpu_data.G_hat_host, mesh3*sizeof(REAL_TYPE), cudaMemcpyHostToDevice);

      cufftPlan3d(&(p3m_gpu_data.fft_plan), mesh, mesh, mesh, CUFFT_PLAN_FLAG);
    }
  }

void p3m_gpu_add_farfield_force() {

  CUDA_particle_data* lb_particle_gpu;
  CUDA_particle_force* lb_particle_force_gpu;
  
  int mesh = p3m_gpu_data.mesh;
  int mesh3 = mesh*mesh*mesh;
  int cao = p3m_gpu_data.cao;
  REAL_TYPE box = p3m_gpu_data.box;

  lb_particle_gpu = gpu_get_particle_pointer();
  lb_particle_force_gpu = gpu_get_particle_force_pointer();

  p3m_gpu_data.npart = gpu_get_global_particle_vars_pointer_host()->number_of_particles;

  if(p3m_gpu_data.npart == 0)
    return;

  //printf("p3m params: mesh %d npart %d cao %d\n", mesh, p3m_gpu_data.npart, cao); //TODO delete

  dim3 gridAssignment(p3m_gpu_data.npart,1,1);
  dim3 threadsAssignment(cao,cao,cao);
  
  dim3 gridConv(mesh,mesh,1);
  dim3 threadsConv(mesh,1,1);

  REAL_TYPE pos_shift = (REAL_TYPE)((cao-1)/2);
  REAL_TYPE hi = mesh/box;
  REAL_TYPE prefactor = 1.0/(box*box*box*2.0);

  cudaMemset( p3m_gpu_data.charge_mesh, 0, mesh3*sizeof(CUFFT_TYPE_COMPLEX));

  KERNELCALL(assign_charges, gridAssignment, threadsAssignment, (lb_particle_gpu,p3m_gpu_data.charge_mesh,mesh,cao,pos_shift,hi));

  cudaThreadSynchronize();

  if (CUFFT_FFT(p3m_gpu_data.fft_plan, p3m_gpu_data.charge_mesh, p3m_gpu_data.charge_mesh, CUFFT_FORWARD) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed\n");
    return;
  }

  KERNELCALL( apply_influence_function, gridConv, threadsConv, (p3m_gpu_data.charge_mesh, mesh, p3m_gpu_data.G_hat));

  KERNELCALL(apply_diff_op<0>, gridConv, threadsConv, (p3m_gpu_data.charge_mesh, mesh, p3m_gpu_data.force_mesh, box));
  
  CUFFT_FFT(p3m_gpu_data.fft_plan, p3m_gpu_data.force_mesh, p3m_gpu_data.force_mesh, CUFFT_INVERSE);

  KERNELCALL(assign_forces, gridAssignment, threadsAssignment, (lb_particle_gpu, p3m_gpu_data.force_mesh, mesh, cao, pos_shift, hi, lb_particle_force_gpu, prefactor, 0));

  KERNELCALL(apply_diff_op<1>, gridConv, threadsConv, (p3m_gpu_data.charge_mesh, mesh, p3m_gpu_data.force_mesh, box));

  CUFFT_FFT(p3m_gpu_data.fft_plan, p3m_gpu_data.force_mesh, p3m_gpu_data.force_mesh, CUFFT_INVERSE);
  
  KERNELCALL(assign_forces, gridAssignment, threadsAssignment, (lb_particle_gpu, p3m_gpu_data.force_mesh, mesh, cao, pos_shift, hi, lb_particle_force_gpu, prefactor, 1));

  KERNELCALL(apply_diff_op<2>, gridConv, threadsConv, (p3m_gpu_data.charge_mesh, mesh, p3m_gpu_data.force_mesh, box));

  CUFFT_FFT(p3m_gpu_data.fft_plan, p3m_gpu_data.force_mesh, p3m_gpu_data.force_mesh, CUFFT_INVERSE);
  
  KERNELCALL(assign_forces, gridAssignment, threadsAssignment, (lb_particle_gpu, p3m_gpu_data.force_mesh, mesh, cao, pos_shift, hi, lb_particle_force_gpu, prefactor, 2));


  // KERNELCALL( add_p3m_farfield_force_gpu, dim_grid, threads_per_block, ( lb_parameters_gpu, lb_particle_gpu, lb_particle_force_gpu ) );
}

}

#endif /* ELECTROSTATICS */
