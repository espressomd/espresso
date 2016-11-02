/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "config.hpp"

#include <stdio.h>
#include <assert.h>
#include "cuda_runtime.h"
#include <device_functions.h>

#include "integrate_sd_cuda.hpp"
#include "integrate_sd_cuda_device.cu"

#ifdef SD

/* *************************************************************************************************************** *
 * ********************************************      CUDA-KERNELS     ******************************************** *
 * *************************************************************************************************************** */


/// This computes the farfield contribution of the mobility in the case of no
/// periodic boundary conditions.
/// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// N is the number of particles
/// self_mobility is 1./(6.*PI*eta*a)
/// a is the particle radius
/// mobility is the mobility matrix which will be retruned
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
__global__ void sd_compute_mobility_matrix(const real * r, int N, real self_mobility, real a, real * mobility){
  real mypos[3];
  const int lda=((3*N+31)/32)*32;
  __shared__ real cachedPos[3*numThreadsPerBlock];
  __shared__ real writeCache[3*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  // get data for myposition - using coalscaled memory access
  for (int l=0;l<3;l++){
    mydebug(" 0x%08x -> 0x%08x  \n",numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x,numThreadsPerBlock*l+threadIdx.x);
    cachedPos[numThreadsPerBlock*l+threadIdx.x] = r[numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x];
  }
  __syncthreads();
  for (int l=0;l<3;l++){
    mypos[l] = cachedPos[threadIdx.x*3+l];
    mydebug("mypos[%d]:  %e\n",l,mypos[l]);
  }

  /*if (i < N){
    // first write the self contribution
#pragma unroll
    for (int k=0; k < DIM; k++){
      //#pragma unroll
      //for (int l=0; l < DIM; l++){
      //mobility[myindex(DIM*i+k,DIM*i+l)]=0;
      //}
      mobility[myindex(DIM*i+k,DIM*i+k)]=self_mobility;
    }
    }*/
  for (int offset=0;offset<N;offset+=numThreadsPerBlock){
    mydebug("offset: %d\n",offset)
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      mydebug("fuu:: 0x%08x  0x%08x  0x%08x  0x%08x %e\n",r, offset*3,numThreadsPerBlock*l,threadIdx.x,cachedPos[numThreadsPerBlock*l+threadIdx.x]);
      cachedPos[numThreadsPerBlock*l+threadIdx.x] = r[offset*3+numThreadsPerBlock*l+threadIdx.x];
    }
    __syncthreads();
    if (i < N){
      for (int j=offset;j<min(offset+numThreadsPerBlock,N);j++){
	real dr[DIM];
	real dr2=0;
#pragma unroll 3
	for (int k=0;k<DIM;k++){
	  dr[k]=mypos[k]-cachedPos[DIM*(j-offset)+k]; // r_ij
	  dr2+=dr[k]*dr[k];
	}
	dr2=max(dr2,0.01);
	real drn= sqrt(dr2); // length of dr
	real b = a/drn;
      
	/*if (0.5 < b){  // drn < 2*a
	  /*real t=3./32./drn/a*self_mobility;
	  real t2=(1-9./32.*drn/a)*self_mobility;
	  for (k=0; k < DIM; k++){
	  for (l=0;l < DIM; l++){
	  mobility[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*t;
	  }
	  mobility[myindex(DIM*i+k,DIM*j+k)]+=t2;
	  }*/ // this should not happen ...
	real t,t2;
	// this also catches the case i == j
	if (0.5 < b ){  // drn < 2*a
	  t=0;
	  t2=0;
	  if (i==j){
	    t2=self_mobility;
	  }
	} else {
	  real b2=(a*a)/dr2;
	  // Rotne Prager
	  //t=(0.75-1.5*b2)*b/dr2*self_mobility;
	  //t2=(0.75+0.5*b2)*b*self_mobility;
	  
	  t=(0.75-1.5*b2)*b/dr2*self_mobility;
	  t2=(0.75+0.5*b2)*b*self_mobility;
	}
	//mobility[threadIdx.x]=3+threadIdx.x;
	real tmp_el13;
#pragma unroll 3
	for (int k=0; k < DIM; k++){
	  if (k ==0){ // these ifs should be removed at compile time ... after unrolling
#pragma unroll 3
	    for (int l=0;l < 3; l++){
	      //mobility[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*t;
	      writeCache[3*threadIdx.x+l]=dr[k]*dr[l]*t;
	    }
	  }
	  else if(k==1){
	    tmp_el13 = writeCache[3*threadIdx.x+2];
	    writeCache[3*threadIdx.x+0]=writeCache[3*threadIdx.x+1];
#pragma unroll 2
	    for (int l=1;l < DIM; l++){
	      //mobility[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*t;
	      writeCache[3*threadIdx.x+l]=dr[k]*dr[l]*t;
	    }	
	  }
	  else{
	    writeCache[3*threadIdx.x+0]=tmp_el13;
	    writeCache[3*threadIdx.x+1]=writeCache[3*threadIdx.x+2];
	    writeCache[3*threadIdx.x+2]=dr[k]*dr[2]*t;
	  }
	  writeCache[3*threadIdx.x+k]+=t2;
	    
	  __syncthreads();
	  //int max = min(blockDim.x, N-(blockIdx.x*blockDim.x));
	  int max = min(blockDim.x,N-blockDim.x*blockIdx.x);
	  for (int l=0;l<3;l++){
	    //mobility[(DIM*j+k)*3*N+blockIdx.x*blockDim.x+threadIdx.x+blockDim.x*l]=writeCache[threadIdx.x+blockDim.x*l];
	    mobility[(DIM*j+k)*lda+blockIdx.x*blockDim.x*3+max*l+threadIdx.x]=writeCache[max*l+threadIdx.x];
	  }
	  //mobility[myindex(DIM*i+k,DIM*j+k)]+=t2;
	}
	// python implementation:
	// T=one*(0.75+0.5*b2)*b+(0.75-1.5*b2)*b*drt*dr/dr2;
	//} // if (j <N)
      } // for (j = ...
    } // if (i < N)
  }// for offset = ...
}
#undef mydebug

/// This computes the farfield contribution of the mobility with ewald summation
/// this version assumes that the real-space cutoff is smaller than the boxlength
/// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// N is the number of particles
/// self_mobility is 1./(6.*PI*eta*a)
/// a is the particle radius
/// L_g is the boxlength
/// mobility is the mobility matrix which will be retruned
/// cutoff the realspace cutoff distance
/// xi the splitting parameter as defined by Beenakker 1986
/// xa  = xi * a
/// xa3 = xa * xa * xa
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
 __global__ void sd_compute_mobility_matrix_real_short(const real * r, int N, real self_mobility, real a, const real * L_g,
						       real * mobility, real cutoff, real xi, real xa, real xa3){
  real mypos[3];
  const int lda=((3*N+31)/32)*32;
  __shared__ real L[3];
  __shared__ real cachedPos[3*numThreadsPerBlock];
  __shared__ real writeCache[3*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
  __syncthreads();
  // get data for myposition - using coalscaled memory access
  for (int l=0;l<3;l++){
    cachedPos[numThreadsPerBlock*l+threadIdx.x] = r[numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x];
  }
  __syncthreads();
  for (int l=0;l<3;l++){
    mypos[l] = cachedPos[threadIdx.x*3+l];
    mydebug("mypos[%d]:  %e\n",l,mypos[l]);
  }

  for (int offset=0;offset<N;offset+=numThreadsPerBlock){
    mydebug("offset: %d\n",offset)
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      mydebug("fuu:: 0x%08x  0x%08x  0x%08x  0x%08x %e\n",r, offset*3,numThreadsPerBlock*l,threadIdx.x,cachedPos[numThreadsPerBlock*l+threadIdx.x]);
      cachedPos[numThreadsPerBlock*l+threadIdx.x] = r[offset*3+numThreadsPerBlock*l+threadIdx.x];
    }
    __syncthreads();
    if (i < N){
      for (int j=offset;j<min(offset+numThreadsPerBlock,N);j++){
	real dr[DIM];
	real dr2=0;
#pragma unroll 3
	for (int k=0;k<DIM;k++){
	  dr[k]=mypos[k]-cachedPos[DIM*(j-offset)+k]; // r_ij
	  dr[k]-=rint(dr[k]/L[k])*L[k]; // fold back
	  dr2+=dr[k]*dr[k];
	}
	dr2=max(dr2,0.01);
	real drn= sqrt(dr2); // length of dr
	real ar = a/drn;
	real xr = xi*drn;
	//real xa is given
	
	real t,t2;
	// this also catches the case i == j
	if (0.5 < ar || drn > cutoff){  // drn < 2*a
	  t=0;
	  t2=0;
	  if (i==j){
	    //t2=1-6./sqrt(M_PI)*xa+40./3./sqrt(M_PI)*xa*xa*xa;
#ifdef SD_USE_FLOAT
	    t2=1-3.385137501286537721688476709364635515064303775f*xa + 7.52252778063675049264105935414363447792067505f*xa3;
#else
	    t2=1-3.385137501286537721688476709364635515064303775*xa + 7.52252778063675049264105935414363447792067505*xa3;
#endif
	  }
	} else {
	  // Rotne Prager
	  real xr2=xr*xr;
	  real ar2=ar*ar;
#ifdef SD_USE_FLOAT
	  t=(0.75f-1.5f*ar2)*ar*erfcf(xr)+(-4.f*xa3*xr2*xr2-3.f*xa*xr2+16.f*xa3*xr2+1.5f*xa-2.f *xa3-3.f*xa*ar2)*0.5641895835477562869480794515607725858440506f*exp(-xr2);
	  t2=(0.75f+0.5f*ar2)*ar*erfcf(xr)+(4.f*xa3*xr2*xr2+3.f*xa*xr2-20.f*xa3*xr2-4.5f*xa+14.f*xa3+    xa*ar2)*0.5641895835477562869480794515607725858440506f*exp(-xr2);
#else
	  //assert(erfc(xr)>1e-15);
	  t=(0.75-1.5*ar2)*ar*erfc(xr);
	  //assert(t>1e-15);
	  real tmp1=(-4.*xa3*xr2*xr2-3.*xa*xr2+16.*xa3*xr2+1.5*xa-2. *xa3-3.*xa*ar2);
	  real tmp2=0.5641895835477562869480794515607725858440506*exp(-xr2);//
	  //assert(abs(tmp1)>1e-15);
	  //assert(abs(tmp2)>1e-15);
	  //assert(abs(tmp1*tmp2)>1e-15);
	  t+=tmp1*tmp2;
	  //assert(abs(t)>1e-15);
	  t2=(0.75+0.5*ar2)*ar*erfc(xr)+(4.*xa3*xr2*xr2+3.*xa*xr2-20.*xa3*xr2-4.5*xa+14.*xa3+   xa*ar2)*0.5641895835477562869480794515607725858440506*exp(-xr2);
#endif
	  //assert(t>1e-15);
	  //assert(t2>1e-15);
	}
	t*=self_mobility;
	t2*=self_mobility;
	t/=dr2;
	real tmp_el13;
#pragma unroll 3
	for (int k=0; k < DIM; k++){
	  if (k ==0){ // these ifs should be removed at compile time ... after unrolling
#pragma unroll 3
	    for (int l=0;l < 3; l++){
	      writeCache[3*threadIdx.x+l]=dr[k]*dr[l]*t;
	    }
	  }
	  else if(k==1){
	    tmp_el13 = writeCache[3*threadIdx.x+2];
	    writeCache[3*threadIdx.x+0]=writeCache[3*threadIdx.x+1];
#pragma unroll 2
	    for (int l=1;l < DIM; l++){
	      writeCache[3*threadIdx.x+l]=dr[k]*dr[l]*t;
	    }	
	  }
	  else{
	    writeCache[3*threadIdx.x+0]=tmp_el13;
	    writeCache[3*threadIdx.x+1]=writeCache[3*threadIdx.x+2];
	    writeCache[3*threadIdx.x+2]=dr[k]*dr[2]*t;
	  }
	  writeCache[3*threadIdx.x+k]+=t2;
	    
	  __syncthreads();
	  int max = min(blockDim.x,N-blockDim.x*blockIdx.x);
	  for (int l=0;l<3;l++){
	    mobility[(DIM*j+k)*lda+blockIdx.x*blockDim.x*3+max*l+threadIdx.x]=writeCache[max*l+threadIdx.x];
	  }
	}
      } // for (j = ...
    } // if (i < N)
  }// for offset = ...
}
#undef mydebug


/// This computes the farfield contribution of the mobility with ewald summation
/// this kernel computes the sines and cosinus.
/// r        : the vector of positions [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// N        : the number of particles
/// vecs *   : a pointer to the k-vectors
/// num      : number of k-vectors
/// sines    : the values for sin(pos_i k_j)
/// cosines  : the values for cos(pos_i k_j)
/// ldd      : the leading dimension of the sines and cosines array
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
__global__ void sd_compute_mobility_sines(const real * r, int N, const real * vecs, int num, real * sines, real * cosines, int ldd){
  real mypos[3];
  __shared__ real cache[3*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  // get data for myposition - using coalscaled memory access
  for (int l=0;l<3;l++){
    mydebug(" 0x%08x -> 0x%08x  \n",numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x,numThreadsPerBlock*l+threadIdx.x);
    cache[numThreadsPerBlock*l+threadIdx.x] = r[numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x];
  }
  __syncthreads();
  for (int l=0;l<3;l++){
    mypos[l] = cache[threadIdx.x*3+l];
    mydebug("mypos[%d]:  %e\n",l,mypos[l]);
  }
  
  for (int offset=0;offset<num;offset+=numThreadsPerBlock){
    mydebug("offset: %d\n",offset)
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      mydebug("fuu:: 0x%08x  0x%08x  0x%08x  0x%08x %e\n",r, offset*3,numThreadsPerBlock*l,threadIdx.x,cachedPos[numThreadsPerBlock*l+threadIdx.x]);
      cache[numThreadsPerBlock*l+threadIdx.x] = vecs[offset*3+numThreadsPerBlock*l+threadIdx.x];
    }
    __syncthreads();
    for (int j=offset;j<offset+numThreadsPerBlock&& j < num;j++){
      real vec_r=0;
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	vec_r+=mypos[k]*cache[3*(j-offset)+k];
      }
      real tmp;
      if (i < N){
	tmp=sin(vec_r);
      } else {
	tmp=0;
      }
      sines[ldd*j+i]  =tmp;
      if (i < N){
	tmp=cos(vec_r);
      }
      cosines[ldd*j+i]=tmp;
    } // for (j = ...
  }// for offset = ...
}
#undef mydebug



/// This computes the farfield contribution of the mobility with ewald summation
/// this kernel computes the sum in eq. (3.34) and (3.35) in the master thesis 
/// "rotne-prager based hydrodynamics" by D. Schwoerer: M^{wave} \sum F_i cos / sin.
/// The number of threads have to be a power of two e.g. 32, 64 or 128.
/// forces is the input vector of form [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
__global__ void sd_wavepart_sum(const real * const forces, const real * const vecs,    int num, const  real * const matrices_d,
				const real * const sines,  const real * const cosines, int ldd, int N, real * sin_sum, real * cos_sum, real max){
  // each block summs over all particles
  // with constant wave vektor
  real myforcesin[3]={0,0,0};
  real myforcecos[3]={0,0,0};
  extern __shared__ real cache[];//3*blockDim.x
  for (int offset=0;offset<N;offset+=blockDim.x){
    // for (int j=offset;j< offset_next;j++){ // is parallel
#pragma unroll 3
    for (int l=0;l<3;l++){
      cache[blockDim.x*l+threadIdx.x] = forces[blockDim.x*l+offset*3+threadIdx.x];
    }
    __syncthreads();
    real sine;
    real cosine;
    if (offset+threadIdx.x < N){
      sine   = sines  [blockIdx.x*ldd+offset+threadIdx.x];
      cosine = cosines[blockIdx.x*ldd+offset+threadIdx.x];
      //} else {
      //sine   = 0;
      //cosine = 0;
      //}
      //assert(!isnan(sine));
      //assert(!isnan(cosine));
      for (int k=0;k<3;k++){
	//assert(abs(cache[threadIdx.x*3+k]) <= max);
	//assert(abs(myforcesin[k]) <= max);
	//assert(abs(myforcecos[k]) <= max);
	//assert(!isnan(myforcesin[k]));
	//assert(!isnan(myforcecos[k]));
	real tmp=  sine*cache[threadIdx.x*3+k];
	//assert(abs(tmp) <= max);
	myforcesin[k]+=tmp;
	tmp     =cosine*cache[threadIdx.x*3+k];
	//assert(abs(tmp) <= max);
	myforcecos[k]+=tmp;
	//assert(!isnan(myforcesin[k]));
	//assert(!isnan(myforcecos[k]));
      } // for k
    } // if offset + threadIdx.x < N
  }// for offset = ...
  #pragma unroll 3
  for (int k=0 ; k < 3; k++){
    cache[threadIdx.x+blockDim.x*k]=myforcesin[k];
  }
  __syncthreads();
  reduce_sum(cache);
  //assert(cache[0] <= max);
  reduce_sum(cache+blockDim.x);
  //assert(cache[blockDim.x] <= max);
  reduce_sum(cache+blockDim.x*2);
  //assert(cache[blockDim.x*2] <= max);
  __syncthreads();
  if (threadIdx.x < 3){
    myforcesin[0]=cache[threadIdx.x*blockDim.x];
  }
  __syncthreads();
  for (int k=0 ; k < 3; k++){
    cache[threadIdx.x+blockDim.x*k]=myforcecos[k];
  }
  __syncthreads();
  reduce_sum(cache);
  //assert(cache[0] <= max);
  reduce_sum(cache+blockDim.x);
  //assert(cache[blockDim.x] <= max);
  reduce_sum(cache+blockDim.x*2);
  //assert(cache[blockDim.x*2] <= max);
  __syncthreads();
  if (threadIdx.x < 3){
    cache[threadIdx.x]   = cache[threadIdx.x*blockDim.x];
    cache[threadIdx.x+3] = myforcesin[0];
  }
  __syncthreads();
  real * mat = cache + 6;
  if (threadIdx.x < 6){
    //assert(cache[threadIdx.x] <= max);
    mat[threadIdx.x] = matrices_d[blockIdx.x*6+threadIdx.x];
    __syncthreads();
    cache[threadIdx.x+12]=cache[threadIdx.x]*mat[threadIdx.x%3];
  }
  __syncthreads();
  if (threadIdx.x < 2){
    cache[threadIdx.x*3+1+12]+=mat[3]*cache[threadIdx.x*3+0];
    cache[threadIdx.x*3+0+12]+=mat[3]*cache[threadIdx.x*3+1];
    cache[threadIdx.x*3+2+12]+=mat[4]*cache[threadIdx.x*3+0];
    cache[threadIdx.x*3+0+12]+=mat[4]*cache[threadIdx.x*3+2];
    cache[threadIdx.x*3+2+12]+=mat[5]*cache[threadIdx.x*3+1];
    cache[threadIdx.x*3+1+12]+=mat[5]*cache[threadIdx.x*3+2];
  }
  __syncthreads();
  if (threadIdx.x < 3){
    //assert(cache[threadIdx.x+12] <= max);
    //assert(cache[threadIdx.x+15] <= max);
    cos_sum[threadIdx.x+blockIdx.x*3]=cache[threadIdx.x+12];
    sin_sum[threadIdx.x+blockIdx.x*3]=cache[threadIdx.x+15];
  }
}



/// This computes the wave space contribution of the matrix vector product of
/// the farfield. This uses the pre computed values from the kernel
/// sd_wavepart_sum
/// num      : number of wave vectors
/// sines    : sinus of the scalar product of position and wave vector
/// cosines  : cosines of the scalar product of position and wave vector
/// sin_sum  : precomputed sum for sinus
/// cos_sum  : precomputed sum for cosinus
/// ldd      : leading dimension of the sin_sum and cos_sum arrays
/// out      : out vector upon which the result gets added
/// max      : maximal value of which the sin_sum and cos_sum can have
/// N        : number of particles
/// factor   : the scalar factor c of the matrix product c A \cdot x
__global__ void sd_wavepart_assemble(const int num, const real * const __restrict__ sines, const real * const __restrict__ cosines,
				     const real * const __restrict__ sin_sum, const real * const __restrict__ cos_sum, 
				     const int ldd, real * __restrict__ out, const real max, const int N, const real factor){
  // particle index of the particle for which we sum:
  int i = threadIdx.x+blockIdx.x*blockDim.x;
  real myout[3]={0,0,0};
  real sine;
  real cosine;
  
  int offset_next;
  for (int offset=0;offset<num;offset=offset_next){
    offset_next=offset+blockDim.x;
    if (offset_next > num){
      offset_next=num;
    }
    // j: wave vector index
    for (int j=offset;j< offset_next ;j++){
      sine   = sines  [i + j*ldd];
      //sine = read_without_caching(sines + i+j*ldd);
#pragma unroll 3
      for (int k=0;k < 3;k++){
	real tmp = sin_sum[k + j*3];
	//assert(tmp <= max);
	myout[k]+=sine*tmp;
      }
      cosine = cosines[i + j*ldd];
      //cosine = read_without_caching(cosines + i+j*ldd);
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	real tmp = cos_sum[k + j*3];
	//assert(tmp <= max);
	myout[k]+=cosine*tmp;
      }
    } // for (j = ...
  }// for offset = ...
  extern __shared__ real cache[];
#pragma unroll 3
  for (int k=0;k<3;k++){
    cache[threadIdx.x*3+k]=myout[k]*factor;;
  }
  __syncthreads();
#pragma unroll 3
  for (int k=0;k<3;k++){
    if (threadIdx.x+(blockIdx.x*3+k)*blockDim.x < 3*N){
      out[threadIdx.x+(blockIdx.x*3+k)*blockDim.x]+=cache[threadIdx.x + k*blockDim.x];
    }
    //cache[threadIdx.x+k*numThreadsPerBlock]+=out[threadIdx.x+(blockIdx.x*3+k)*numThreadsPerBlock];
    //out[threadIdx.x+(blockIdx.x*3+k)*numThreadsPerBlock]=cache[threadIdx.x+k*numThreadsPerBlock];
  }
}
/* this is slower (if we use a block per particle instead of a thread
__global__ void sd_wavepart_assemble_block(const int num, const real * const sines, const real * const cosines, const real * const sin_sum,
				     const real * const cos_sum, const int ldd, real * out, real max, int N, const real factor){
  // particle index of the particle for which we sum:
  int i = blockIdx.x;
  real myout[3]={0,0,0};
  real sine;
  real cosine;
  
  extern __shared__ real cache[];
  int offset_next;
  for (int offset=0;offset<num;offset=offset_next){
    offset_next=offset+blockDim.x;
    // j: wave vector index
    //for (int j=offset;j< offset_next ;j++){
    int j=offset+threadIdx.x;
    if (j < num){
      sine   = sines  [i + j*ldd];
#pragma unroll 3
      for (int k=0;k < 3;k++){
	real tmp = sin_sum[k + j*3];
	//assert(tmp <= max);
	myout[k]+=sine*tmp;
      }
      cosine = cosines[i + j*ldd];
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	real tmp = cos_sum[k + j*3];
	//assert(tmp <= max);
	myout[k]+=cosine*tmp;
      }
    } // for (j = ...
  }// for offset = ...
#pragma unroll 3
  for (int k=0;k<3;k++){
    cache[threadIdx.x]=myout[k]*factor;
    __syncthreads();
    reduce_sum(cache);
    __syncthreads();
    if (threadIdx.x == k){
      myout[0]=cache[0];
    }
  }
  if (threadIdx.x < 3){
    out[blockIdx.x * 3 + threadIdx.x]+= myout[0];
  }
}*/



/// Add the wavepart contribution of the mobility to the matrix
/// num      : number of wave space vectors
/// matrices_d : pointer to the matrices of the wave space contribution in the
///            ewald sum
/// sines    : sines of the scalar product of wave space vector and particle
///            positions
/// cosines  : cosines of the scalar product of wave space vector and
///            particle positions
/// ldd      : leading dimension of sines/cosines
/// N        : number of particles
/// mobility_d : mobility matrix to which the wave space contribution is added
/// mat_ldd  : leading dimension of the mobility matrix
__global__ void sd_wavepart_addto_matrix(const int num, const  real * const matrices_d, const real * const sines,  const real * const cosines,
					 const int ldd, const int N, real * const mobility_d, const int mat_ldd){
  const int i= threadIdx.x + blockDim.x*blockIdx.x;
  const int j= blockIdx.y;
  
  real mymat[6]={0,0,0,
		 0,0,0};
  if (i < N){
    for (int k = 0; k < num; k++){
      real mycos = cosines[i+k*ldd];
      mycos *= cosines[j+k*ldd];
      mycos += sines[i+k*ldd] * sines[j+k*ldd];
      for (int l=0;l<6;l++){
	mymat[l]+=matrices_d[6*k+l]*mycos;
      }
    }
  }
  extern __shared__ real share[];//3*blockDim.x
  #pragma unroll 3
  for (int d=0;d<3;d++){
    for (int l=0;l<3;l++){
      if (threadIdx.x+3*blockDim.x*blockIdx.x+l*blockDim.x<3*N)
	share[threadIdx.x+l*blockDim.x]=mobility_d[threadIdx.x+3*blockDim.x*blockIdx.x+l*blockDim.x+mat_ldd*(j*3+d)];
#ifdef SD_DEBUG
      else
	share[threadIdx.x+l*blockDim.x]=0;
#endif
    }
    __syncthreads();
    if (d==0){
      share[threadIdx.x*3+0]+=mymat[0];
      share[threadIdx.x*3+1]+=mymat[3];
      share[threadIdx.x*3+2]+=mymat[4];
    } else if (d==1){
      share[threadIdx.x*3+0]+=mymat[3];
      share[threadIdx.x*3+1]+=mymat[1];
      share[threadIdx.x*3+2]+=mymat[5];
    } else {
      share[threadIdx.x*3+0]+=mymat[4];
      share[threadIdx.x*3+1]+=mymat[5];
      share[threadIdx.x*3+2]+=mymat[2];
    }
    __syncthreads();
    for (int l=0;l<3;l++){
      if (threadIdx.x+3*blockDim.x*blockIdx.x+l*blockDim.x<3*N)
	mobility_d[threadIdx.x+3*blockDim.x*blockIdx.x+l*blockDim.x+mat_ldd*(j*3+d)]=share[threadIdx.x+l*blockDim.x];
    }
    __syncthreads();
  }
}



#define mydebug(str,...)
// if (threadIdx.x < 3 && blockIdx.x < 2){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
/// this computes the near field as a  ResistanceMatrix
/// pos           : is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// N             : is the number of particles
/// self_mobility : is 1./(6.*PI*eta*a)
/// a             : is the particle radius
/// L_g           : is the boxlength
/// resistance    : is the resistance matrix which will be retruned
/// myInfo        : contains infos about the operation:
///                myInfo[0] : number of overlapping particles
///                myInfo[1] : number of interacting particles (via nf)
///                myInfo[2] : max number of interacting particles per particle
__global__ void sd_compute_resistance_matrix(const real * pos, int N, real self_mobility, real a, const real * L_g, real * resistance, int * myInfo){
  int interactions=0;
  real mypos[3];
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  const int lda=(((N*3)+31)/32)*32;
  real myresistance[6]={0,0,0,0,0,0};
  int i = blockIdx.x*blockDim.x + threadIdx.x;
#ifndef SD_NOT_PERIODIC
  __shared__ real L[3];
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
#endif
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
  }

  __syncthreads();
  for (int d=0;d<3;d++){
    mypos[d] = cachedPos[threadIdx.x*3+d];
  }
  
  for (int offset=0;offset<N;offset+=numThreadsPerBlock){
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*numThreadsPerBlock+offset*3];
    }
    __syncthreads();
    for (int j=offset;j<min(offset+numThreadsPerBlock,N);j++){
      real dr[DIM];
      real dr2=0;
#pragma unroll
      for (int k=0;k<DIM;k++){
	dr[k]=mypos[k]-cachedPos[3*(j-offset)+k]; // r_ij
#ifndef SD_NOT_PERIODIC
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
#endif
	dr2+=dr[k]*dr[k];
	mydebug("dr[%d]: %f\n",k,dr[k]);
      }
#ifdef SD_RESISTANCE_CORRECT
      real r2bcorr_diag_self     = 0;
      real r2bcorr_diag_mix      = 0;
      real r2bcorr_offdiag_self  = 0;
      real r2bcorr_offdiag_mix   = 0;
#else
      real offdiag_fac=0;
      real diag_fac=0;
#endif
      if (i >= N || i ==j || j >= N){
	;
      }
      else if (dr2 < 4*a*4*a){
	if (!(2*a*2*a < dr2 )){
	  atomicAdd(myInfo,1); // count overlapping particles
	}
	else {// 2*a < drn < 4*a 
	  interactions++;
	  // python code:
	  // # Use only singular therms, namely to order O(s_ij^0)                                                                  
	  // T=(1./4./s-1/4-9./40.*ls)*dr*drt/dr2
	  // #           ^ this additonal constant is so that the mobility is smooth
	  // # c.f. N.-Q. Nguyen and A. J. C. Ladd, PHYSICAL REVIEW E 66, 046708 (2002) equation (34)                               
	  // T+=1./6.*ls*(-one+dr*drt/dr2)
	  // R[3*i:3*i+3,3*j:3*j+3]=-T
	  // R[3*i:3*i+3,3*i:3*i+3]+=T
	  real drn= sqrt(dr2); // length of dr
	  real s = drn/a-2;
	  real ls = log(s);
	  
	  real const para_fac_c=-0.125+(9./40.)*log(2.)+3./112.*2.*log(2.);
#ifdef SD_RESISTANCE_CORRECT
	  real para_fac    =(-0.25/s+(9./40.)*ls+(3./112.)*s*ls);
	  real perp_fac    =((1./6.)*ls);
#else
	  real const perp_fac_c=1./6.*log(2.);
	  real para_fac    =(-0.25/s+(9./40.)*ls+(3./112.)*s*ls-para_fac_c)/dr2/self_mobility;
	  real perp_fac    =((1./6.)*ls-perp_fac_c)/self_mobility;
	  diag_fac    = perp_fac;
	  offdiag_fac = para_fac-perp_fac;
#endif
#ifdef SD_RESISTANCE_CORRECT
	  real dr4=dr2*dr2;
	  real dr6=dr4*dr2;
	  // constants for correction at cutoff
	  const real dr_c1 = 4;
	  const real dr_c2 = 4*4;
	  const real dr_c3 = 4*4*4;
	  const real dr_c4 = 4*4*4*4;
	  const real dr_c5 = 4*4*4*4*4;
	  const real dr_c6 = 4*4*4*4*4*4;
	  const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	  const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	  const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)                 + 1./6.*log(2.);
	  const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1)      + 1./6.*log(2.);
	  // TODO: Make sure to use (real) and not double ...
	  // real computation
	  real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	  real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	  real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)      - r2bcorr_perp_self_c ));
	  real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)     - r2bcorr_perp_mix_c  ));
	  //printf("%d %d   show  [ %e, %e,  %e, %e ]\n",i,j,r2bcorr_para_self,r2bcorr_perp_self,r2bcorr_para_mix,r2bcorr_perp_mix);
	    
	  r2bcorr_diag_self     = (r2bcorr_perp_self)/self_mobility;
	  r2bcorr_diag_mix      = (r2bcorr_perp_mix )/self_mobility;
	  r2bcorr_offdiag_self  = (r2bcorr_para_self - r2bcorr_perp_self) /self_mobility/dr2;
	  r2bcorr_offdiag_mix   = (r2bcorr_para_mix  - r2bcorr_perp_mix ) /self_mobility/dr2;
#endif
	}
      }
      if (i < N){
#pragma unroll 3
	for (int k=0; k < DIM; k++){
#pragma unroll 3
	  for (int l=0;l < DIM; l++){
#ifdef SD_RESISTANCE_CORRECT
	    resistance[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*r2bcorr_offdiag_mix;
#else
	    resistance[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*offdiag_fac;
#endif
	  }
#ifdef SD_RESISTANCE_CORRECT
	  myresistance[k]-=dr[k]*dr[k]*r2bcorr_offdiag_self;
	  resistance[myindex(DIM*i+k,DIM*j+k)]+=r2bcorr_diag_mix;
	  myresistance[k]-=r2bcorr_diag_self;
#else
	  myresistance[k]-=dr[k]*dr[k]*offdiag_fac;
	  resistance[myindex(DIM*i+k,DIM*j+k)]+=diag_fac;
	  myresistance[k]-=diag_fac;
#endif
	}
      }
#ifdef SD_RESISTANCE_CORRECT
      myresistance[3]-=r2bcorr_offdiag_self*dr[0]*dr[1];
      myresistance[4]-=r2bcorr_offdiag_self*dr[0]*dr[2];
      myresistance[5]-=r2bcorr_offdiag_self*dr[1]*dr[2];
#else
      myresistance[3]-=offdiag_fac*dr[0]*dr[1];
      myresistance[4]-=offdiag_fac*dr[0]*dr[2];
      myresistance[5]-=offdiag_fac*dr[1]*dr[2];
#endif
      // python implementation:
      //T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
    }
  }
  if ( i < N){
#pragma unroll
    for (int k=0;k<3;k++){
      resistance[myindex(DIM*i+k,DIM*i+k)]=myresistance[k];
    }
    resistance[myindex(DIM*i+0,DIM*i+1)]=myresistance[3];
    resistance[myindex(DIM*i+1,DIM*i+0)]=myresistance[3];
    resistance[myindex(DIM*i+0,DIM*i+2)]=myresistance[4];
    resistance[myindex(DIM*i+2,DIM*i+0)]=myresistance[4];
    resistance[myindex(DIM*i+1,DIM*i+2)]=myresistance[5];
    resistance[myindex(DIM*i+2,DIM*i+1)]=myresistance[5];
  }
  __syncthreads();
  int * sharedInteractions = (int *) cachedPos; // reuse shared memory
  int * maxInteractions    = sharedInteractions + blockDim.x*2;
  sharedInteractions[threadIdx.x]=interactions;
  sharedInteractions[threadIdx.x+blockDim.x]=0;
  maxInteractions[threadIdx.x]   =interactions;
  maxInteractions[threadIdx.x+blockDim.x]   =0;
  for (int t=(blockDim.x+1)/2;t>1;t=(t+1)/2){
    if (threadIdx.x < t){
      sharedInteractions[threadIdx.x]+=sharedInteractions[threadIdx.x+t];
      sharedInteractions[threadIdx.x+t]=0;
      maxInteractions[threadIdx.x]=max(maxInteractions[threadIdx.x+t],maxInteractions[threadIdx.x]);
    }
    __syncthreads();
  }
  if (threadIdx.x==0){
    sharedInteractions[0]+=sharedInteractions[1];
    atomicAdd(myInfo+1, sharedInteractions[0]);
    maxInteractions[0]=max(maxInteractions[0],maxInteractions[1]);
    atomicMax(myInfo+2, maxInteractions[0]);
  }
}

/// Find particle pairs which are within a certain cutoff range
/// pos           : is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// L_g           : is the boxlength
/// N             : is the number of particles
/// interacting   : list of interaction particles. First entry is the first
///                 particle close enough to particle 1. Followed by the first
///                 one close enough to the second and so on.
///                 The second interaction particle of the first is at
///                 lda_short=((N+31)/32)*32 stored.
/// num_interacting : for each particle the number of interacting particles.
/// interaction_range : the cutoff range for interaction. Has to be smaller
///                 than half the boxlength.
__global__ void sd_find_interacting_particles(const real * pos,const real * L_g, const int N, int * interacting, int * num_interacting,
					      const real interaction_range){
  const int lda_short=((N+31)/32)*32;
  const real interaction_range_squared=SQR(interaction_range);
  __shared__ real cachedPos[3*numThreadsPerBlock];
#ifndef SD_NOT_PERIODIC
  __shared__ real L[3];
  if (threadIdx.x < 3){
    L[threadIdx.x]=L_g[threadIdx.x];
  }
#endif
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*blockDim.x] = pos[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
  }
  __syncthreads();
  real mypos[3];
  for (int l=0;l<3;l++){
    mypos[l]=cachedPos[threadIdx.x*3+l];
  }
  if (i < N){
    int interactions=0;
    int j0next;
    for (int j0=0;j0<N;j0+=j0next){
      j0next=min(j0+blockDim.x,N);
      j0next-=j0;
      for (int l=0;l<3;l++){
	cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*blockDim.x+j0*3];
      }
      for (int j=0;j<j0next;j++){
	real dr2=0;
	for (int l=0;l<3;l++){
	  real dr=mypos[l]-cachedPos[j*3+l];
#ifndef SD_NOT_PERIODIC
	  dr-=L[l]*rint(dr/L[l]);
#endif
	  dr2+=dr*dr;
	}
	if (dr2 < interaction_range_squared && i != j+j0) {
	  interacting[i+interactions*lda_short]=j+j0;
	  interactions++;
	}
      }
    }
    num_interacting[i]=interactions;
  }
}


/// Find particle pairs which are within a certain cutoff range using buckets
/// pos           : the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// _L            : the boxlength
/// N             : the number of particles
/// interacting   : list of interaction particles. First entry is the first
///                 particle close enough to particle 1. Followed by the first
///                 one close enough to the second and so on.
///                 The second interaction particle of the first is at
///                 lda_short=((N+31)/32)*32 stored.
/// num_interacting : for each particle the number of interacting particles.
/// particle_count : the number of particles within a certain bucket
/// particle_sorted_list : the list of the particles in the buckets
/// bucket_size_  : real[3] array of the dimensions of a bucket
/// bucket_num_   : int[3] array of the number of buckets in each direction
/// particle_to_bucket_list : lookup table what particle is in what bucket
/// interaction_range : the cutoff range for interaction. Has to be smaller
///                 than half the boxlength.
/// total_bucket_num : number of buckets
__global__ void sd_find_interacting_particles(const real * pos,const real * _L, const int N, int * interacting, int * num_interacting,
					      const int * particle_count, const int * particle_sorted_list, const real * bucket_size_,
					      const int * bucket_num_, const int * particle_to_bucket_list, const real interaction_range,
					      const int total_bucket_num){
  const int lda_short=((N+31)/32)*32;
  const real interaction_range_squared=SQR(interaction_range);
  __shared__ real bucket_size[3];
  __shared__ int bucket_num[6];
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  // copy small vectors to shared memory
#ifndef SD_NOT_PERIODIC
  __shared__ real L[3];
  if (threadIdx.x < 3){
    L[threadIdx.x]=_L[threadIdx.x];
  }
#endif
  if (threadIdx.x < 3){
    bucket_size[threadIdx.x]=bucket_size_[threadIdx.x];
    bucket_num[threadIdx.x]=bucket_num_[threadIdx.x];
  }
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
  }
  __syncthreads();
  int interactions=0;
  if ( i < N){
    int my_bucket=particle_to_bucket_list[i];
    int bucket_x=my_bucket%bucket_num[0];
    int bucket_y=(my_bucket/bucket_num[0])%bucket_num[1];
    int bucket_z=my_bucket/(bucket_num[0]*bucket_num[1]);
    int search_x=__real2int_ru(interaction_range/bucket_size[0]);
    int search_y=__real2int_ru(interaction_range/bucket_size[1]);
    int search_z=__real2int_ru(interaction_range/bucket_size[2]);
    int cz = bucket_z-search_z;
    while (cz < 0){
      cz+=bucket_num[2];
    }
    int cy0 = bucket_y-search_y;
    while (cy0 < 0){
      cy0 += bucket_num[1];
    }
    int cx0=bucket_x-search_x;
    while (cx0 < 0){
      cx0 += bucket_num[0];
    }
    int current_bucket  = cz*bucket_num[1];
    current_bucket += cy0;
    current_bucket *= bucket_num[0];
    current_bucket += cx0;
    for (int z = -search_z;z<=search_z;z++){
      int cy=cy0;
      int z_bucket=current_bucket;
      for (int y=-search_y;y<=search_y;y++){
	if (cy == bucket_num[1]){
	  cy = 0;
	  current_bucket -= bucket_num[0]*bucket_num[1];
	}
	int y_bucket = current_bucket;
	int cx=cx0;
	for (int x=-search_x;x<=search_x;x++){
	  if (cx == bucket_num[2]){
	    cx = 0;
	    current_bucket -= bucket_num[0];
	  }
	
	  //assert((current_bucket>=0 && current_bucket < total_bucket_num));
	  int in_loop=particle_count[current_bucket];
	  for (int j_loop=0;j_loop<in_loop;j_loop++){
	    int j=particle_sorted_list[current_bucket+j_loop*total_bucket_num];
	    real dr2=0;
#pragma unroll 3
	    for (int d=0;d<3;d++){
	      real dr=cachedPos[threadIdx.x*3+d]-pos[3*j+d];
#ifndef SD_NOT_PERIODIC
	      dr-=L[d]*rint(dr/L[d]);
#endif
	      //#warning folding back
	      dr2+=dr*dr;
	    }
	    if (dr2<interaction_range_squared){
	      interacting[i+interactions*lda_short]=j;
	      interactions++;
	    } 
	  }
	  cx++;
	  if (cx == bucket_num[0]){
	    cx = 0;
	    current_bucket -= bucket_num[0];
	  }
	  current_bucket++;
	} // end of x loop
	current_bucket = y_bucket;
	cy++;
	if (cy == bucket_num[1]){
	  cy = 0;
	  current_bucket -= bucket_num[0]*bucket_num[1];
	}
	current_bucket += bucket_num[0];
      }
      current_bucket  = z_bucket;
      cz++;
      if (cz == bucket_num[2]){
	cz = 0;
	current_bucket -= total_bucket_num;
      }
      current_bucket += bucket_num[0]*bucket_num[1];
    }
  }
  num_interacting[i]=interactions;
  // sort the interactions
  /*int loop = reduce_max((int *)cachedPos,interactions); // buggy
  if ( i < N){
    for (int j=0;j < loop;j++){
      if (j>=interactions){
	interacting[i+j*lda_short]=N+1;
      }
    }
    bool swapped;
    do {
      swapped=false;
      int last=interacting[i];
      int cur;
      for (int j=1;j<loop;j++){
	cur=interacting[i+lda_short*j];
	if (last > cur){
	  //swap
	  interacting[i+lda_short*(j-1)]=cur;
	  interacting[i+lda_short*j]    =last;
	  swapped=true;
	} else {
	  last=cur;
	}
      }
      loop--;
    } while (swapped);
  }*/
  if ( i < N){
    int loop = interactions;
    do {
      int last=interacting[i];
      int cur;
      int last_j=0;
      for (int j=1;j<loop;j++){
	cur=interacting[i+lda_short*j];
	if (last > cur){
	  //swap
	  interacting[i+lda_short*(j-1)]=cur;
	  interacting[i+lda_short*j]    =last;
	  last_j=j;
	} else {
	  last=cur;
	}
      }
      loop=last_j;
    } while (loop>0);
  }
  
}

/// Compute the lubrication correction in the resistance formulation.
/// pos           : is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
/// N             : is the number of particles
/// self_mobility : is 1./(6.*PI*eta*a)
/// a             : is the particle radius
/// _L            : is the boxlength
/// resistance    : is the sparse resistance matrix which will be returned
/// col_idx       : list of interacting partilces
/// row_l         : number of col_idx entries per particle
/// myInfo        : contains infos about the operation:
///                myInfo[0] : number of overlapping particles
///                myInfo[1] : number of interacting particles (via nf)
///                myInfo[2] : max number of interacting particles per particle
__global__ void sd_compute_resistance_matrix_sparse(const real * pos, const int N, const real self_mobility,const real a,const real * _L,
  						    real * resistance,const int * col_idx,const int * row_l, int * myInfo){
   const int lda_short=((N+31)/32)*32;
   const int lda = ((N*3+31)/32)*32;
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  // copy small vectors to shared memory
#ifndef SD_NOT_PERIODIC
  __shared__ real L[3];
  if (threadIdx.x < 3){
    L[threadIdx.x]=_L[threadIdx.x];
  }
#endif
  real myresistance[6]={0,0,0,0,0,0};
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  // get data for myposition - but coalscaled
#pragma unroll 3
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
  }
  
  __syncthreads();
  
  int num_interactions=i<N?row_l[i]:0;
  int i_index;
  for (int j_loop=0;j_loop<num_interactions;j_loop++){
    int j=col_idx[i+j_loop*lda_short];
    real dr[DIM];
    real dr2=0;
#pragma unroll
    for (int d=0;d<DIM;d++){
      dr[d]=cachedPos[threadIdx.x*3+d]-pos[3*j+d]; // r_ij
#ifndef SD_NOT_PERIODIC
      dr[d]-=L[d]*rint(dr[d]/L[d]); // fold back
#endif
      dr2+=dr[d]*dr[d];
      mydebug("dr[%d]: %f\n",d,dr[d]);
    }
#ifdef SD_RESISTANCE_CORRECT
    real r2bcorr_diag_self     = 0;
    real r2bcorr_diag_mix      = 0;
    real r2bcorr_offdiag_self  = 0;
    real r2bcorr_offdiag_mix   = 0;
#else
    real offdiag_fac=0;
    real diag_fac=0;
#endif
    if (i >= N || i ==j || j >= N){
      ;
    }
    else if (dr2 < 4*a*4*a){
      if (!(2*a*2*a < dr2 )){
	if (i==j){
	  i_index=j_loop;
	} else {
	  atomicAdd(myInfo,1); // count overlapping particles
	}
      }
      else {// 2*a < drn < 4*a 
	// python code:
	// # Use only singular therms, namely to order O(s_ij^0)                                                                  
	// T=(1./4./s-1/4-9./40.*ls)*dr*drt/dr2
	// #           ^ this additonal constant is so that the mobility is smooth
	// # c.f. N.-Q. Nguyen and A. J. C. Ladd, PHYSICAL REVIEW E 66, 046708 (2002) equation (34)                               
	// T+=1./6.*ls*(-one+dr*drt/dr2)
	// R[3*i:3*i+3,3*j:3*j+3]=-T
	// R[3*i:3*i+3,3*i:3*i+3]+=T
	real drn= sqrt(dr2); // length of dr
	real s = drn/a-2;
	real ls = log(s);
	
	real const para_fac_c=-0.125+(9./40.)*log(2.)+3./112.*2.*log(2.);
#ifdef SD_RESISTANCE_CORRECT
	real para_fac    =(-0.25/s+(9./40.)*ls+(3./112.)*s*ls);
	real perp_fac    =((1./6.)*ls);
#else
	real const perp_fac_c=1./6.*log(2.);
	real para_fac    =(-0.25/s+(9./40.)*ls+(3./112.)*s*ls-para_fac_c)/dr2/self_mobility;
	real perp_fac    =((1./6.)*ls-perp_fac_c)/self_mobility;
	diag_fac    = perp_fac;
	offdiag_fac = para_fac-perp_fac;
#endif
#ifdef SD_RESISTANCE_CORRECT
	real dr4=dr2*dr2;
	real dr6=dr4*dr2;
	// constants for correction at cutoff
	const real dr_c1 = 4;
	const real dr_c2 = 4*4;
	const real dr_c3 = 4*4*4;
	const real dr_c4 = 4*4*4*4;
	const real dr_c5 = 4*4*4*4*4;
	const real dr_c6 = 4*4*4*4*4*4;
	
	const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)                 + 1./6.*log(2.);
	const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1)      + 1./6.*log(2.);
	// TODO: Make sure to use (real) and not double ...
	// real computation
	real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)      - r2bcorr_perp_self_c ));
	real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)     - r2bcorr_perp_mix_c  ));
	//printf("%d %d   show  [ %e, %e,  %e, %e ]\n",i,j,r2bcorr_para_self,r2bcorr_perp_self,r2bcorr_para_mix,r2bcorr_perp_mix);
	
	r2bcorr_diag_self     = (r2bcorr_perp_self)/self_mobility;
	r2bcorr_diag_mix      = (r2bcorr_perp_mix )/self_mobility;
	r2bcorr_offdiag_self  = (r2bcorr_para_self - r2bcorr_perp_self) /self_mobility/dr2;
	r2bcorr_offdiag_mix   = (r2bcorr_para_mix  - r2bcorr_perp_mix ) /self_mobility/dr2;
#endif
      }
    }
    if (i < N){
#pragma unroll 3
      for (int k=0; k < DIM; k++){
#pragma unroll 3
	for (int l=0;l < DIM; l++){
#ifdef SD_RESISTANCE_CORRECT
	  resistance[DIM*i+k+(DIM*j_loop+l)*lda]=dr[k]*dr[l]*r2bcorr_offdiag_mix;
#else
	  resistance[DIM*i+k+(DIM*j_loop+l)*lda]=dr[k]*dr[l]*offdiag_fac;
#endif
	}
#ifdef SD_RESISTANCE_CORRECT
	myresistance[k]-=dr[k]*dr[k]*r2bcorr_offdiag_self;
	resistance[DIM*i+k+(DIM*j_loop+k)*lda]+=r2bcorr_diag_mix;
	myresistance[k]-=r2bcorr_diag_self;
#else
	myresistance[k]-=dr[k]*dr[k]*offdiag_fac;
	resistance[DIM*i+k+(DIM*j_loop+k)*lda]+=diag_fac;
	myresistance[k]-=diag_fac;
#endif
      }
    }  
#ifdef SD_RESISTANCE_CORRECT
    myresistance[3]-=r2bcorr_offdiag_self*dr[0]*dr[1];
    myresistance[4]-=r2bcorr_offdiag_self*dr[0]*dr[2];
    myresistance[5]-=r2bcorr_offdiag_self*dr[1]*dr[2];
#else
    myresistance[3]-=offdiag_fac*dr[0]*dr[1];
    myresistance[4]-=offdiag_fac*dr[0]*dr[2];
    myresistance[5]-=offdiag_fac*dr[1]*dr[2];
#endif
    // python implementation:
    //T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
  }
  
  /*else{ // set the block to zero
  // it might be faster to set everything in the beginning to zero ...
  // or use sparse matrices ...
  #pragma unroll 3
  for (int k=0; k < DIM; k++){
  #pragma unroll 3
  for (int l=0;l < DIM; l++){
  resistance[myindex(DIM*i+k,DIM*j+l)]=0;
  }
  }  
  }*/
  
  if ( i < N){
    resistance[DIM*i+0+(DIM*i_index+0)*lda]=myresistance[0];
    resistance[DIM*i+0+(DIM*i_index+1)*lda]=myresistance[3];
    resistance[DIM*i+0+(DIM*i_index+2)*lda]=myresistance[4];
    resistance[DIM*i+1+(DIM*i_index+0)*lda]=myresistance[3];
    resistance[DIM*i+1+(DIM*i_index+1)*lda]=myresistance[1];
    resistance[DIM*i+1+(DIM*i_index+2)*lda]=myresistance[5];
    resistance[DIM*i+2+(DIM*i_index+0)*lda]=myresistance[4];
    resistance[DIM*i+2+(DIM*i_index+1)*lda]=myresistance[5];
    resistance[DIM*i+2+(DIM*i_index+2)*lda]=myresistance[2];
  }
  /*
  __syncthreads();
  int * sharedInteractions = (int *) cachedPos; // reuse shared memory
  int * maxInteractions    = sharedInteractions + blockDim.x*2;
  sharedInteractions[threadIdx.x]=num_iterations-1;
  sharedInteractions[threadIdx.x+blockDim.x]=0;
  maxInteractions[threadIdx.x]   =num_iterations-1;
  maxInteractions[threadIdx.x+blockDim.x]   =0;
  for (int t=(blockDim.x+1)/2;t>1;t=(t+1)/2){
    if (threadIdx.x < t){
      sharedInteractions[threadIdx.x]+=sharedInteractions[threadIdx.x+t];
      sharedInteractions[threadIdx.x+t]=0;
      maxInteractions[threadIdx.x]=max(maxInteractions[threadIdx.x+t],maxInteractions[threadIdx.x]);
    }
    __syncthreads();
  }
  if (threadIdx.x==0){
    sharedInteractions[0]+=sharedInteractions[1];
    atomicAdd(myInfo+1, sharedInteractions[0]);
    maxInteractions[0]=max(maxInteractions[0],maxInteractions[1]);
    atomicMax(myInfo+2, maxInteractions[0]);
    }*/
}


#ifndef SD_RESISTANCE_CORRECT
#warning "SD Brownian motion only support corrected resistance calculation ..."
#endif
/// Computing the brownian forces which arise from the nearfield contribution
/// of the mobility matrix
/// r        : position of the particles
/// gaussian : pointer to a sufficent amount of gaussian distributed random
///            numbers 
/// N        : number of particles
/// L_g      : size of periodic box
/// a        : hydrodynamic particle radius
/// self_mobility : is 1./(6.*PI*eta*a)
/// brownian_force_nf : computed brownian force
__global__ void sd_compute_brownian_force_nearfield(const real * r, const real * gaussian,int N,const real * L_g,
						    real a, real self_mobility,real * brownian_force_nf){
  const int gaussian_ldd=((N+31)/32)*32;
  int interactions=0;
  real mypos[3];
  real writeCache[6];
  //real otherWriteCache[3];
  __shared__ real cachedPos[3*numThreadsPerBlock];
  __shared__ real choleskyCache[12*numThreadsPerBlock];
  //const int lda=(((N*3)+31)/32)*32;
  //__shared__ real myresistance[6*numThreadsPerBlock];
  //real myresistance[6];
  //__shared__ real otherresistance[6*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
#ifndef SD_NOT_PERIODIC
  __shared__ real L[3];
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
#endif
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = r[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
    writeCache[l]= 0;
  }
  __syncthreads();
  for (int d=0;d<3;d++){
    mypos[d] = cachedPos[threadIdx.x*3+d];
  }
  
  for (int offset=0;offset<N;offset+=numThreadsPerBlock){
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      cachedPos[threadIdx.x+l*numThreadsPerBlock] = r[threadIdx.x+l*numThreadsPerBlock+offset*3];
    }
    __syncthreads();
for (int j=max(offset, blockIdx.x*numThreadsPerBlock+1);j<min(offset+numThreadsPerBlock,N);j++){
      real dr[DIM];
      real dr2=0;
#pragma unroll
      for (int k=0;k<DIM;k++){
	dr[k]=mypos[k]-cachedPos[3*(j-offset)+k]; // r_ij
#ifndef SD_NOT_PERIODIC
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
#endif
	dr[k]/=a;
	dr2+=dr[k]*dr[k];
      }
      real r2bcorr_diag_self     = 0;
      real r2bcorr_diag_mix      = 0;
      real r2bcorr_offdiag_self  = 0;
      real r2bcorr_offdiag_mix   = 0;
      
#ifdef SD_DEBUG
      for (int t=0;t<12;t++){
	choleskyCache[threadIdx.x+ (t)*numThreadsPerBlock]=0;
      }
#endif
      int wasInLoop = 0;
      if (i >= N || i >= j || j >= N){
	writeCache[3]=0;
	writeCache[4]=0;
	writeCache[5]=0;
      }
      // j > i
      else if (dr2 < 16  && 4 < dr2 ){// 2*a < drn < 4*a 
	wasInLoop = 1;
	{
	  real drn= sqrt(dr2); // length of dr
	  real s = drn-2;
	  real ls = log(s);
	  
	  real const para_fac_c=-0.125+(9./40.)*log(2.)+3./112.*2.*log(2.);
	  real para_fac    =(-0.25/s+(9./40.)*ls+(3./112.)*s*ls);
	  real perp_fac    =((1./6.)*ls);
	  
	  real dr4=dr2*dr2;
	  real dr6=dr4*dr2;
	  // constants for correction at cutoff
	  const real dr_c1 = 4;
	  const real dr_c2 = 4*4;
	  const real dr_c3 = 4*4*4;
	  const real dr_c4 = 4*4*4*4;
	  const real dr_c5 = 4*4*4*4*4;
	  const real dr_c6 = 4*4*4*4*4*4;
	  
	  const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	  const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	  const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)                 + 1./6.*log(2.);
	  const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1)      + 1./6.*log(2.);
	  // TODO: Make sure to use (real) and not double ...
	  // real computation
	  real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	  real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	  real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)      - r2bcorr_perp_self_c ));
	  real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)     - r2bcorr_perp_mix_c  ));
	  //printf("%d %d   show  [ %e, %e,  %e, %e ]\n",i,j,r2bcorr_para_self,r2bcorr_perp_self,r2bcorr_para_mix,r2bcorr_perp_mix);
	  
	  r2bcorr_diag_self     = (r2bcorr_perp_self)/self_mobility;
	  r2bcorr_diag_mix      = (r2bcorr_perp_mix )/self_mobility;
	  r2bcorr_offdiag_self  = (r2bcorr_para_self - r2bcorr_perp_self) /self_mobility/dr2;
	  r2bcorr_offdiag_mix   = (r2bcorr_para_mix  - r2bcorr_perp_mix ) /self_mobility/dr2;
	  //printf("%d %d   shoz     [%e, %e, %e, %e]\n",i,j,r2bcorr_diag_self,r2bcorr_offdiag_self,r2bcorr_diag_mix,r2bcorr_offdiag_mix);
	  //r2bcorr_offdiag_self /=  dr2;
	  //r2bcorr_offdiag_mix  /=  dr2;
	}
	// This is the cholesky decomposition.
	// note that we try to avoid the usage of registers, so we use shared mem
	// myCC is a makro, defined here to shorten the lines:
#define myCC(pos) choleskyCache[threadIdx.x+ (pos)*numThreadsPerBlock]
	// without it would look more like this:
	//choleskyCache[threadIdx.x+ 0*numThreadsPerBlock] = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[0]*dr[0]);
	//choleskyCache[threadIdx.x+ 1*numThreadsPerBlock] = r2bcorr_offdiag_self*dr[0]*dr[1] / choleskyCache[threadIdx.x+ 0*numThreadsPerBlock];
	// L_{1,1} to L_{6,1}
	myCC(0)  = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[0]*dr[0]);
	myCC(1)  =                        r2bcorr_offdiag_self*dr[0]*dr[1] / myCC(0);
	myCC(2)  =                        r2bcorr_offdiag_self*dr[0]*dr[2] / myCC(0);
	myCC(3)  =    (r2bcorr_diag_mix + r2bcorr_offdiag_mix *dr[0]*dr[0])/ myCC(0);
	myCC(4)  =                        r2bcorr_offdiag_mix *dr[0]*dr[1] / myCC(0);
	myCC(5)  =                        r2bcorr_offdiag_mix *dr[0]*dr[2] / myCC(0);
	
	writeCache[0]+=myCC(0)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[1]+=myCC(1)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[2]+=myCC(2)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[3] =myCC(3)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[4] =myCC(4)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5] =myCC(5)  * gaussian[0*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 6
	// L_{2,2} to L_{6,2}
	myCC(0)  = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[1]*dr[1] - SQR(    myCC(1)));
	myCC(6)  =                       (r2bcorr_offdiag_self*dr[1]*dr[2] - myCC(2)*myCC(1))/myCC(0);
	myCC(7)  =                       (r2bcorr_offdiag_mix *dr[1]*dr[0] - myCC(3)*myCC(1))/myCC(0);
	myCC(8)  =     (r2bcorr_diag_mix +r2bcorr_offdiag_mix *dr[1]*dr[1] - myCC(4)*myCC(1))/myCC(0);
	myCC(9)  =                       (r2bcorr_offdiag_mix *dr[1]*dr[2] - myCC(5)*myCC(1))/myCC(0);
	writeCache[1]+=myCC(0)  * gaussian[1*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[2]+=myCC(6)  * gaussian[1*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[3]+=myCC(7)  * gaussian[1*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[4]+=myCC(8)  * gaussian[1*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(9)  * gaussian[1*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 11 - 1
	// L_{3,3} to L_{6,3}
	myCC(0)  = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[2]*dr[2] - SQR(    myCC(2))- SQR(    myCC(6)));
	myCC(1)  =                       (r2bcorr_offdiag_mix *dr[2]*dr[0] - myCC(3)*myCC(2) - myCC(7)*myCC(6))/myCC(0);
	myCC(10) =                       (r2bcorr_offdiag_mix *dr[2]*dr[1] - myCC(4)*myCC(2) - myCC(8)*myCC(6))/myCC(0);
	myCC(11) =     (r2bcorr_diag_mix +r2bcorr_offdiag_mix *dr[2]*dr[2] - myCC(5)*myCC(2) - myCC(9)*myCC(6))/myCC(0);
	writeCache[2]+=myCC(0)  * gaussian[2*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[3]+=myCC(1)  * gaussian[2*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[4]+=myCC(10) * gaussian[2*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(11) * gaussian[2*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 15 - 3
	// L_{4,4} to L_{6,4}
	myCC(0)  = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[0]*dr[0] - SQR(    myCC(3))- SQR(    myCC(7))
			- SQR(     myCC(1)));
	if (isnan(myCC(0))){
	  printf("%d %d : [%e %e %e]\n",i,j,dr[0],dr[1],dr[2]);
	}
	myCC(2)  =                       (r2bcorr_offdiag_self*dr[0]*dr[1] - myCC(4)*myCC(3) - myCC(8)*myCC(7) 
					  - myCC(10)*myCC(1))/myCC(0);
	myCC(6)  =                       (r2bcorr_offdiag_self*dr[0]*dr[2] - myCC(5)*myCC(3) - myCC(9)*myCC(7) 
					  - myCC(11)*myCC(1))/myCC(0);
	writeCache[3]+=myCC(0)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[4]+=myCC(2)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(6)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 18 - 6
	// L_{5,5} and L_{6,5}
	myCC(0)  = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[1]*dr[1] - SQR(    myCC(4))- SQR(    myCC(8))
			- SQR(     myCC(10))- SQR(    myCC(2)));
	if (isnan(myCC(0))){
	  printf("%d %d : [%e %e %e] \n",i,j,dr[0],dr[1],dr[2]);
	}
	myCC(3)  =                       (r2bcorr_offdiag_self*dr[1]*dr[2] - myCC(5)*myCC(4) - myCC(9)*myCC(8) 
					  - myCC(11)*myCC(10) - myCC(6)*myCC(2))/myCC(0);
	writeCache[4]+=myCC(0)  * gaussian[4*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(3)  * gaussian[4*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 20 - 10
	// L_{6,6} would be:
	myCC(0) = sqrt(r2bcorr_diag_self+r2bcorr_offdiag_self*dr[2]*dr[2] - SQR(myCC(5))    - SQR(myCC(9))     
		       - SQR(myCC(11)) - SQR(myCC(6)) - SQR(myCC(3)));
	if (isnan(myCC(0))){
	  printf("%d %d : [%e %e %e] \n",i,j,dr[0],dr[1],dr[2]);
	}
	writeCache[5]+=myCC(0)  * gaussian[5*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used 21 - 15
	interactions++;
	for (int l=0;l<6;l++){
	  if (isnan(writeCache[l])){
	    printf("%d %d, %d: dr=[%e, %e, %e]\n",i,j,l,dr[0],dr[1],dr[2]);
	  }
	}
      }
      // for the particle j (writeCache[3-5]) we can reduce localy:
      
      int * haveInteraction = (int *) choleskyCache+6*numThreadsPerBlock; // reuse shared memory
      choleskyCache[threadIdx.x+0*numThreadsPerBlock]=writeCache[3];
      choleskyCache[threadIdx.x+1*numThreadsPerBlock]=0;
      choleskyCache[threadIdx.x+2*numThreadsPerBlock]=writeCache[4];
      choleskyCache[threadIdx.x+3*numThreadsPerBlock]=0;
      choleskyCache[threadIdx.x+4*numThreadsPerBlock]=writeCache[5];
      choleskyCache[threadIdx.x+5*numThreadsPerBlock]=0;
      haveInteraction[threadIdx.x]=wasInLoop;
      haveInteraction[threadIdx.x+numThreadsPerBlock]=0;
      for (int t=(blockDim.x+1)/2;t>1;t=(t+1)/2){
	if (threadIdx.x < t){
	  choleskyCache[threadIdx.x]+=choleskyCache[threadIdx.x+t];
	  choleskyCache[threadIdx.x+2*numThreadsPerBlock]+=choleskyCache[threadIdx.x+t +2*numThreadsPerBlock];
	  choleskyCache[threadIdx.x+4*numThreadsPerBlock]+=choleskyCache[threadIdx.x+t +2*numThreadsPerBlock];
	  haveInteraction[threadIdx.x]|=haveInteraction[threadIdx.x+t];
	  choleskyCache[threadIdx.x+t]=0;
	  choleskyCache[threadIdx.x+t +2*numThreadsPerBlock]=0;
	  choleskyCache[threadIdx.x+t +4*numThreadsPerBlock]=0;
	  haveInteraction[threadIdx.x+t]=0;
	}
	__syncthreads();
      }
      if (threadIdx.x==0){
	if (haveInteraction[0] || haveInteraction[1]){
	  choleskyCache[0]+=choleskyCache[1];
	  choleskyCache[2*numThreadsPerBlock]+=choleskyCache[1+2*numThreadsPerBlock];
	  choleskyCache[4*numThreadsPerBlock]+=choleskyCache[1+4*numThreadsPerBlock];
	  real tmp=atomicAdd(brownian_force_nf+j*3,   choleskyCache[0]);
	  bool error=isnan(tmp);
	  tmp = atomicAdd(brownian_force_nf+j*3+1, choleskyCache[2*numThreadsPerBlock]);
	  error|=isnan(tmp);
	  tmp = atomicAdd(brownian_force_nf+j*3+2, choleskyCache[4*numThreadsPerBlock]);
	  error|=isnan(tmp);
	  if (error){
	    printf("%d %d: dr=[%e, %e, %e]\n",i,j,dr[0],dr[1],dr[2]);
	  }
	}
      }
    }
  }
  if ( i < N){
#pragma unroll 3
    for (int k=0;k<3;k++){
      real tmp = atomicAdd(brownian_force_nf+i*3+k, writeCache[k]);
      //assert(!isnan(tmp));
    }
  }
}

/// computing the matrix vector y = c * A * x product with a square sparse
/// matrix A and a scalar factor c
/// size     : size of the matrix and of the input and output vectors
/// factor   : factor $c$ is the given prefactor
/// matrix   : data of the sparse matrix
/// ldd      : leading dimension of the matrix
/// ldd_short : leading dimension of col_idx
/// col_idx  : indices of the entries in the matrix
/// row_l    : number of entris per row
/// vec_in   : the input vecor $x$
/// vec_out  : the output vector $y$
__global__ void sd_multiply_sparse_matrix_vector(const int size, const real factor, const real * matrix, const int ldd, const int ldd_short,
						 const int * col_idx, const int * row_l, const real * vec_in, real * vec_out ){
  int id = threadIdx.x + blockIdx.x*blockDim.x;
  if (id<size){
    int i=id/3;
    int rows=row_l[i];
    real erg=0;
    for (int row=0;row<rows;row++){
      int colid=col_idx[i+row*ldd_short];
      for (int k=0;k<3;k++){
	// asm version is about 33% faster
	erg+=read_without_caching(matrix+id+(row*3+k)*ldd)*vec_in[colid*3+k];
	//erg+=matrix[id+(row*3+k)*ldd]*vec_in[colid*3+k];
      }
    }
    vec_out[id]=erg*factor;
  }
}

/// this adds the identity matrix to a given matrix of ld=size
/// matrix   : pointer to the given matrix
/// size     : the size of the matrix (in the example below 3N)
/// block    : (ignored) the number of elements to process per thread
///            if this is e.g. 3 and the matrix is 3Nx3N, than N threads have
///            to be started 
__global__ void sd_add_identity_matrix(real * matrix, int size, int lda){
  //int lda=((size+31)/32)*32;
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //for (int i = idx*block; i< (idx+1)*block; i++){
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    matrix[i+i*lda]+=1;
  }
}

/// this sets a block to zero
/// matrix   : pointer to the given matrix
/// size     : the size of the matrix (in the example below 3N)
/// ldd      : the leading dimension of the matrix
__global__ void sd_set_zero_matrix(real * matrix, int size, int ldd){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  int matsize=ldd;
  matsize*=size;
  for (int i = idx;i< matsize; i+=blockDim.x*gridDim.x){
    matrix[i]=0;
  }
}


/// this sets a block to zero
/// data     : pointer to the given data
/// size     : the size of the data
__global__ void sd_set_zero(real * data, int size){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    data[i]=0;
  }
}

/// this sets a block to an given integer
/// data     : pointer to the given data
/// size     : the size of the data
/// value    : the value written to the data block
__global__ void sd_set_value(int * data, int size, int value){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    data[i]=value;
  }
}

/// this sets a block to an given real
/// data     : pointer to the given data
/// size     : the size of the data
/// value    : the value written to the data block
__global__ void sd_set_value(real * data, int size, real value){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    data[i]=value;
  }
}




#define DIST (2+1e-1)
#define DISP_MAX (1000)
/// Checking before the actual integration whether any particles moves more
/// than DISP_MAX and if so, reducing the displacement to DISP_MAX
/// r_d      : position of the particles
/// disp_d   : displacements of the particles
/// L        : box dimension
/// a        : particle radius
/// N        : number of particles
__global__ void sd_real_integrate_prepare( const real * r_d , real * disp_d, const real * L, real a, int N){
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  i*=3;
  real disp2;
#pragma unroll
  for (int d=0;d<3;d++){
    disp2+=disp_d[i+d]*disp_d[i+d];
  }
  if (disp2> DISP_MAX*DISP_MAX){
    disp2=DISP_MAX/sqrt(disp2);
#pragma unroll
    for (int d=0;d<3;d++){
      disp_d[i+d]*=disp2;
    }
  }
}

/// if SD_AVOID_COLLISION is set avoid that particles are overlapping after
/// the displacement, else only perform the integration step
/// r_d      : position of the particles
/// disp_d   : displacement of the particles
/// L        : box dimensions
/// a        : particle radius
/// N        : number of particles
__global__ void sd_real_integrate( real * r_d , const real * disp_d, const real * L, real a, int N)
{
  int idx =  blockIdx.x*blockDim.x + threadIdx.x;
  // t is the factor how far of disp_d we will move.
  // in case everything is fine, we will move t, if there is some trouble,
  // we will move less to avoid collision
  real t=1;
#ifdef SD_AVOID_COLLISON
  real rnew[DIM];
  for (int d=0;d<DIM;d++){
    rnew[d]=r_d[DIM*idx+d]+disp_d[DIM*idx+d];
  }
  const real distmin=(3*a)*(3*a);
  for (int i=0;i<N;i++){
    if (idx==i){
      i++;
      if (i >N){
	continue;
      }
    }
    real dr2=0;
    for (int d=0;d<DIM;d++){
      real tmp=r_d[i*DIM+d]-rnew[d];
#ifndef SD_NOT_PERIODIC
      tmp-=L[d]*rint(tmp/L[d]);
#endif
      dr2+=tmp*tmp;
    }
    if (dr2 <distmin){ // possible colision - check better
      dr2=0;
      //real dr2o=0; // or do we need old distance?
      for (int d=0;d<DIM;d++){
	real tmp=r_d[i*DIM+d]+disp_d[i*DIM+d]-rnew[d];
#ifndef SD_NOT_PERIODIC
	tmp-=L[d]*rint(tmp/L[d]);
#endif
	dr2+=tmp*tmp;
	//tmp=r_d[i*DIM+d]-r_d[idx*DIM+d];
	//tmp-=L*rint(tmp/L);
	//dr2o+=tmp*tmp;
      }
      if (dr2 < DIST*DIST*a*a){ // do they collide after the step?
	// ideal: the motion which is responsible for the crash: avoid it.
	// just move them that far that they nearly touch each other.
	// therefore we need the soluten of an quadratic equation
	// in case they are already closer than DIST*a this will move them appart.
	// first: get the coefficents
	real alpha=0,beta=0,gamma=0;
	for (int d=0;d<DIM;d++){
	  real t1=r_d[i*DIM+d]-r_d[idx*DIM+d];
#ifndef SD_NOT_PERIODIC
	  t1-=L[d]*rint(t1/L[d]);
#endif
	  real t2=disp_d[i*DIM+d]-disp_d[idx*DIM+d];
	  //t2-=L*rint(t2/L); // we would have a problem if we would need to fold back these ...
	  alpha +=t2*t2;
	  beta  +=2*t1*t2;
	  gamma +=t1*t1;
	} 
	// now we want to solve for t: alpha*t**2+beta*t+gamma=DIST*a
	// we want the solution with the minus in the 'mitternachtsformel'
	// because the other solution is when the particles moved through each other
	real tnew = (-beta-sqrt(beta*beta-4*alpha*gamma))/(2*alpha);
	if (tnew < t){ // use the smallest t
	  t=tnew;
	}
      }
    }
  }
#endif
  for (int d=0;d<DIM;d++){ // actually do the integration
    r_d[DIM*idx+d]+=disp_d[DIM*idx+d]*t;
  }
  //#warning "Debug is still enabaled"
    //pos_d[DIM*N+idx]=t;
}

/// sorting the particles in buckets
/// pos      : position of the particles
/// bucketSize : dimension of bucket
/// bucketNum : number of buckets in each direction
/// N        : number of particles
/// particleCount : number of particles within a bucket
/// particleList : list of the particles within a given bucket
/// totalBucketNum : number of buckets, product of the three values of
///            bucketNum
/// particleToBucketList : array containg the bucket number for each particle
__global__ void sd_bucket_sort( const real * pos , const real * bucketSize, const int * bucketNum, int N, int * particleCount,
				int * particleList, int maxParticlePerCell, int totalBucketNum, int * particleToBucketList){
  for (int i = blockIdx.x*blockDim.x + threadIdx.x;
       i<N ;
       i+=blockDim.x*gridDim.x){
    int3 bucket;
#pragma unroll 3
    for (int d =0; d<3; d++){
      real tmp;
      // no asm version:
      // tmp = pos[i*3+d];
      // asm version avoids caching
#ifdef SD_USE_FLOAT
      asm("ld.global.cs.f32 %0,[%1];\n"
	: "=f"(tmp) : "l"(pos+i*3+d) : );
#else
      asm("ld.global.cs.f64 %0,[%1];\n"
	: "=d"(tmp) : "l"(pos+i*3+d) : );
#endif
      tmp/=bucketSize[d];
      int x;
      // this should work with compute capability >= 2.0
      x=__real2int_rd(tmp);
      //x=floor(tmp);
      x%=bucketNum[d];
      // avoid negativ numbers
      x= (x < 0)?x+bucketNum[d]: x;
      //x+=bucketNum[d];
      //x%=bucketNum[d];
      switch (d){
      case 0:
	bucket.x = x;
	break;
      case 1:
	bucket.y = x;
	break;
      case 2:
	bucket.z = x;
	break;
      }
    }
    int myBucket = bucket.x + bucket.y*bucketNum[0] + bucket.z*bucketNum[0]*bucketNum[1];
    int num = atomicAdd(particleCount+myBucket, 1);
    if (num < maxParticlePerCell){
      particleList[myBucket+num*totalBucketNum]=i;
      particleToBucketList[i]=myBucket;
    }else{
      // Note: printf in device code works only with c.c.>=2.0 
#if (__CUDA_ARCH__>=200)
      printf("error: overflow in grid cell (%i,%i,%i)\n",bucket.x,bucket.y,bucket.z);
#endif
    }
  }
}

/// This function reads the diagonal element entries of a square matrix
/// and stores them in the vector diag and the inverse of them in diag_i
/// size     : size of the matrix
/// mat_a    : the matrix
/// lda      : leading dimension of the matrix mat_a
/// diag     : the diagonal (OUT)
/// diag_i   : the inverse of the diagonal (OUT)
__global__ void sd_get_diag(int size, const real * mat_a,int lda,real * diag,real * diag_i){
  const int stepw=lda+1;
  for (int l=threadIdx.x+blockDim.x*blockIdx.x;l<size;l+=blockDim.x*gridDim.x){
    real tmp;
    tmp = mat_a[stepw*l];
    diag[l]=tmp;
    diag_i[l]=1/tmp;
  }
}

/// computing the 1 norm of a given vector
/// size     : size of the vector
/// vec      : data of the vector
/// res      : the 1 norm (OUT)
__global__ void sd_nrm1(const int size, const real * const vec, real * res){
  __shared__ real cache[numThreadsPerBlock];
  if (blockIdx.x == 0){ // only use one block
    //assert(blockDim.x == numThreadsPerBlock);
    cache[threadIdx.x]=0;
    for (int i = threadIdx.x; i < size ; i+=blockDim.x){
      cache[threadIdx.x]+=abs(vec[i]);
    }
    reduce_sum(cache);
    if (threadIdx.x == 0 ){
      res[0]=cache[0];
    }
  }
}

/// computing the 1 norm of a given vector
/// size     : size of the vector
/// vec      : data of the vector
/// res      : the 1 norm (OUT)
__global__ void sd_nrm_inf(const int size, const int * const vec, int * res){
  extern __shared__ int cache2[]; // blockDim.x * sizeof(int)
  //if (blockIdx.x == 0){ // only use one block
  cache2[threadIdx.x]=vec[threadIdx.x];
  for (int i = threadIdx.x+blockDim.x; i < size ; i+=blockDim.x){
    if ( cache2[threadIdx.x] < vec[i])
      cache2[threadIdx.x]=vec[i];
  }
  __syncthreads();
  reduce_max(cache2);
  if (threadIdx.x == 0 ){
    res[blockIdx.x]=cache2[0];
  }
  //}
}

#endif /* SD */
