#define CUDA

#include <stdio.h>
#include <assert.h>
#include "cuda_runtime.h"
#include <device_functions.h>

#include "integrate_sd_cuda.cuh"
#include "integrate_sd_cuda_device.cu"

/* *************************************************************************************************************** *
 * ********************************************      CUDA-KERNELS     ******************************************** *
 * *************************************************************************************************************** */


// This computes the farfield contribution of the mobility
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L_d is the boxlength
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
__global__ void sd_compute_mobility_matrix(real * r, int N, real self_mobility, real a, real * L_g, real * mobility){
  real mypos[3];
  const int lda=((3*N+31)/32)*32;
  __shared__ real L[3];
  __shared__ real cachedPos[3*numThreadsPerBlock];
  __shared__ real writeCache[3*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (threadIdx.x < 3){ // copy L to shared memory
    //mydebug("0x%08x  \n",L_g + threadIdx.x);
    L[threadIdx.x]=L_g[threadIdx.x];
  }
  __syncthreads();
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
#warning "Disabled fold back to avoid negative eigenvalues!"
	  //dr[k]-=rint(dr[k]/L[k])*L[k]; // fold back
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
	// python implementation:
	//T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
	//}
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
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L_d is the boxlength
// cutoff the realspace cutoff distance
// xi the splitting parameter as defined by Beenakker 1986
// xa  = xi * a
// xa3 = xa * xa * xa
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
 __global__ void sd_compute_mobility_matrix_real_short(real * r, int N, real self_mobility, real a, real * L_g, real * mobility, real cutoff, real xi, real xa, real xa3){
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
	if (0.5 < ar || drn < cutoff){  // drn < 2*a
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
	  t=(0.75f-1.5f*ar2)*ar*erfcf(xr)+(-4.f*xa3*xr2*xr2-3.f*xa*xr2+16f*xa3*xr2+1.5f*xa-2.f *xa3-3.f*xa*ar2)*0.5641895835477562869480794515607725858440506f*exp(-xr2);
	  t2=(0.75f+0.5f*ar2)*ar*erfcf(xr)+(4.f*xa3*xr2*xr2+3.f*xa*xr2-20f*xa3*xr2-4.5f*xa+14.f*xa3+    xa*ar2)*0.5641895835477562869480794515607725858440506f*exp(-xr2);
#else
	  t=(0.75-1.5*ar2)*ar*erfc(xr)+(-4.*xa3*xr2*xr2-3.*xa*xr2+16.*xa3*xr2+1.5*xa-2. *xa3-3.*xa*ar2)*0.5641895835477562869480794515607725858440506*exp(-xr2);
	  t2=(0.75+0.5*ar2)*ar*erfc(xr)+(4.*xa3*xr2*xr2+3.*xa*xr2-20.*xa3*xr2-4.5*xa+14.*xa3+   xa*ar2)*0.5641895835477562869480794515607725858440506*exp(-xr2);
#endif
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
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned (in/out)
// L_d is the boxlength
// cutoff the wavespace cutoff distance
// xi the splitting parameter as defined by Beenakker 1986
// xa  = xi * a
// xa3 = xa * xa * xa
#define mydebug(str,...)
// if (threadIdx.x < 3 && (blockIdx.x == 0 || blockIdx.x == 1)){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
__global__ void sd_compute_mobility_sines(real * r, int N, real * vecs, int num, real * sines, real * cosines, int ldd){
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


// This computes the farfield contribution of the mobility with ewald summation
// this kernel computes the sines and cosinus.
// forces is the vector of [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
// a is the particle radius
// mobility is the mobility matrix which will be retruned (in/out)
// L_d is the boxlength
// cutoff the wavespace cutoff distance
__global__ void sd_wavepart_sum_old(const real * const forces, const real * const vecs,    int num, const  real * const matrices_d,
				const real * const sines,  const real * const cosines, int ldd, int N, real * sin_sum, real * cos_sum, real max){
  // particle index i
  int i = numThreadsPerBlock*(blockIdx.x)+threadIdx.x;
  // in vector on this particle
  real myforce[3];
  real sine;
  real cosine;
  __shared__ real cache[6*numThreadsPerBlock];
  __shared__ real matrices[3*numThreadsPerBlock];
  // get data for myforce - using coalscaled memory access
  for (int l=0;l<3;l++){
    cache[numThreadsPerBlock*l+threadIdx.x] = forces[numThreadsPerBlock*(l+blockIdx.x*3)+threadIdx.x];
  }
  __syncthreads();
  for (int l=0;l<3;l++){
    myforce[l] = cache[threadIdx.x*3+l];
  }
  if (i >= N){
    for (int l=0;l<3;l++){
      myforce[l] = 0;
    }
  }
  int offset_next;
  for (int offset=0;offset<num;offset=offset_next){
    offset_next=offset+numThreadsPerBlock/2;
    // copy matrices to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      matrices[numThreadsPerBlock*l+threadIdx.x] = matrices_d[offset*6+numThreadsPerBlock*l+threadIdx.x];
    }
    __syncthreads();
    for (int j=offset;j< offset_next && j < num;j++){
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	cache[threadIdx.x+k*numThreadsPerBlock]=matrices[(j-offset)*6+k]*myforce[k];
      }
      for (int k=0;k<3;k++){
	assert(!isnan(myforce[k]));
	assert(!isnan(matrices[(j-offset)*6+k]));
	assert(!isnan(cache[threadIdx.x+k*numThreadsPerBlock]));
      }
      sine   = sines  [i + j*ldd];
      cosine = cosines[i + j*ldd];
      assert(!isnan(sine));
      assert(!isnan(cosine));
      cache[threadIdx.x+0*numThreadsPerBlock]+=matrices[(j-offset)*6+3]*myforce[1];
      cache[threadIdx.x+1*numThreadsPerBlock]+=matrices[(j-offset)*6+3]*myforce[0];
      cache[threadIdx.x+0*numThreadsPerBlock]+=matrices[(j-offset)*6+4]*myforce[2];
      cache[threadIdx.x+2*numThreadsPerBlock]+=matrices[(j-offset)*6+4]*myforce[0];
      cache[threadIdx.x+2*numThreadsPerBlock]+=matrices[(j-offset)*6+5]*myforce[1];
      cache[threadIdx.x+1*numThreadsPerBlock]+=matrices[(j-offset)*6+5]*myforce[2];
      for (int k=0;k<6;k++){
	assert(!isnan(cache[threadIdx.x+k*numThreadsPerBlock]));
      }
      for (int k=0;k<DIM;k++){
	cache[threadIdx.x+(k+3)*numThreadsPerBlock]=cosine*cache[threadIdx.x+k*numThreadsPerBlock];
	cache[threadIdx.x+k*numThreadsPerBlock]*=sine;
      }
      for (int k=0;k<6;k++){
	assert(!isnan(cache[threadIdx.x+k*numThreadsPerBlock]));
      }
      reduce_sum(cache);
      reduce_sum(cache+numThreadsPerBlock);
      reduce_sum(cache+numThreadsPerBlock*2);
      reduce_sum(cache+numThreadsPerBlock*3);
      reduce_sum(cache+numThreadsPerBlock*4);
      reduce_sum(cache+numThreadsPerBlock*5);
      if (threadIdx.x < 3){
	real tmp = atomicAdd(sin_sum+threadIdx.x+j*3,cache[numThreadsPerBlock*threadIdx.x]);
	assert(abs(tmp) < max);
	assert(abs(tmp) + cache[numThreadsPerBlock*threadIdx.x] < max);
	tmp = atomicAdd(cos_sum+threadIdx.x+j*3,cache[numThreadsPerBlock*(threadIdx.x+3)]);
	assert(abs(tmp) < max);
	assert(abs(tmp) + cache[numThreadsPerBlock*(threadIdx.x+3)] < max);
      }
    } // for (j = ...
  }// for offset = ...
}

/// This computes the farfield contribution of the mobility with ewald summation
/// this kernel computes the sines and cosinus.
/// The number of threads have to be a power of two! (e.g. 32,64,128 ... )
// forces is the vector of [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
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
      sine   = sines  [blockIdx.x*ldd+threadIdx.x];
      cosine = cosines[blockIdx.x*ldd+threadIdx.x];
      //} else {
      //sine   = 0;
      //cosine = 0;
      //}
      assert(!isnan(sine));
      assert(!isnan(cosine));
      for (int k=0;k<3;k++){
	assert(abs(cache[threadIdx.x*3+k]) < max);
	assert(abs(myforcesin[k]) < max);
	assert(abs(myforcecos[k]) < max);
	assert(!isnan(myforcesin[k]));
	assert(!isnan(myforcecos[k]));
	real tmp=sine*cache[threadIdx.x*3+k];
	assert(abs(tmp) < max);
	myforcesin[k]+=tmp;
	tmp=cosine*cache[threadIdx.x*3+k];
	assert(abs(tmp) < max);
	myforcecos[k]+=tmp;
	assert(!isnan(myforcesin[k]));
	assert(!isnan(myforcecos[k]));
      } // for k
    } // if offset + threadIdx.x < N
  }// for offset = ...
  #pragma unroll 3
  for (int k=0 ; k < 3; k++){
    cache[threadIdx.x+blockDim.x*k]=myforcesin[k];
  }
  reduce_sum(cache);
  reduce_sum(cache+blockDim.x);
  reduce_sum(cache+blockDim.x*2);
  if (threadIdx.x < 3){
    myforcesin[0]=cache[threadIdx.x];
  }
  for (int k=0 ; k < 3; k++){
    cache[threadIdx.x+blockDim.x*k]=myforcecos[k];
  }
  reduce_sum(cache);
  assert(cache[0] < max);
  reduce_sum(cache+blockDim.x);
  assert(cache[blockDim.x] < max);
  reduce_sum(cache+blockDim.x*2);
  assert(cache[blockDim.x*2] < max);
  if (threadIdx.x < 3){
    cache[threadIdx.x]   = cache[threadIdx.x*blockDim.x];
    cache[threadIdx.x+3] = myforcesin[0];
  }
  real * mat = cache + 6;
  if (threadIdx.x < 6){
    assert(cache[threadIdx.x] < max);
    mat[threadIdx.x] = matrices_d[blockIdx.x*6+threadIdx.x];
    cache[threadIdx.x+12]=cache[threadIdx.x]*mat[threadIdx.x%3];
  }
  if (threadIdx.x < 2){
    cache[threadIdx.x*3+1+12]+=mat[3]*cache[threadIdx.x*3+0];
    cache[threadIdx.x*3+0+12]+=mat[3]*cache[threadIdx.x*3+1];
    cache[threadIdx.x*3+2+12]+=mat[4]*cache[threadIdx.x*3+0];
    cache[threadIdx.x*3+0+12]+=mat[4]*cache[threadIdx.x*3+2];
    cache[threadIdx.x*3+2+12]+=mat[5]*cache[threadIdx.x*3+1];
    cache[threadIdx.x*3+1+12]+=mat[5]*cache[threadIdx.x*3+2];
  }
  if (threadIdx.x < 3){
    assert(cache[threadIdx.x+12] < max);
    assert(cache[threadIdx.x+15] < max);
    cos_sum[threadIdx.x+blockIdx.x*3]=cache[threadIdx.x+12];
    sin_sum[threadIdx.x+blockIdx.x*3]=cache[threadIdx.x+15];
  }
}



/// This computes the 
// forces is the vector of [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
// a is the particle radius
// mobility is the mobility matrix which will be retruned (in/out)
// L_d is the boxlength
// cutoff the wavespace cutoff distance
__global__ void sd_wavepart_assemble(const int num, const real * const sines, const real * const cosines, const real * const sin_sum,
				     const real * const cos_sum, const int ldd, real * out, real max){
  int i = threadIdx.x+blockIdx.x*blockDim.x;
  real myout[3];
  real sine;
  real cosine;
#pragma unroll 3
  for (int l=0;l<3;l++){
    myout[l] = 0;
  }
  int offset_next;
  for (int offset=0;offset<num;offset=offset_next){
    offset_next=offset+numThreadsPerBlock;
    if (offset_next > num){
      offset_next=num;
    }
    // copy positions to shared memory
    for (int j=offset;j< offset_next ;j++){
      sine   = sines  [i + j*ldd];
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	real tmp = sin_sum[k+3*j];
	assert(tmp < max);
	myout[k]=sine*tmp;
      }
      cosine = cosines[i + j*ldd];
#pragma unroll 3
      for (int k=0;k<DIM;k++){
	real tmp = cos_sum[k+3*j];
	assert(tmp < max);
	myout[k]+=cosine*tmp;
      }
    } // for (j = ...
  }// for offset = ...
  __shared__ real cache[3*numThreadsPerBlock];
#pragma unroll 3
  for (int k=0;k<3;k++){
    cache[threadIdx.x*3+k]=myout[k];
  }
  __syncthreads();
#pragma unroll 3
  for (int k=0;k<3;k++){
    cache[threadIdx.x+k*numThreadsPerBlock]+=out[threadIdx.x+(blockIdx.x*3+k)*numThreadsPerBlock];
    out[threadIdx.x+(blockIdx.x*3+k)*numThreadsPerBlock]=cache[threadIdx.x+k*numThreadsPerBlock];
  }
}



#define mydebug(str,...)
// if (threadIdx.x < 3 && blockIdx.x < 2){printf("line: %d thread: %2d, block: %2d "str,__LINE__,threadIdx.x,blockIdx.x,__VA_ARGS__);}
// this computes the near field as a  ResistanceMatrix
// r             : is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N             : is the number of particles
// self_mobility : is 1./(6.*PI*eta*a)
// a             : is the particle radius
// L_d           : is the boxlength
// resistance    : is the resistance matrix which will be retruned
// myInfo        : contains infos about the operation:
//                myInfo[0] : number of overlapping particles
//                myInfo[1] : number of interacting particles (via nf)
//                myInfo[2] : max number of interacting particles per particle
__global__ void sd_compute_resistance_matrix(real * pos, int N, real self_mobility, real a, real * L_g, real * resistance, int * myInfo){
  //__shared__ real myPos[3*numThreadsPerBlock];
  int interactions=0;
  real mypos[3];
  __shared__ real L[3];
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  const int lda=(((N*3)+31)/32)*32;
  //__shared__ real myresistance[6*numThreadsPerBlock];
  real myresistance[6]={0,0,0,0,0,0};
  //__shared__ real otherresistance[6*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
  //__syncthreads();
  // get data for myposition - but coalscaled
  /*for (int l=0;l<3;l++){
    myPos[threadIdx.x+l*numThreadsPerBlock] = r[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
    }*/
  for (int l=0;l<3;l++){
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = pos[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
  }

  __syncthreads();
  for (int d=0;d<3;d++){
    mypos[d] = cachedPos[threadIdx.x*3+d];
  }
  
  //for (int i = idx; i < N; i+=blockDim.x*gridDim.x){
  /*if (i < N){
#pragma unroll 3
    for (int k=0; k < DIM; k++){
#pragma unroll 3
      for (int l=0;l < DIM; l++){
	resistance[myindex(DIM*i+k,DIM*i+l)]=0; // we will add some terms on the diagonal, so set it to zero before
      }
    }
  }*/
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
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
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
	  const real safety_fac=0.666; /* this should be one.
					* but setting it to one will result in negative eigenvalues.
					* The Problem is that the correction for the perpendicular
					* entries is larger then the actual value.
					* Reducing the correction is not good - but solves the Problem.
					* This comes mainly from the colesky-decomposition, but is
					* also here to be at least selfconsitant ...
					*/
	  const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	  const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	  const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)            *safety_fac     + 1./6.*log(2.);
	  const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1) *safety_fac     + 1./6.*log(2.);
	  // TODO: Make sure to use (real) and not double ...
	  // real computation
	  real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	  real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	  real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)   * safety_fac   - r2bcorr_perp_self_c ));
	  real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)  * safety_fac   - r2bcorr_perp_mix_c  ));
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
	    
	    //resistance[myindex(DIM*i+k,DIM*i+l)]-=dr[k]*dr[l]*t;
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


__global__ void sd_find_interacting_particles(const real * pos,const real * _L, const int N, int * interacting, int * num_interacting,
					      const int * particle_count, const int * particle_sorted_list, const real * bucket_size_,
					      const int * bucket_num_, const int * particle_to_bucket_list, const real interaction_range,
					      const int total_bucket_num){
  const int lda_short=((N+31)/32)*32;
  const real interaction_range_squared=SQR(interaction_range);
  __shared__ real L[3];
  __shared__ real bucket_size[3];
  __shared__ int bucket_num[6];
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  // copy small vectors to shared memory
  /*{
    const real * tmp;
    if (threadIdx.x < 3){
      tmp=_L;
    } else if (threadIdx.x < 6){
      tmp=bucket_size_-3;
    } else if (threadIdx.x < 9){
      tmp=(real *)bucket_num_;
      tmp-=6;
    }
    if (threadIdx.x < 9){ // wait only once for data
      cachedPos[threadIdx.x]=tmp[threadIdx.x];
    }
    if (threadIdx.x < 3){
      L[threadIdx.x]=cachedPos[threadIdx.x];
      bucket_size[threadIdx.x]=cachedPos[threadIdx.x+3];
      int * tmp2=(int *)(cachedPos+6);
      bucket_num[threadIdx.x]=tmp2[threadIdx.x];
    }
    }*/
  if (threadIdx.x < 3){
    L[threadIdx.x]=_L[threadIdx.x];
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
	
	  assert((current_bucket>=0 && current_bucket < total_bucket_num));
	  int in_loop=particle_count[current_bucket];
	  for (int j_loop=0;j_loop<in_loop;j_loop++){
	    int j=particle_sorted_list[current_bucket+j_loop*total_bucket_num];
	    real dr2=0;
#pragma unroll 3
	    for (int d=0;d<3;d++){
	      real dr=cachedPos[threadIdx.x*3+d]-pos[3*j+d];
	      //dr-=L[d]*rint(dr/L[d]);
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


__global__ void sd_compute_resistance_matrix_sparse(const real * pos, const int N, const real self_mobility,const real a,const real * _L,
  						    real * resistance,const int * col_idx,const int * row_l, int * myInfo){
   const int lda_short=((N+31)/32)*32;
   const int lda = ((N*3+31)/32)*32;
  __shared__ real L[3];
#ifdef SD_USE_FLOAT
  __shared__ real cachedPos[4*numThreadsPerBlock];
#else
  __shared__ real cachedPos[3*numThreadsPerBlock];
#endif
  // copy small vectors to shared memory
  if (threadIdx.x < 3){
    L[threadIdx.x]=_L[threadIdx.x];
  }
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
      dr[d]-=L[d]*rint(dr[d]/L[d]); // fold back
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
	const real safety_fac=0.666; /* this should be one.
				      * but setting it to one will result in negative eigenvalues.
				      * The Problem is that the correction for the perpendicular
				      * entries is larger then the actual value.
				      * Reducing the correction is not good - but solves the Problem.
				      * This comes mainly from the colesky-decomposition, but is
				      * also here to be at least selfconsitant ...
				      */
	const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)            *safety_fac     + 1./6.*log(2.);
	const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1) *safety_fac     + 1./6.*log(2.);
	// TODO: Make sure to use (real) and not double ...
	// real computation
	real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)   * safety_fac   - r2bcorr_perp_self_c ));
	real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)  * safety_fac   - r2bcorr_perp_mix_c  ));
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
__global__ void sd_compute_brownian_force_nearfield(real * r,real * gaussian,int N,real * L_g, real a, real self_mobility,real * brownian_force_nf){
  const int gaussian_ldd=((N+31)/32)*32;
  int interactions=0;
  real mypos[3];
  real writeCache[6];
  //real otherWriteCache[3];
  __shared__ real L[3];
  __shared__ real cachedPos[3*numThreadsPerBlock];
  __shared__ real choleskyCache[12*numThreadsPerBlock];
  //const int lda=(((N*3)+31)/32)*32;
  //__shared__ real myresistance[6*numThreadsPerBlock];
  //real myresistance[6];
  //__shared__ real otherresistance[6*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
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
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
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
      else if (dr2 < 4*a*4*a  && 2*a*2*a < dr2 ){// 2*a < drn < 4*a 
	wasInLoop = 1;
	{
	  real drn= sqrt(dr2); // length of dr
	  real s = drn/a-2;
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
	  const real safety_fac=0.666; /* this should be one.
					* but setting it to one will result in negative eigenvalues.
					* The Problem is that the correction for the perpendicular
					* entries is larger then the actual value.
					* Reducing the correction is not good - but solves the Problem
					*/
	  const real r2bcorr_para_self_c  =                  1./(1.-9./4./dr_c2+3./dr_c4-1./dr_c6) + para_fac_c;
	  const real r2bcorr_para_mix_c   = (6.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)   + para_fac_c;
	  const real r2bcorr_perp_self_c  =                  1./(1.-25./16./dr_c2)            *safety_fac     + 1./6.*log(2.);
	  const real r2bcorr_perp_mix_c   =                  1./(16./20.*dr_c1-25./20./dr_c1) *safety_fac     + 1./6.*log(2.);
	  // TODO: Make sure to use (real) and not double ...
	  // real computation
	  real r2bcorr_para_self     =-( para_fac  + (                      1./(1.-9./4./dr2+3./dr4-1./dr6)  - r2bcorr_para_self_c ));
	  real r2bcorr_para_mix      = ( para_fac  + ( (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.)    - r2bcorr_para_mix_c  ));
	  real r2bcorr_perp_self     =-( perp_fac  + (                      1./(1.-((real)25./16.)/dr2)   * safety_fac   - r2bcorr_perp_self_c ));
	  real r2bcorr_perp_mix      = ( perp_fac  + (                      1./(16./20.*drn-25./20./drn)  * safety_fac   - r2bcorr_perp_mix_c  ));
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
      assert(!isnan(tmp));
    }
  }
}

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

// this adds the identity matrix to a given matrix of ld=size
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
// block : (ignored) the number of elements to process per thread
//         if this is e.g. 3 and the matrix is 3Nx3N, than N threads have to be started
__global__ void sd_add_identity_matrix(real * matrix, int size, int lda){
  //int lda=((size+31)/32)*32;
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //for (int i = idx*block; i< (idx+1)*block; i++){
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    matrix[i+i*lda]+=1;
  }
}

// this sets a block to zero
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
__global__ void sd_set_zero_matrix(real * matrix, int size){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  int matsize=((size+31)/32)*32;
  matsize*=size;
  for (int i = idx;i< matsize; i+=blockDim.x*gridDim.x){
    matrix[i]=0;
  }
}


// this sets a block to zero
// data  : pointer to the given data
// size  : the size of the data
__global__ void sd_set_zero(real * data, int size){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    data[i]=0;
  }
}

// this sets a block to zero
// data  : pointer to the given data
// size  : the size of the data
// value : the value written to the data block
__global__ void sd_set_int(int * data, int size, int value){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    data[i]=value;
  }
}




#define DIST (2+1e-1)
#define DISP_MAX (100)

__global__ void sd_real_integrate_prepare( real * r_d , real * disp_d, real * L, real a, int N){
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
__global__ void sd_real_integrate( real * r_d , real * disp_d, real * L, real a, int N)
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
      tmp-=L[d]*rint(tmp/L[d]);
      dr2+=tmp*tmp;
    }
    if (dr2 <distmin){ // possible colision - check better
      dr2=0;
      //real dr2o=0; // or do we need old distance?
      for (int d=0;d<DIM;d++){
	real tmp=r_d[i*DIM+d]+disp_d[i*DIM+d]-rnew[d];
	tmp-=L[d]*rint(tmp/L[d]);
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
	  t1-=L[d]*rint(t1/L[d]);
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

__global__ void sd_bucket_sort( real * pos , real * bucketSize, int * bucketNum, int N,	int * particleCount,
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

// This function reads the diagonal element entries of a matrix
// and stores them in the vector diag
// and the inverse of them in diag_i
__global__ void sd_get_diag(int size,real * mat_a,int lda,real * diag,real * diag_i){
  const int stepw=lda+1;
  for (int l=threadIdx.x+blockDim.x*blockIdx.x;l<size;l+=blockDim.x*gridDim.x){
    real tmp;
    tmp = mat_a[stepw*l];
    diag[l]=tmp;
    diag_i[l]=1/tmp;
  }
}


__global__ void sd_nrm1(const int size, const real * const vec, real * erg){
  __shared__ real cache[numThreadsPerBlock];
  if (blockIdx.x == 0){ // only use one block
    assert(blockDim.x == numThreadsPerBlock);
    cache[threadIdx.x]=0;
    for (int i = threadIdx.x; i < size ; i+=blockDim.x){
      cache[threadIdx.x]+=abs(vec[i]);
    }
    reduce_sum(cache);
    if (threadIdx.x == 0 ){
      erg[0]=cache[0];
    }
  }
}

