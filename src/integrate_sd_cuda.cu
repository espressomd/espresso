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

/* This is where the hydro dynamic interaction is implemented */

// TODO:
// * use preconditioner in iterative solver
// * implement matrix-free farfield (via fft)
// * add brownian motion
// * add bucket versions


#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "curand.h"
#include <device_functions.h>
// C++98:
#include <algorithm>
// C++11:
//#include <utility>

#include "assert.h"
#include "integrate_sd_cuda_debug.cuh"
#include "integrate_sd.hpp" // this includes magma and cublas
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "global.hpp"

#ifndef SD_USE_FLOAT
#warning using double
typedef double real;
#define cublasRgemm(...)              cublasDgemm( __VA_ARGS__)
#define cublasRdot(...)               cublasDdot( __VA_ARGS__)
#define cublasRgemv(...)              cublasDgemv( __VA_ARGS__)
#define cublasRcopy(...)              cublasDcopy( __VA_ARGS__)
#define cublasRaxpy(...)              cublasDaxpy( __VA_ARGS__)
#define cublasRscal(...)              cublasDscal( __VA_ARGS__)
#define cublasRnrm2(...)              cublasDnrm2( __VA_ARGS__)
#define curandGenerateNormalReal(...) curandGenerateNormalDouble(__VA_ARGS__)
#define __real2int_rd(...)            __double2int_rd(__VA_ARGS__)
#define rnaupd(...)                   dnaupd_(__VA_ARGS__)
#else // use float
#warning Using Float
typedef float real;
#define cublasRgemm(...)              cublasSgemm( __VA_ARGS__)
#define cublasRdot(...)               cublasSdot( __VA_ARGS__)
#define cublasRgemv(...)              cublasSgemv( __VA_ARGS__)
#define cublasRcopy(...)              cublasScopy( __VA_ARGS__)
#define cublasRaxpy(...)              cublasSaxpy( __VA_ARGS__)
#define cublasRscal(...)              cublasSscal( __VA_ARGS__)
#define cublasRnrm2(...)              cublasSnrm2( __VA_ARGS__)
#define curandGenerateNormalReal(...) curandGenerateNormal(__VA_ARGS__)
#define __real2int_rd(...)            __float2int_rd(__VA_ARGS__)
#define rnaupd(...)                   snaupd_(__VA_ARGS__)
#endif //#ifndef SD_USE_FLOAT


extern double temperature; // this is defined in thermostat.cpp

const int numThreadsPerBlock = 32;

void _cudaCheckError(const char *msg, const char * file, const int line);
#define cudaCheckError(msg)  _cudaCheckError((msg),__FILE__,__LINE__)

#define myindex(i,j) ((i)*(lda)+(j))

#define SQR(x) (x)*(x)

#define cublasCall(call) { cublasStatus_t stat=(call);	\
    assert(stat==CUBLAS_STATUS_SUCCESS);		\
  }
#define curandCall(call) { curandStatus_t stat =(call);	\
    assert(stat == CURAND_STATUS_SUCCESS);		\
  }

// headers for ARPACK-NG: http://forge.scilab.org/index.php/p/arpack-ng/
extern "C"
{
  void dnaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, double* TOL, double RESID[], int* NCV, double V[], int* LDV, int IPARAM[],
	       int IPNTR[], double WORKD[], double WORKL[], int* LWORKL, int* INFO);
  void snaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, float* TOL, float RESID[], int* NCV, float V[], int* LDV, int IPARAM[],
	       int IPNTR[], float WORKD[], float WORKL[], int* LWORKL, int* INFO);
  // this is to compute eigenvectors - but we want only eigenvalues
  //void dneupd_(int* RVEC, char* HOWMNY, int SELECT[], double DR[], double DI[], double Z[], int* LDZ, double* SIGMAR, double* SIGMAI, double WORKEV[],
  //	       char* BMAT, int* N, char WHICH[], int* NEV, double* TOL, double RESID[], int* NCV, double V[], int* LDV, int IPARAM[], int IPNTR[],
  //	       double WORKD[], double WORKL[], int* LWORKL, int* INFO);
}

#define SD_RESISTANCE_CORRECT

/* ************************************* *
 * *******   private functions   ******* *
 * ************************************* */
void sd_compute_displacement(cublasHandle_t cublas, real * r_d, int N, real eta, real a, real * L_d, 
			     real * total_mobility_d, real * force_d, real * disp_d, int * myInfo);


// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
real sd_iterative_solver(cublasHandle_t cublas, const real * mobility, const real * resistance, const real * force, int size, real * disp);

// This computes the farfield contribution.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L is the boxlength
__global__ void sd_compute_mobility_matrix(real * r, int N, real self_mobility, real a, real * L, real * mobility);

// adds to each of the diagonal elemnts of the sizse*size matrix matrix
// with lda lda 1
__global__ void sd_add_identity_matrix(real * matrix, int size, int lda);
void _cudaCheckError(const char *msg, const char * file, const int line);
// this computes the near field
// it calculates the ResistanceMatrix
__global__ void sd_compute_resistance_matrix(real * r, int N, real self_mobility, real a, real * L, real * resistance, int * myInfo);
// TODO: make the order of arguments uniform (and logical?)
// TODO: description here
__global__ void sd_compute_brownian_force_nearfield(real * r,real * gaussian_nf,int N,real * L, real a, real self_mobility,real * brownian_force_nf);

// make sure to have one thread per particle
__global__ void sd_real_integrate_prepare( real * r_d , real * disp_d, real * L, real a, int N);
__global__ void sd_real_integrate( real * r_d , real * disp_d, real * L, real a, int N);


// this sets a block to zero
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
__global__ void sd_set_zero_matrix(real * matrix, int size);



// this sets a block to zero
// data  : pointer to the given data
// size  : the size of the data
__global__ void sd_set_zero(real * data, int size);

// this sets a block to zero
// data  : pointer to the given data
// size  : the size of the data
// value : the value written to the data block
__global__ void sd_set_int(int * data, int size, int value);


// implementation of a bucket sort algorithm
// puts all the N particles with given position pos 
// and particle radius a within the periodic boundary 
// conditions of boxSize L_i = bucketSize_i * bucketNum_i
// puts them in the list particleList
// pos                device array of particle position xyz
// bucketSize         device array with the number of buckets in x y and z direction
// bucketNum          device array with the size of a bucket in x y and z direction
// N                  number of particles
// particleCount      device array of the numbers of particles per bucket. must be initalized to zero
// particleList       device array of the partilces in each bucket
// maxParticlePerCell maximum particles per cell
// totalBucketNUm     bucketNum[0]*bucketNum[1]*bucketNum[2] - the total number of buckets
__global__ void sd_bucket_sort( real * pos , real * bucketSize, int * bucketNum, int N,
				int * particleCount, int * particleList, int maxParticlePerCell, int totalBucketNum);

// BICGSTAB-Solver
// implimented as given in `Numerik linearer Gleichungssysteme` by Prof. Dr. Andreas Meister
// this solves A*x=b
// cublas a handle for cublas
// size   the size n of the matrix
// A      the given n*n matrix (in)
// lda    the leading demension of A
// b      the given solution vector (in)
// tol    requested tolerance of the solution
// maxit  maximum number of iterations
// x      the requested solution with an initial guess (in/out)
// returns 0 on success, else error code
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, real * A,int lda, real * b, real tol, int maxit, real * x, real * res);



// calculates the largest and snalles eigenvalue of the matrix
// size        : size of the eigenvector / the matrix           (IN)
// mobility_d  : handle of the mobility matrix (on the device)  (IN)
// lambda_min  : smalles eigenvalue                            (OUT)
// lambda_max  : largest eigenvalue                            (OUT)
void calculate_maxmin_eigenvalues(int size,real *mobility_d,real * lamba_min,real * lambda_max);


// this function should be fast, as the data should fit (mostly) in L1
// lambda_min   : the lower boundery
// lambda_max   : the upper boundery of the interval
// tol          : the given tollerance which should be achieved
// coefficents  : the pointer where the data will be stored
real calculate_chebyshev_coefficents(real lambda_min, real lambda_max, real tol,real ** coefficents);

typedef unsigned long long ull;
// atomicAdd implementation for double
__device__ double atomicAdd(double * address, double inc);

// global variables for usage in this file
cublasHandle_t cublas=NULL;
/* *************************************************************************************************************** *
 * ********************************************     implementation    ******************************************** *
 * *************************************************************************************************************** */
/* *************************************************************************************************************** *
 * *******     III MM   MM PPP  L     EEEE MM   MM EEEE NN    N TTTTTTT  AAA  TTTTTTT III  OOO  NN    N    ******* *
 * *******      I  M M M M P  P L     E    M M M M E    N N   N    T    A   A    T     I  O   O N N   N    ******* *
 * *******      I  M  M  M PPP  L     EEE  M  M  M EEE  N  N  N    T    AAAAA    T     I  O   O N  N  N    ******* *
 * *******      I  M     M P    L     E    M     M E    N   N N    T    A   A    T     I  O   O N   N N    ******* *
 * *******     III M     M P    LLLL  EEEE M     M EEEE N    NN    T    A   A    T    III  OOO  N    NN    ******* *
 * *************************************************************************************************************** */
/* *************************************************************************************************************** */

/* *************************************************************************************************************** *
 * ********************************************     HOST-Functions    ******************************************** *
 * *************************************************************************************************************** */

// this calls all the functions to:
//  * generate the mobility matrix (farfield and nearfield)
//  * compute the displacements
//  * add the displacements to the positions
// TODO: add brownian motion, which is currently missing
// PARAMTERS:
// box_l_h : the size of the box in x,y and z-direction, on the host (in)
// N       : Number of particles (in)
// pos_h   : position of the particles, simple* array on host (in and out)
// force_h : forces on the particles, simple* array on host (in)
// velo_h  : velocities of the particles, simple* array on host (in and out)
// * : a simple array is e.g. [x_1, y_1, z_1, x_2, y_2, z_2, ...]
void propagate_pos_sd_cuda(real * box_l_h, int N,real * pos_h, real * force_h, real * velo_h){
  //printVectorHost(pos_h,3*N,"pos after call");
  real viscosity=sd_viscosity;
  real radius   =sd_radius;
  if (viscosity  < 0){
    std::cerr << "The viscosity for SD was not set\n";
    errexit();
  }
  if (radius  < 0){
    std::cerr << "The particle radius for SD was not set\n";
    errexit();
  }
  if (time_step < 0){
    std::cerr << "The timestep was not set\n";
    errexit();
  }
  
  int lda=((3*N+31)/32)*32;
  
  //static cublasHandle_t cublas=NULL;
  if (cublas==NULL){
    cublasStatus_t stat = cublasCreate(&cublas);
    if (stat != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "CUBLAS initialization failed in " << __FILE__ << " l. " << __LINE__ <<"\n\t"  ;
      if (stat == CUBLAS_STATUS_NOT_INITIALIZED){
	std::cerr << "the CUDA Runtime initialization failed.\n";
      } else if (stat == CUBLAS_STATUS_ALLOC_FAILED) {
	std::cerr << "the resources could not be allocated\n";
      } else {
	std::cerr << "unknown error\n";
      }
      errexit();
    }
    //magma_init();
  }

  real * box_l_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&box_l_d, 3*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(box_l_d,box_l_h,3*sizeof(real),cudaMemcpyHostToDevice));
  real * pos_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&pos_d, (DIM)*N*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(pos_d,pos_h,N*DIM*sizeof(real),cudaMemcpyHostToDevice));
  //printVectorDev(pos_d,3*N,"pos after copy");
  real * force_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&force_d, DIM*N*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(force_d,force_h,N*DIM*sizeof(real),cudaMemcpyHostToDevice));
  real * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&mobility_d, lda*N*3*sizeof(real)));
  real * disp_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&disp_d, DIM*N*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(disp_d,velo_h,N*DIM*sizeof(real),cudaMemcpyHostToDevice));
  int myInfo_h[]={0,0,0};
  int * myInfo_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&myInfo_d, 3*sizeof(int)));
  cuda_safe_mem(cudaMemcpy(myInfo_d,myInfo_h,3*sizeof(int),cudaMemcpyHostToDevice));
  // rescale forces - this should not be done somewhere else ...
  real alpha=time_step;
  cublasCall(cublasRscal( cublas, 3*N, &alpha, force_d, 1));
  //alpha=1/time_step;
  //cublasCall(cublasRscal(cublas, DIM*N, &alpha, disp_d, 1));  
  sd_compute_displacement(cublas, pos_d, N, viscosity, radius, box_l_d, mobility_d, force_d, disp_d, myInfo_d);
  cuda_safe_mem(cudaMemcpy(myInfo_h,myInfo_d,3*sizeof(int),cudaMemcpyDeviceToHost));
  //std::cerr <<"MyInfo: "<< myInfo_h[0] <<"\t" << myInfo_h[1] <<"\t" << myInfo_h[2] <<"\n";
  
  //int numThreadsPerBlock = 3;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  //stat = cublasRaxpy(cublas, DIM*N, &alpha, v_d, 1, xr_d, 1);
  //assert(stat==CUBLAS_STATUS_SUCCESS);
  alpha=time_step;
  cublasCall(cublasRscal( cublas, 3*N, &alpha, disp_d, 1));  
  sd_real_integrate_prepare<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
  sd_real_integrate<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
  
  // copy back the positions
  cuda_safe_mem(cudaMemcpy(pos_h,pos_d,N*DIM*sizeof(real),cudaMemcpyDeviceToHost));
  // save the displacements as velocities (maybe somebody is interested)
  alpha=1/time_step;
  cublasCall(cublasRscal(cublas, DIM*N, &alpha, disp_d, 1));
  cuda_safe_mem(cudaMemcpy(velo_h,disp_d,N*DIM*sizeof(real),cudaMemcpyDeviceToHost));

  if (myInfo_h[0]){
    ;
    // this needs to be done later, after the data is put where it was in the beginning ...
    //sd_set_particles_apart();
  }
  
  
  cuda_safe_mem(cudaFree((void*)box_l_d));
  cuda_safe_mem(cudaFree((void*)pos_d));
  cuda_safe_mem(cudaFree((void*)force_d));
  cuda_safe_mem(cudaFree((void*)mobility_d));
  cuda_safe_mem(cudaFree((void*)disp_d));
  cuda_safe_mem(cudaFree((void*)myInfo_d));
}



// calculate the farfield and the nearfield and add them
// PARAMETERS:
// cublas : a valid handle of cublas (in)
// r_d    : position of the particles on the device, size 3*N (in)
//          the form has to be [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N      : Number of particles (in)
// eta    : viscositiy of the fluid (in)
// a      : Particle radius (in)
// L_d    : boxsize in x y and z-directions (in)
// total_mobility_d: matrix of the computed total mobility, size 3*3*N*N (in/out, is overwritten)
void sd_compute_displacement(cublasHandle_t cublas, real * r_d, int N, real eta, real a, real * L_d, 
			     real * total_mobility_d, real * force_d, real * disp_d, int * myInfo_d)
{
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("START");
  int lda=((3*N+31)/32)*32;
  //int numThreadsPerBlock = 32;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  
  // compute the mobility Matrix
  real * helper_d=NULL;
  int ressize=max(4000,lda*DIM*N); // ressize has to be larger for small matrizes (else magma complains)
  cuda_safe_mem(cudaMalloc( (void**)&helper_d, ressize*sizeof(real) ));
  assert(helper_d);
  real * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mobility_d, lda*DIM*N*sizeof(real) ));
  assert(mobility_d);
  sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(mobility_d,3*N);
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("sd set zero");
  if(N>32){
    printVectorDev(r_d,96,"pos a: ");
    printVectorDev(r_d+96,N*3-96,"pos b: ");
  }
  else{
    printVectorDev(r_d,3*N,"pos: ");
  }
  sd_compute_mobility_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, mobility_d);
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("compute mobility error");
  // compute the resistance matrix
  real * resistance_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&resistance_d, ressize*sizeof(real) )); //this needs to be bigger for matrix inversion
  assert(resistance_d !=NULL);
  sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(resistance_d,3*N);
  cudaThreadSynchronize(); // debug
  cudaCheckError("sd_set_zero");
  sd_compute_resistance_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, resistance_d, myInfo_d);
  cudaThreadSynchronize(); // we need both matrices to continue;
  cudaCheckError("compute resistance error");
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mobility_d,N*3*lda));
  assert(!hasAnyNanDev(resistance_d,N*3*lda));
  assert(isSymmetricDev(resistance_d,lda,N*3));
#endif
  
  double err = sd_iterative_solver(cublas, mobility_d, resistance_d, force_d, 3*N,disp_d);
#ifdef SD_DEBUG
  if (hasAnyNanDev(disp_d,N*3)){
    printVectorDev(disp_d,N*3,"disp");
    printVectorDev(force_d,N*3,"force");
    printMatrixDev(resistance_d,lda,N*3,"resistance produces nans?");
  }
  assert(!hasAnyNanDev(disp_d,N*3));
#endif
  
  // brownian part
  if (temperature > 0){
    int myInfo_h[3];
    cuda_safe_mem(cudaMemcpy(myInfo_h,myInfo_d,3*sizeof(int),cudaMemcpyDeviceToHost));
    int N_ldd = ((N+31)/32)*32;
    int num_of_rands = N_ldd*myInfo_h[2]*2*DIM+N_ldd*DIM;
    
    real * brownian_force_nf = NULL;
    cuda_safe_mem(cudaMalloc( (void**)&brownian_force_nf, (3*N)*sizeof(real) ));     assert(brownian_force_nf != NULL);
    real * brownian_force_ff = NULL;
    cuda_safe_mem(cudaMalloc( (void**)&brownian_force_ff, (3*N)*sizeof(real) ));     assert(brownian_force_ff != NULL);
    real * gaussian = NULL;
    cuda_safe_mem(cudaMalloc( (void**)&gaussian, (num_of_rands)*sizeof(real) ));     assert(gaussian != NULL);
    real * gaussian_ff = gaussian;
    real * gaussian_nf = gaussian+N_ldd*DIM;
    static curandGenerator_t generator = NULL;
    static int               sd_random_generator_offset=0;
    if (generator == NULL){
      curandCall(curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_DEFAULT));
      curandCall(curandSetPseudoRandomGeneratorSeed(generator, (unsigned long long)('E'+'S'+'P'+'R'+'e'+'s'+'S'+'o')));
      curandCall(curandSetGeneratorOrdering( generator, CURAND_ORDERING_PSEUDO_BEST));
      curandCall(curandSetGeneratorOffset( generator, sd_random_generator_offset));
    }
    //#ifdef FLATNOISE
    // this does not work yet:
    //curandCall(curandGenerateUniformReal(generator, gaussian_d, num_of_rands, 0, sqrt(24.*temperature/time_step)));
    //#else
    curandCall(curandGenerateNormalReal(generator, gaussian, num_of_rands, 0, sqrt(2.*temperature/time_step)));
    //#endif
    if (myInfo_h[2]){
      std::cerr << "Bla!\n";
      int * gamma_index = NULL;
      cuda_safe_mem(cudaMalloc( (void**)&gamma_index, (myInfo_h[2]*N)*sizeof(int) ));    assert(gamma_index != NULL);
      sd_set_zero<<<64,192>>>(brownian_force_nf,3*N);
      sd_set_int<<<64,192>>>(gamma_index,myInfo_h[2]*N, -1);
      sd_compute_brownian_force_nearfield<<<numBlocks, numThreadsPerBlock>>>(r_d, gaussian_nf, N, L_d, a,
									     1./(6.*M_PI*eta*a), brownian_force_nf);
    }// end of near field
    static real * cheby_coefficents=NULL;
    static int N_chebyshev;
    static bool recalc_ew = true;
    int size=3*N;
    real lambda_min, lambda_max;
    if (recalc_ew){
      calculate_maxmin_eigenvalues(3*N,mobility_d,&lambda_min, &lambda_max);
      N_chebyshev = calculate_chebyshev_coefficents(lambda_min, lambda_max,1e-3,&cheby_coefficents);
      recalc_ew=false;
    }
    if (lambda_min < 0){
      printMatrixDev(mobility_d,lda,size,"Mobility has negative eigenvalues!\n");
      errexit();
    }
    real * chebyshev_vec_curr, * chebyshev_vec_last, * chebyshev_vec_next;
    real gaussian_ff_norm;
    sd_set_zero<<<192,192>>>(brownian_force_ff,size);
    cudaThreadSynchronize(); // just for debugging
    cudaCheckError("set zero");
    cublasCall(cublasRnrm2(cublas, size, gaussian_ff, 1, &gaussian_ff_norm));
    chebyshev_vec_curr=gaussian_ff;
    cublasCall(cublasRaxpy( cublas, size, cheby_coefficents+0, chebyshev_vec_curr, 1, brownian_force_ff, 1 ));
    printVectorDev(brownian_force_ff, min(15,size), "l. ?458: ");
    //chebyshev_vec_last=chebyshev_vec_curr;
    chebyshev_vec_last=NULL;
    cuda_safe_mem(cudaMalloc( (void**)&chebyshev_vec_last, size*sizeof(real) ));    assert(chebyshev_vec_last != NULL);
    chebyshev_vec_next=NULL;
    cuda_safe_mem(cudaMalloc( (void**)&chebyshev_vec_next, size*sizeof(real) ));    assert(chebyshev_vec_next != NULL);
    //sd_set_zero<<<192,192>>>(chebyshev_vec_????,size);
    real lambda_minus=lambda_max-lambda_min;
    for (int i=1;i<=N_chebyshev;i++){
      real alpha=2./lambda_minus, beta =0;
      if (i==1){
	alpha=1./lambda_minus;
      }
      cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility_d, lda, chebyshev_vec_curr, 1, &beta, chebyshev_vec_next , 1));
      alpha=-2*(lambda_min+lambda_max)/lambda_minus;
      if (i==1){
	alpha=-(lambda_min+lambda_max)/lambda_minus;
      }
      cublasCall(cublasRaxpy( cublas, size, &alpha, chebyshev_vec_curr, 1, chebyshev_vec_next, 1 ));
      if (i>1){
	alpha=-1;
	cublasCall(cublasRaxpy( cublas, size, &alpha, chebyshev_vec_last, 1, chebyshev_vec_next, 1 ));
      }
      std::swap(chebyshev_vec_curr,chebyshev_vec_next);
      std::swap(chebyshev_vec_last,chebyshev_vec_next);
      cublasCall(cublasRaxpy( cublas, size, cheby_coefficents+i, chebyshev_vec_curr, 1, brownian_force_ff, 1 ));
      real tmp;
      cublasCall(cublasRnrm2(cublas, size, brownian_force_ff, 1 , &tmp));
      fprintf(stderr,"norm brownian force: %e ", tmp);
      cublasCall(cublasRnrm2(cublas, size, chebyshev_vec_last, 1 , &tmp));
      fprintf(stderr,"norm cheby last: %e ", tmp);
      cublasCall(cublasRnrm2(cublas, size, chebyshev_vec_curr, 1 , &tmp));
      fprintf(stderr,"norm cheby curr: %e ", tmp);
      cublasCall(cublasRnrm2(cublas, size, chebyshev_vec_next, 1 , &tmp));
      fprintf(stderr,"norm cheby next: %e \n", tmp);
      //printVectorDev(brownian_force_ff, 15,  "l. ?467:      ");
      //printVectorDev(chebyshev_vec_last, 15, "cheb467: last ");
      //printVectorDev(chebyshev_vec_curr, 15, "cheb467: curr ");
      //printVectorDev(chebyshev_vec_next, 15, "cheb467: next ");
    }
    // errorcheck of chebyshev polynomial
    assert(isSymmetricDev(mobility_d,lda,size));
    real zMz;
    real alpha = 1, beta = 0;
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, mobility_d, lda, brownian_force_ff, 1, &beta,  chebyshev_vec_last, 1));
    cublasCall(cublasRdot(cublas, size, chebyshev_vec_last, 1, brownian_force_ff, 1, &zMz));
    real E_cheby = sqrt(abs(zMz-gaussian_ff_norm*gaussian_ff_norm))/gaussian_ff_norm;
    fprintf(stderr, "The error of the Chebyshev-approximation was %7.3f%% nrm y:%e,  zMz:%e \n",E_cheby*100,gaussian_ff_norm, zMz);
    
    
    
    
  }// end of brownian motion
  cudaCheckError("brownian motion error");
  
  // free everything
  cudaFree((void*)resistance_d);
  cudaFree((void*)mobility_d);
  cudaFree((void*)helper_d);
  cudaCheckError("in mobility");
}

// this calls magma functions to solve the problem: 
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
real sd_iterative_solver(cublasHandle_t cublas, const real * mobility, const real * resistance, const real * force, int size, real * disp)
{
  int lda = ((size+31)/32)*32;
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mobility,size*lda));
  assert(!hasAnyNanDev(resistance,size*lda));
  assert(!hasAnyNanDev(force,size));
#endif
  real * mat_a = NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mat_a, lda*size*sizeof(real) ));       assert(mat_a != NULL);
  real * mat_a_bak = NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mat_a_bak, lda*size*sizeof(real) ));   assert(mat_a_bak != NULL);
  sd_set_zero_matrix<<<192,32>>>(mat_a,size);
  real * mob_force=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mob_force, size*sizeof(real) ));       assert(mob_force !=NULL);
  real * result_checker=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&result_checker, size*sizeof(real) ));  assert(result_checker !=NULL);
  // vars for cuBLAS calls
  real alpha=1;
  real beta=0;
  // mat_a = (1+resistance*mobility)
  cublasCall(cublasRgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, size , size ,size, &alpha, mobility, lda,resistance, lda, &beta,mat_a, lda));
  sd_add_identity_matrix<<<128,10>>>(mat_a,size,lda);// TODO: FIXME:  calculate something to set better values ...
  cuda_safe_mem(cudaMemcpy(mat_a_bak, mat_a, lda*size*sizeof(real),cudaMemcpyDeviceToDevice));
  // mob_force = mobility * force
  cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility, lda, force, 1, &beta, mob_force, 1));
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mat_a,size*lda));
  assert(!hasAnyNanDev(mob_force,size));
#endif
  int info;
  real res;
  //printVectorDev((real *)force,6,"Kraft");
  //printVectorDev(disp,6,"before");
  info = sd_bicgstab_solver(cublas ,size, mat_a,lda, mob_force, 1e-4, 10*size+100, disp, &res);
  //printVectorDev(disp,6,"after");
  // compary to expected result
  //cuda_safe_mem(cudaMemcpy(mat_a, mat_a_bak, lda*size*sizeof(real),cudaMemcpyDeviceToDevice));
  
  if (info != 0){
    if (info == 1){
      if (warnings>1) fprintf(stderr, "Iterative solver did not fully converge ... the residuum was %6e\nWe will continue anyway ...\n",res);
    }
    else{ // info == 2 || info == 4
      // try again with reseted displacement vector as initial guess
      sd_set_zero<<<192,16>>>(disp,size);
      info = sd_bicgstab_solver(cublas ,size, mat_a,lda, mob_force, 1e-4, 10*size+100, disp, &res);
      //printVectorDev(disp,6,"after zeroing");
      if (info == 1){
	if (warnings>1) fprintf(stderr, "Iterative solver did not fully converge ... the residuum was %6e\nWe will continue anyway ...\n",res);
      }
      else if (info == 2){
	if(warnings) fprintf(stderr, "Iterative solver failed ... the residuum was %6e\nWe will continue but the results may be problematic ...\n",res);
      }
    }
    // dgetrs is not better - the contrary: results are worse ...
    /*int ipiv[size];
      magma_dgetrf_gpu( size, size,mat_a, lda, ipiv, &info);
      assert(info==0);
      magma_dgetrs_gpu('N', size, 1,
      mat_a, lda, ipiv,
      disp, size, &info);
      assert(info==0);
      // compary to expected result
      cuda_safe_mem(cudaMemcpy(mat_a, mat_a_bak, lda*size*sizeof(real),cudaMemcpyDeviceToDevice));
      cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, disp, 1, &beta, result_checker, 1));
      alpha=-1;
      cublasCall(cublasRaxpy( cublas, size, &alpha, mob_force, 1, result_checker, 1));
      alpha=1;
      cublasCall(cublasRdot( cublas, size, result_checker, 1, result_checker, 1,&res));
      if (res > 1e-1){
      fprintf(stderr, "All methods failed :(. The residuum from getrs was %e\n",res);
      //cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, disp, 1, &beta, result_checker, 1));
      //printVectorDev(mob_force, size, "mob_force");
      //printVectorDev(result_checker, size, "result_checker");
      //printVectorDev(disp, size, "disp");
      //printMatrixDev((real *)mobility,lda,size,"mobility");
      //printMatrixDev((real *)resistance,lda,size,"res");
      //printMatrixDev((real *)mat_a,lda,size,"mat_a");
      }*/
    //magma_int_t magma_dgetrs_gpu( magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    //				  double *dA, magma_int_t ldda, magma_int_t *ipiv,
    //				  double *dB, magma_int_t lddb, magma_int_t *info);
  }
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(disp,size));
#endif
  //assert(info==0);
  cuda_safe_mem(cudaFree((void*)mat_a));
  cuda_safe_mem(cudaFree((void*)mat_a_bak));
  cuda_safe_mem(cudaFree((void*)mob_force));
  cuda_safe_mem(cudaFree((void*)result_checker));
  return res;
}
// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
void sd_iterative_solver_cg(cublasHandle_t cublas, const real * mobility, const real * resistance, const real * force, int size, real * disp)
{
  int lda = ((size+31)/32)*32;
  real * mat_a = NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mat_a, lda*size*sizeof(real) ));  assert(mat_a != NULL);
  sd_set_zero_matrix<<<192,32>>>(mat_a,size);
  real * mob_force=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mob_force, size*sizeof(real) ));  assert(mob_force !=NULL);
  real * resid=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&resid, size*sizeof(real) ));      assert(resid !=NULL);
  real * p=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&p, size*sizeof(real) ));          assert(p !=NULL);
  real * Ap=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&Ap, size*sizeof(real) ));         assert(Ap !=NULL);
  real rs_old;
  // count how many iterations we need
  int counter=0;
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mobility,size*lda));
  assert(!hasAnyNanDev(resistance,size*lda));
#endif
  // vars for cuBLAS calls
  real alpha=1;
  real beta=0;
  // mat_a = (1+resistance*mobility)
  cublasCall(cublasRgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, size , size ,size, &alpha, mobility, lda,resistance, lda, &beta,mat_a, lda));
  sd_add_identity_matrix<<<128,10>>>(mat_a,size,lda);// TODO: FIXME:  calculate something to set better values ...
  // mob_force = mobility * force
  cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility, lda, force, 1, &beta, mob_force, 1));
  //printMatrixDev(mat_a,lda,size,"A");
  // use mob_force as initial guess
  cublasCall(cublasRcopy(cublas, size,mob_force,1,disp, 1));
  //resid = mob_force-mat_a * disp; //r = b-A*x
  alpha = -1;
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, disp, 1, &beta, resid, 1));
  //printVectorDev(resid,size,"-A*disp");
  //printVectorDev(mob_force,size,"solution");
  alpha = 1;
  cublasCall(cublasRaxpy(cublas, size, &alpha, mob_force, 1, resid, 1));
  //printVectorDev(resid,size,"residuum");
  
  // p = resid;                     //p=r
  cublasCall(cublasRcopy(cublas, size,resid,1,p, 1));
  // rsquare_old = r * r;           //rsold=r*r
  cublasCall(cublasRdot( cublas, size, resid, 1, resid, 1, &rs_old));
  std::cerr << counter <<" iterations in integrate_sd::inversion, residuum is "<<rs_old<<std::endl;
  const real req_prec=1e-4;
  if (sqrt(rs_old) < req_prec){
    printf("Converged immediatly\n");
    return;
  }
  while (true){
    // Ap = A * p
    beta = 0;  alpha = 1; cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, p, 1, &beta, Ap, 1));
    real pAp;
    cublasCall(cublasRdot( cublas, size, p, 1, Ap, 1, &pAp));
    assert(!isnan(pAp));
    //                              //alpha=rsold/pAp
    alpha = rs_old / pAp;
    // disp += alpha * p            // x=x+alpha * p
    cublasCall(cublasRaxpy(cublas, size, &alpha,  p, 1, disp, 1));
    // resid -= alpha * Ap;         // r=r-alpha * Ap
    real minusalpha=-alpha;
    cublasCall(cublasRaxpy(cublas, size, &minusalpha, Ap, 1, resid, 1));
    real rs_new;
    // rs_new = r * r;              // rsnew = r*r
    cublasCall(cublasRdot( cublas, size, resid, 1, resid, 1, &rs_new));
    if (sqrt(rs_new) < req_prec || counter > 2000){
      break;
    }
    // p=resid+rs_new/rs_old*p      // p = r+rsnew/rsold*p
    alpha = rs_new/rs_old;
    cublasCall(cublasRscal( cublas, size, &alpha, p, 1));
    alpha=1;
    cublasCall(cublasRaxpy( cublas, size, &alpha, resid, 1, p, 1));
    //                              // rsold=rsnew;
    rs_old=rs_new;
    counter++;
    if (counter % 100 == 0){
      std::cerr << counter <<" iterations in integrate_sd::inversion, residuum is "<<rs_new<<std::endl;
    }
  }
  printf("Converged after %d iterations\n",counter);
  cuda_safe_mem(cudaFree((void*)mat_a));
  cuda_safe_mem(cudaFree((void*)mob_force));
  cuda_safe_mem(cudaFree((void*)resid));
  cuda_safe_mem(cudaFree((void*)p));
  cuda_safe_mem(cudaFree((void*)Ap));
}

// BICGSTAB-Solver
// implimented as given in Numerik linearer Gleichungssysteme by Prof. Dr. Andreas Meister
// this solves A*x=b
// cublas a handle for cublas
// size   the size n of the matrix
// A      the given n*n matrix (in)
// lda    the leading demension of A
// b      the given solution vector (in)
// tol    requested tolerance of the solution
// maxit  maximum number of iterations
// x      the requested solution with an initial guess (in/out)
// returns 0 on success, else error code
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, real * A,int lda, real * b, real tol, int maxit, real * x, real * res){
  // vector malloc
  real * r0=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&r0, size*sizeof(real) ));       assert(r0 != NULL);
  real * r=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&r, size*sizeof(real) ));        assert(r != NULL);
  real * p=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&p, size*sizeof(real) ));        assert(p != NULL);
  real * v=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&v, size*sizeof(real) ));        assert(v != NULL);
  real * t=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&t, size*sizeof(real) ));        assert(t != NULL);
  real * test=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&test, size*sizeof(real) ));     assert(test != NULL);
  // constants
  real eps;
  if (sizeof(real) == sizeof(double)){
    eps = 1e-15;
  } else {
    eps = 1e-7;
  }
  eps = min(eps,tol*1e-2);
  // other variables
  real alpha=1;
  real beta=0;
  real tolb;
  // compute the norm of b
  real normb;
  cublasCall(cublasRdot( cublas, size, b, 1, b, 1, &normb));
  normb=sqrt(normb);
  //tolb=min(tol*size, tol*normb); // tol is not realy usefull as this wont be reached ... at least without preconditioning
  //tolb=max(normb*eps, tolb);
  tolb=tol*normb;
  // r0 = b-A*x
  alpha = -1;
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, A, lda, x, 1, &beta, r0, 1));
  alpha = 1;
  cublasCall(cublasRaxpy(cublas, size, &alpha, b, 1, r0, 1));
  // r = r0
  cublasCall(cublasRcopy(cublas, size,r0,1,r, 1));
  // rr0 = r*r0
  real rr0;
  cublasCall(cublasRdot( cublas, size, r0, 1, r0, 1, &rr0));
  // p =r
  cublasCall(cublasRcopy(cublas, size,r0,1,p, 1));
  // normr=norm(r)
  real normr=sqrt(rr0);
  int iteration=0;
  real lastnorm=normr;
  real initnorm=normr;
  // check for conversion or max iterations
  while (iteration < maxit && normr >= tolb){
    // v=A*p
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, A, lda, p, 1, &beta, v, 1));
    // vr0 = v*r0
    real vr0;
    cublasCall(cublasRdot( cublas, size, v, 1, r0, 1, &vr0));
    if (fabs(vr0) < eps || rr0 == 0){
      if (fabs(vr0) < eps){
	if (warnings > 1) fprintf(stderr, "BICGSTAB break-down.\n");
      }else{
	if (warnings > 1) fprintf(stderr, "BICGSTAB solution stagnates.\n");
      }
      cuda_safe_mem(cudaFree((void*)r0));cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));cuda_safe_mem(cudaFree((void*)test));
      res[0] = normr;
      if (tolb > normr)     { return 0;}
      if (tolb*100 > normr) { return 1;}
      else                  { return 2;}
    }
    // alpha = rr0/vr0
    real myAlpha=rr0/vr0;
    real minusMyAlpha = -myAlpha;
    // s = r - alpha v
    //cublasCall(cublasRcopy(cublas, size,r,1,s, 1));
    cublasCall(cublasRaxpy(cublas, size, &minusMyAlpha, v, 1, r, 1)); //s->r
    // t = A * s
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, A, lda, r, 1, &beta, t, 1));// s->r
    // ts = s * t
    real ts;
    cublasCall(cublasRdot( cublas, size, t, 1, r, 1, &ts));// s->r
    // tt = t * t
    real tt;
    cublasCall(cublasRdot( cublas, size, t, 1, t, 1, &tt));
    if (abs(tt)<eps || ts == 0){
      fprintf(stderr,"Exit: abs(tt)<eps || ts == 0\n");
      if (warnings > 1) fprintf(stderr, "BICGSTAB break-down.\n");
      cuda_safe_mem(cudaFree((void*)r0));cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));cuda_safe_mem(cudaFree((void*)test));
      res[0] = normr;
      if (tolb > normr)     { 
	//fprintf(stderr, "0: tolb: %e normr: %e \n",tolb, normr);
	return 0;}
      if (tolb*100 > normr) { 
	//fprintf(stderr, "1: tolb: %e normr: %e \n",tolb, normr);
	return 1;}
      else                  { 
	//fprintf(stderr, "2: tolb: %e normr: %e \n",tolb, normr);
	return 2;}
    }
    // omega = ts/tt
    real myOmega=ts/tt;
    // x = x + alpha p + omega s
    cublasCall(cublasRaxpy(cublas, size, &myAlpha, p, 1, x, 1));
    cublasCall(cublasRaxpy(cublas, size, &myOmega, r, 1, x, 1));
    // copyback of s to r
    // r = s - omega t
    real minusMyOmega=-1*myOmega;
    cublasCall(cublasRaxpy(cublas, size, &minusMyOmega, t, 1, r, 1));
    //myOmega*=-1;
    // r1r0 = r * r0
    real r1r0;
    cublasCall(cublasRdot( cublas, size, r, 1, r0, 1, &r1r0));
    // beta = (alpha * r1r0 ) / (omega rr0)
    real myBeta = (myAlpha*r1r0)/(myOmega*rr0);
    if (abs(myBeta)>1/eps){
      fprintf(stderr,"Exit: abs(myBeta)<1/eps\n");
      cuda_safe_mem(cudaFree((void*)r0));cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));cuda_safe_mem(cudaFree((void*)test));
      res[0] = normr;
      if (tolb > normr)     { return 0;}
      if (tolb*100 > normr) { return 1;}
      else                  { return 2;}
    }
    // p = r + beta ( p - omega v)= beta p + r - beta omega v
    cublasCall(cublasRscal(cublas, size, &myBeta, p, 1));
    cublasCall(cublasRaxpy(cublas, size, &alpha, r, 1, p, 1));
    alpha=-myBeta*myOmega;
    cublasCall(cublasRaxpy(cublas, size, &alpha, v, 1, p, 1));
    alpha=1;
    rr0=r1r0;
    real r1r1;
    cublasCall(cublasRdot( cublas, size, r, 1, r, 1, &r1r1));
    normr=sqrt(r1r1);
    iteration++;
    if (lastnorm*sqrt(eps) > normr){ // restart
      //fprintf(stderr, "recalculation r\n");
      cublasCall(cublasRcopy(cublas, size,b,1,r, 1));
      alpha=-1;beta=1;
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, A, lda, x, 1, &beta, r, 1));
      alpha= 1;beta=0;
      cublasCall(cublasRdot( cublas, size, r, 1, r, 1, &rr0));
      normr=sqrt(rr0);
      lastnorm = normr;
      // r = r0
      cublasCall(cublasRcopy(cublas, size,r,1,r0, 1));
      // p =r
      cublasCall(cublasRcopy(cublas, size,r,1,p, 1));
    }
    if (iteration%500000 == 0){ // enable debugging by setting this to a lower value
      real realnorm;
      {// recalculate normr
	cublasCall(cublasRcopy(cublas, size,b,1,test, 1));
	alpha=-1;beta=1;
	cublasCall(cublasRgemv(cublas, CUBLAS_OP_T, size, size, &alpha, A, lda, x, 1, &beta, test, 1));
	alpha= 1;beta=0;
	cublasCall(cublasRdot( cublas, size, test, 1, test, 1, &realnorm));
	realnorm=sqrt(realnorm);
      }
      fprintf(stderr,"  Iteration: %6d Residuum: %12f RealResiduum: %12f\n",iteration, normr, realnorm);
    }
    if (initnorm*1e10 < normr){ // somehow our solution explodes ...
      if (warnings) fprintf(stderr, "BICGSTAB did not converge, residuum exploded. Aborting.\n");
      cuda_safe_mem(cudaFree((void*)r0));cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));cuda_safe_mem(cudaFree((void*)test));
      res[0]= normr;
      return 4;
    }
  }
  res[0]=normr;
  if (normr > tolb*1.01){
    fprintf(stderr, "BICGSTAB solution did not converge after %d iterations. Error was %e1 %% to high.\n",iteration,(normr/tolb-1)*100);
    cuda_safe_mem(cudaFree((void*)r0));cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));cuda_safe_mem(cudaFree((void*)test));
    if (tolb*100 > normr){
      fprintf(stderr, "1: tolb: %e normr: %e \n",tolb, normr);
      return 1;
    } else {
      fprintf(stderr, "2: tolb: %e normr: %e \n",tolb, normr);
      return 2;
    }
  }
  //fprintf(stderr, "BICGSTAB solution did converge after %d iterations.\n",iteration);
  
  cuda_safe_mem(cudaFree((void*)r0));
  cuda_safe_mem(cudaFree((void*)r));
  cuda_safe_mem(cudaFree((void*)p));
  cuda_safe_mem(cudaFree((void*)v));
  cuda_safe_mem(cudaFree((void*)t));
  cuda_safe_mem(cudaFree((void*)test));
  return 0;
}

// calculates the largest and snalles eigenvalue of the matrix
// size        : size of the eigenvector / the matrix           (IN)
// mobility_d  : handle of the mobility matrix (on the device)  (IN)
// lambda_min  : smalles eigenvalue                            (OUT)
// lambda_max  : largest eigenvalue                            (OUT)
void calculate_maxmin_eigenvalues(int size,real *mobility_d,real * lambda_min,real * lambda_max){
  int lda = ((size+31)/32)*32;
  int maxit=max(500,size);
  int IDO;
  char BMAT='I'; // standard eigenvalue problem
  char WHICH[]="SR"; // start with largest eigenvalue
  int NEV = 1; // only one eigenvalue
  // TUNING: these could be adjusted?
  real TOL=1e-1;
  if (sizeof(double) == sizeof(real)){
    TOL=max(1e-12,TOL);
  } else {
    TOL=max(1e-6,TOL);
  }
  // TUNING: make some tests to find good value ...
  int NCV=min(size, 6); // must be at least 3, but a bit bigger should be better ...
  int LDV=lda;
  int LWORKL=3*NCV*(NCV + 2);
  int mode=1;
  int IPARAM[11] = {1,0,maxit,1,0,0,mode,0,0,0,0};
  int IPNTR[14];
  int INFO=0;
  real RESID[size];
  real V[LDV*NCV];
  real WORKD[3*size];
  real WORKL[LWORKL];
  real * vec_in_d;
  real * vec_out_d;
  cuda_safe_mem(cudaMalloc((void**)&vec_in_d , lda*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**)&vec_out_d, lda*sizeof(real)));
  for (int minmax=0;minmax<2;minmax++){
    IDO=0;
    if (minmax){
      sprintf(WHICH,"LR");
      INFO=1;
      IPARAM[2]=maxit;
      TOL=1e-3;
    }
    while (IDO != 99){
      //dnaupd_(&IDO,&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD,WORKL,&LWORKL,&INFONAUP);
      rnaupd(&IDO,&BMAT,&size, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
      switch (IDO){
      case 1:
	cuda_safe_mem(cudaMemcpy(vec_in_d,WORKD+IPNTR[0]-1,size*sizeof(real),cudaMemcpyHostToDevice));
	{
	  real alpha=1, beta=0;
	  cublasCall(cublasRgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility_d, lda, vec_in_d, 1, &beta, vec_out_d, 1));
	}
	cuda_safe_mem(cudaMemcpy(WORKD+IPNTR[1]-1,vec_out_d,size*sizeof(real),cudaMemcpyDeviceToHost));
	break;
      case -1:
      case 2:
      case 3:
      case 4:
	fprintf(stderr,"Error in %s l. %d: unexpected work from rnaupd: %d: Not Implemented!\n",__FILE__,__LINE__,IDO);
	break;
      case 99: //we are done
	break;
      default:
	fprintf(stderr,"Error in %s l. %d: unexpected work from rnaupd: %d: Not Understood!\n",__FILE__,__LINE__,IDO);
	break;
      }
    } 
    fprintf(stderr,"calculationg eigenvalue needed %d iterations and %d gemv operations (tolerance is %e, EW is %e).\n"
	    ,IPARAM[2], IPARAM[8], TOL,WORKL[IPNTR[5]-1]);
    if (INFO){
      fprintf(stderr,"Unexpected return value in %s l. %d from rnaupd_: %d\n",__FILE__,__LINE__,INFO);
    }
    if (WORKL[IPNTR[5]-1]<0 && TOL > 1e-3){
      minmax--;
      TOL=1e-4;
      INFO=1;
      IPARAM[2]=maxit*100;
    }
    if (minmax){ // make them a bit larger/smaller to be sure that we are in the interval of interrest ...
      *lambda_max=WORKL[IPNTR[5]-1]*(1+TOL);
    } else {
      *lambda_min=WORKL[IPNTR[5]-1]*(1-TOL);
    }
  }
  /* FORTRAN Comments ...
     c          IPNTR(6): pointer to the real part of the ritz value array     
     c                    RITZR in WORKL.                                          
     c          IPNTR(7): pointer to the imaginary part of the ritz value array    
     c                    RITZI in WORKL.                                          
     c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
     c                    with the Ritz values located in RITZR and RITZI in WORK
  */
}

// lambda_min   : the lower boundery
// lambda_max   : the upper boundery of the interval
// tol          : the given tollerance which should be achieved
// coefficents  : the pointer where the data will be stored
real calculate_chebyshev_coefficents(real lambda_min, real lambda_max, real tol,real ** coefficents){
  // use chebyshev-gausquadrature: https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature
  const int steps=1024*128; // with 1024 it should fit in L1
  unsigned int N=1024; // guess the number of coefficents we need, if more are needed -> realocate
  if (*coefficents==NULL){
    *coefficents = (real *)malloc(N*sizeof(real));
  } else {
    *coefficents = (real *)realloc(*coefficents,N*sizeof(real));
  }
  real * current_polynome=NULL;
  real * last_polynome=NULL;
  last_polynome    = (real *) malloc(steps * sizeof(real));
  current_polynome = (real *) malloc(steps * sizeof(real));
  real x[steps];
  real weight_and_func[steps];
  real lambda_m=lambda_max-lambda_min;
  real lambda_p=lambda_max+lambda_min;
  //fprintf(stderr,"lambda_min: %e lambda_max: %e  ", lambda_min, lambda_max);
  //fprintf(stderr,"lambda_minus: %e lambda_plusminus: %e\n", lambda_m, lambda_pm);
  // init
  real fac=2./steps;
  //fprintf(stderr,"fac: %e\n", fac);
  double ai=0;
  double a1=0;
  for (int i=0;i<steps;i++){
    last_polynome[i]=1; //T0(x)=1
    //real x = -1.+(i*2.+1.)/steps;
    x[i]=cos((2.*i+1)/2./steps*M_PI);
    current_polynome[i]= x[i];  //T1(x)=x
    weight_and_func[i]=fac*1./sqrt(x[i]*lambda_m/2.+lambda_p/2.);// /sqrt(1-x*x);// this could be big, but should not be inf
    ai+=weight_and_func[i];
    a1+=weight_and_func[i]*x[i];
    if (i%100 == 0){
      ;//fprintf(stderr,"i: %d a0: %e  a1: %e \n",i,ai,a1);
    }
  }
  real error;
  int loop=0;
  //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,0,ai);
  //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,1,a1);
  (*coefficents)[loop]=ai/2.;
  loop++;
  (*coefficents)[loop]=a1;
  //double sumfacmax=lambda_max;
  //double totalsum=abs(ai)+abs(a1)*sumfacmax;
  double totalsum=abs(ai)+abs(a1);
  //sumfacmax*=lambda_max;
  const int miniloop=10;
  do{
    error=0;
    do {
      std::swap(last_polynome,current_polynome);
      //{	real * tmp       = last_polynome; last_polynome    = current_polynome; current_polynome = tmp; }
      //printf("addresses: 0x%08x  0x%08x ",last_polynome,current_polynome);
      ai=0;
      for (int i=0;i<steps;i++){
	//real x = -1.+(i*2.+1.)/steps;
	current_polynome[i]=-1.*current_polynome[i]+((real)2)*x[i]*last_polynome[i];
	//printf ("%e %e %e\n", x, weight_and_func[i], current_polynome[i]);
	ai+=current_polynome[i]*weight_and_func[i];
      }
      //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,loop,ai);
      //printf("\n\n");
      //sumfacmax*=lambda_max;
      error+=abs(ai);//*sumfacmax;
      totalsum+=abs(ai);//*sumfacmax;
      (*coefficents)[loop]=ai;
      loop++;
    } while (loop%miniloop);
    if (loop+miniloop > N){
      N*=2;
      *coefficents=(real *)realloc(*coefficents,N*sizeof(real));
    }
  } while ((error > tol*totalsum || loop < 20 ) && loop < sqrt(steps));
  if (loop >=steps/10 -1 ){
    fprintf(stderr,"to few steps to get sufficent results in %s l. %d\n",__FILE__,__LINE__);
  }
  error=0;
  while (error < tol*totalsum){ // approximate error
    loop--;
    error+=abs((*coefficents)[loop]);//*sumfacmax;
    //sumfacmax/=lambda_max;
  }
  fprintf(stderr,"sum: %e   error: %e",totalsum,error);
  loop++;
  free(last_polynome);
  free(current_polynome);
  fprintf(stderr,"loops: %d\n",loop);
  return loop;
}


// check whether there was any cuda error so far.
// do not use this function directly but use the macro cudaCheckError(const char *msg);
// which requires only the first paramter
// PARAMTERS:
// msg   : the message which should be printed in case of an error
// file  : the file in which the function is called
// line  : the line in which the function is called
void _cudaCheckError(const char *msg, const char * file, const int line)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err)
    {
      std::cerr <<  "Cuda error:" <<  msg << ": '" <<  cudaGetErrorString( err) << "' in "<<file << " l. "<<line<<"\n";
      errexit();
    }
}


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
	  dr[k]-=rint(dr[k]/L[k])*L[k]; // fold back
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
	if (0.5 < b){  // drn < 2*a
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
#warning "Wrong Mobility in Farfield - to assure positive definiteness"
	  // the /5. is to much ... but like this it seems to work
	  t=(0.75-1.5*b2)*b/dr2*self_mobility/5.;
	  t2=(0.75+0.5*b2)*b*self_mobility/5.;
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
    //mydebug("pos: 0x%010x   offset: 0x%010x\n",pos,threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3);
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
      mydebug("dr2: %f\n",dr2);
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
	  
	  mydebug("ls: %e \n",ls);
#ifdef SD_RESISTANCE_CORRECT
	  real const t_c=-0.125+9./40.*log(2.)+3./112.*2.*log(2.);
	  real offdiag_fac =(-0.25/s+9./40.*ls+3./112.*s*ls-t_c)/dr2;
	  real diag_fac    =(1./6.*ls);
#else
	  real const t_c=-0.125+9./40.*log(2.)+3./112.*2.*log(2.);
	  real const t2_c=2./6.*log(2.);
	  offdiag_fac =(-0.25/s+9./40.*ls+3./112.*s*ls-t_c)/dr2/self_mobility;
	  diag_fac    =(1./6.*ls-t2_c)/self_mobility;
#endif
#ifdef SD_RESISTANCE_CORRECT
	  real dr4=dr2*dr2;
	  real dr6=dr4*dr2;
	  // constants for correction
	  const real dr_c1 = 4;
	  const real dr_c2 = 4*4;
	  const real dr_c3 = 4*4*4;
	  const real dr_c4 = 4*4*4*4;
	  const real dr_c5 = 4*4*4*4*4;
	  const real dr_c6 = 4*4*4*4*4*4;
	  const real r2bcorr_diag_self_c    = (4.*dr_c6)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)         ;
	  const real r2bcorr_diag_mix_c     = (9.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.);
	  const real r2bcorr_offdiag_self_c = 16.*dr_c2 /(16.*dr_c2-25)                            - 2./6.*log(2.);
	  const real r2bcorr_offdiag_mix_c  = 20.*dr_c1 /(16.*dr_c2-25)                            - 2./6.*log(2.);
	  // real computation
	  r2bcorr_diag_self     = diag_fac    - 1./(1-9./4./dr2+3./dr4-1./dr6)                     + r2bcorr_diag_self_c;
	  r2bcorr_diag_mix      = diag_fac    - (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.) + r2bcorr_diag_mix_c;
	  r2bcorr_offdiag_self  = offdiag_fac - 1./(1.-25./16./dr2)                                + r2bcorr_offdiag_self_c;
	  r2bcorr_offdiag_mix   = offdiag_fac - 1./(16./20.*drn-25./20./drn)                       + r2bcorr_offdiag_mix_c;
	  r2bcorr_diag_self    /= self_mobility;
	  r2bcorr_diag_mix     /= self_mobility;
	  r2bcorr_offdiag_self /= self_mobility;
	  r2bcorr_offdiag_mix  /= self_mobility;
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
	  mydebug(""
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
    for (int j=offset;j<min(offset+numThreadsPerBlock,N);j++){
      real dr[DIM];
      real dr2=0;
#pragma unroll
      for (int k=0;k<DIM;k++){
	dr[k]=mypos[k]-cachedPos[3*(j-offset)+k]; // r_ij
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
	dr2+=dr[k]*dr[k];
      }
#ifndef SD_RESISTANCE_CORRECT
#warning "SD Brownian motion only support corrected resistance calculation ..."
#endif
      real r2bcorr_diag_self     = 0;
      real r2bcorr_diag_mix      = 0;
      real r2bcorr_offdiag_self  = 0;
      real r2bcorr_offdiag_mix   = 0;

      int wasInLoop = 0;
      if (i >= N || i >= j || j >= N){
	writeCache[3]=0;
	writeCache[4]=0;
	writeCache[5]=0;
      }
      // j > i
      else if (dr2 < 4*a*4*a  && 2*a*2*a < dr2 ){// 2*a < drn < 4*a 
	wasInLoop = 1;
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
	
	real const t_c=-0.125+9./40.*log(2.)+3./112.*2.*log(2.);
	real offdiag_fac =(-0.25/s+9./40.*ls+3./112.*s*ls-t_c)/dr2;
	real diag_fac    =(1./6.*ls);
	
	real dr4=dr2*dr2;
	real dr6=dr4*dr2;
	// constants for correction
	const real dr_c1 = 4;
	const real dr_c2 = 4*4;
	const real dr_c3 = 4*4*4;
	const real dr_c4 = 4*4*4*4;
	const real dr_c5 = 4*4*4*4*4;
	const real dr_c6 = 4*4*4*4*4*4;
	const real r2bcorr_diag_self_c    = (4.*dr_c6)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.)         ;
	const real r2bcorr_diag_mix_c     = (9.*dr_c5-4.*dr_c3)/(4.*dr_c6-9.*dr_c4+12.*dr_c2-4.);
	const real r2bcorr_offdiag_self_c = 16.*dr_c2 /(16.*dr_c2-25)                            - 2./6.*log(2.);
	const real r2bcorr_offdiag_mix_c  = 20.*dr_c1 /(16.*dr_c2-25)                            - 2./6.*log(2.);
	// real computation
	r2bcorr_diag_self     = diag_fac    - 1./(1-9./4./dr2+3./dr4-1./dr6)                     + r2bcorr_diag_self_c;
	r2bcorr_diag_mix      = diag_fac    - (6.*dr4*drn-4.*dr2*drn)/(4.*dr6-9.*dr4+12.*dr2-4.) + r2bcorr_diag_mix_c;
	r2bcorr_offdiag_self  = offdiag_fac - 1./(1.-25./16./dr2)                                + r2bcorr_offdiag_self_c;
	r2bcorr_offdiag_mix   = offdiag_fac - 1./(16./20.*drn-25./20./drn)                       + r2bcorr_offdiag_mix_c;
	r2bcorr_diag_self    /= self_mobility;
	r2bcorr_diag_mix     /= self_mobility;
	r2bcorr_offdiag_self /= self_mobility;
	r2bcorr_offdiag_mix  /= self_mobility;
	
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
	myCC(0)  = sqrt(r2bcorr_diag_mix +r2bcorr_offdiag_mix *dr[0]*dr[0] - SQR(    myCC(3))- SQR(    myCC(7))
			- SQR(     myCC(1)));
	myCC(2)  =                       (r2bcorr_offdiag_mix *dr[0]*dr[1] - myCC(4)*myCC(3) - myCC(8)*myCC(7) 
					  - myCC(10)*myCC(1))/myCC(0);
	myCC(6)  =                       (r2bcorr_offdiag_mix *dr[0]*dr[2] - myCC(5)*myCC(3) - myCC(9)*myCC(7) 
					  - myCC(11)*myCC(1))/myCC(0);
	writeCache[3]+=myCC(0)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[4]+=myCC(2)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(6)  * gaussian[3*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 18 - 6
	// L_{5,5} and L_{6,5}
	myCC(0)  = sqrt(r2bcorr_diag_mix +r2bcorr_offdiag_mix *dr[1]*dr[1] - SQR(    myCC(4))- SQR(    myCC(8))
			- SQR(     myCC(10))- SQR(    myCC(2)));
	myCC(3)  =                       (r2bcorr_offdiag_mix *dr[1]*dr[2] - myCC(5)*myCC(4) - myCC(9)*myCC(8) 
					  - myCC(11)*myCC(10) - myCC(6)*myCC(2))/myCC(0);
	writeCache[4]+=myCC(0)  * gaussian[4*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	writeCache[5]+=myCC(3)  * gaussian[4*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used: 20 - 10
	// L_{6,6} would be:
	myCC(0) = sqrt(r2bcorr_diag_mix +r2bcorr_offdiag_mix *dr[2]*dr[2] - SQR(myCC(5))    - SQR(myCC(9))     
		       - SQR(myCC(11)) - SQR(myCC(6)) - SQR(myCC(3)));
	writeCache[5]+=myCC(0)  * gaussian[5*gaussian_ldd+threadIdx.x+blockDim.x*blockIdx.x+6*gaussian_ldd*interactions];
	// used 21 - 15
	interactions++;
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
	  atomicAdd(brownian_force_nf+j*3,   choleskyCache[0]);
	  atomicAdd(brownian_force_nf+j*3+1, choleskyCache[2*numThreadsPerBlock]);
	  atomicAdd(brownian_force_nf+j*3+2, choleskyCache[4*numThreadsPerBlock]);
	}
      }
    }
  }
  if ( i < N){
#pragma unroll 3
    for (int k=0;k<3;k++){
      atomicAdd(brownian_force_nf+i*3+k, writeCache[k]);
    }
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
#define DISP_MAX (10000)

__global__ void sd_real_integrate_prepare( real * r_d , real * disp_d, real * L, real a, int N){
  /*for (int idx = blockIdx.x*blockDim.x + threadIdx.x;
       idx<N ;
       idx+=blockDim.x*gridDim.x){
    real disp2=0;
#pragma unroll
    for (int d=0;d<DIM;d++){
      disp2+=disp_d[idx*DIM+d]*disp_d[idx*DIM+d];
    }
    if (disp2 > DISP_MAX*DISP_MAX){
      real fac=DISP_MAX/sqrt(disp2);
#pragma unroll
      for (int d=0;d<DIM;d++){
	disp_d[idx*DIM+d]*=fac;
      }
    }
  }*/
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
  
  //for (int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //     idx<N ;
  //     idx+=blockDim.x*gridDim.x){
  int idx =  blockIdx.x*blockDim.x + threadIdx.x;
  // t is the factor how far of disp_d we will move.
  // in case everything is fine, we will move t, if there is some trouble,
  // we will move less to avoid collision
  real t=1;
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
  for (int d=0;d<DIM;d++){ // actually do the integration
    r_d[DIM*idx+d]+=disp_d[DIM*idx+d]*t;
  }
  //#warning "Debug is still enabaled"
    //pos_d[DIM*N+idx]=t;
}

__global__ void sd_bucket_sort( real * pos , real * bucketSize, int * bucketNum, int N,
				int * particleCount, int * particleList, int maxParticlePerCell, int totalBucketNum){
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
      // this should work - but somehow it does not compile
      x=__real2int_rd(tmp);
      // the following code is an replacement ...
      // but with this the loop is not getting unrolled
      //asm("cvt.rmi.s32.f64 %0, %1;\n"
      //    : "=r"(x) : "d"(tmp) : );
      // this should also work.
      // but the corresponding ptx code first rounds, and then converts in a second step ...
      // this could lead to rounding errors ...
      //x=floor(tmp);
      //x%=bucketNum[d];
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
    int num = atomicAdd(particleList+myBucket, 1);
    if (num < maxParticlePerCell){ // every thread should do this - so this is not a branch ...
      particleList[myBucket+num*totalBucketNum]=i;
    }else{
      // Note: printf in device code works only with cc>=2.0 //
#if (__CUDA_ARCH__>=200)
      printf("error: overflow in grid cell (%i,%i,%i)\n",bucket.x,bucket.y,bucket.z);
#endif
    }
  }
}


/* *************************************************************************************************************** *
 * ********************************************    DEVICE-Functions   ******************************************** *
 * *************************************************************************************************************** */

__device__ double atomicAdd(double * address, double inc){

ull *addressUll = (ull*) address;
ull oldValue=*addressUll;
ull assumedValue;
do {
assumedValue=oldValue;
ull newValue = __double_as_longlong (__longlong_as_double(assumedValue)+inc);
oldValue = atomicCAS(addressUll,assumedValue,newValue);
}
  while (oldValue != assumedValue);
return __longlong_as_double(oldValue);
}

#endif
