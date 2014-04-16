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
//#include "inversion/GPUGausSeidel.h"
//#include "integrate_sd_cuda.cuh"



#include <stdio.h>
#include <iostream>

#include "assert.h"
#include "integrate_sd_cuda_debug.cuh"
#include "integrate_sd.hpp" // this includes magma and cublas
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "global.hpp"

const int numThreadsPerBlock = 32;

void _cudaCheckError(const char *msg, const char * file, const int line);
#define cudaCheckError(msg)  _cudaCheckError((msg),__FILE__,__LINE__)

#define myindex(i,j) ((i)*(lda)+(j))

#define cublasCall(call) stat=(call);assert(stat==CUBLAS_STATUS_SUCCESS)

/* ************************************* *
 * *******   private functions   ******* *
 * ************************************* */
void sd_compute_displacement(cublasHandle_t cublas, double * r_d, int N, double eta, double a, double * L_d, 
			     double * total_mobility_d, double * force_d, double * disp_d);
//void sd_compute_mobility(cublasHandle_t cublas, double * r_d, int N, double eta, double a, double * L_d, double * total_mobility_d);


// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
void sd_iterative_solver(cublasHandle_t cublas, const double * mobility, const double * resistance, const double * force, int size, double * disp);

// This computes the farfield contribution.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L is the boxlength
__global__ void sd_compute_mobility_matrix(double * r, int N, double self_mobility, double a, double * L, double * mobility);


__global__ void sd_add_identity_matrix(double * matrix, int size, int lda);
void _cudaCheckError(const char *msg, const char * file, const int line);
// this computes the near field
// it calculates the ResistanceMatrix
__global__ void sd_compute_resistance_matrix(double * r, int N, double self_mobility, double a, double * L, double * resistance);


__global__ void sd_real_integrate( double * r_d , double * disp_d, double * L, double a, int N);


// this sets a block to zero
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
__global__ void sd_set_zero_matrix(double * matrix, int size);


/* ************************************* *
 * *******     implementation    ******* *
 * ************************************* */

// this calls all the functions to:
//  * generate the mobility matrix (farfield and nearfield)
//  * compute the displacements
//  * add the displacements to the positions
// TODO: add brownian motion, which is current missing
// PARAMTERS:
// box_l_h : the size of the box in x,y and z-direction, on the host (in)
// N       : Number of particles (in)
// pos_h   : position of the particles, simple* array on host (in and out)
// force_h : forces on the particles, simple* array on host (in)
// * : a simple array is e.g. [x_1, y_1, z_1, x_2, y_2, z_2, ...]
void propagate_pos_sd_cuda(double * box_l_h, int N,double * pos_h, double * force_h, double * velo_h){
  //printVectorHost(pos_h,3*N,"pos after call");
  double viscosity=sd_viscosity;
  double radius   =sd_radius;
  if (viscosity  < 0){
    std::cerr << "The viscosity for SD was not set\n";
    errexit();
  }
  if (radius  < 0){
    std::cerr << "The particle radius for SD was not set\n";
    errexit();
  }
  
  int lda=((3*N+31)/32)*32;
  
  static cublasHandle_t cublas=NULL;
  if (cublas==NULL){
    if (cublasCreate(&cublas) != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "CUBLAS initialization failed\n";
      errexit();
    }
    magma_init();
  }

  double * box_l_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&box_l_d, 3*sizeof(double)));
  cuda_safe_mem(cudaMemcpy(box_l_d,box_l_h,3*sizeof(double),cudaMemcpyHostToDevice));
  double * pos_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&pos_d, (DIM)*N*sizeof(double)));
  cuda_safe_mem(cudaMemcpy(pos_d,pos_h,N*DIM*sizeof(double),cudaMemcpyHostToDevice));
  //printVectorDev(pos_d,3*N,"pos after copy");
  double * force_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&force_d, DIM*N*sizeof(double)));
  cuda_safe_mem(cudaMemcpy(force_d,force_h,N*DIM*sizeof(double),cudaMemcpyHostToDevice));
  double * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&mobility_d, lda*N*3*sizeof(double)));
  double * disp_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&disp_d, DIM*N*sizeof(double)));
  
  sd_compute_displacement(cublas, pos_d, N, viscosity, radius, box_l_d, mobility_d, force_d, disp_d);
  
  
  
  //int numThreadsPerBlock = 3;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  //stat = cublasDaxpy(cublas, DIM*N, &alpha, v_d, 1, xr_d, 1);
  //assert(stat==CUBLAS_STATUS_SUCCESS);
  sd_real_integrate<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
  
  // copy back the positions
  cuda_safe_mem(cudaMemcpy(pos_h,pos_d,N*DIM*sizeof(double),cudaMemcpyDeviceToHost));
  // save the displacements as velocities (maybe somebody is interested)
  double alpha=1/time_step;
  cublasStatus_t stat = cublasDscal(cublas, DIM*N, &alpha, disp_d, 1);
  assert(stat == CUBLAS_STATUS_SUCCESS);
  cuda_safe_mem(cudaMemcpy(velo_h,disp_d,N*DIM*sizeof(double),cudaMemcpyDeviceToHost));
  
  
  cuda_safe_mem(cudaFree((void*)box_l_d));
  cuda_safe_mem(cudaFree((void*)pos_d));
  cuda_safe_mem(cudaFree((void*)force_d));
  cuda_safe_mem(cudaFree((void*)mobility_d));
  cuda_safe_mem(cudaFree((void*)disp_d));
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
void sd_compute_displacement(cublasHandle_t cublas, double * r_d, int N, double eta, double a, double * L_d, 
			     double * total_mobility_d, double * force_d, double * disp_d)
{
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("");
  int lda=((3*N+31)/32)*32;
  //int numThreadsPerBlock = 32;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  
  // compute the mobility Matrix
  double * helper_d=NULL;
  int ressize=max(4000,lda*DIM*N); // ressize has to be larger for small matrizes (else magma complains)
  cuda_safe_mem(cudaMalloc( (void**)&helper_d, ressize*sizeof(double) ));
  assert(helper_d);
  double * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mobility_d, lda*DIM*N*sizeof(double) ));
  assert(mobility_d);
  //printMatrixDev(mobility_d,3*N,3*N,"before mobility:");
  //printVectorDev(r_d,3*N,"positions");
  sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(mobility_d,3*N);
  sd_compute_mobility_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, mobility_d);
  cudaThreadSynchronize(); // just for debugging
  //printMatrixDev(mobility_d,3*N,3*N,"mobility_d");
  cudaCheckError("compute mobility error");
  //printMatrixDev(mobility_d,3*N,3*N,"early mobility:");
  // compute the resistance matrix
  double * resistance_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&resistance_d, ressize*sizeof(double) )); //this needs to be bigger for matrix inversion
  assert(resistance_d !=NULL);
  sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(resistance_d,3*N);
  sd_compute_resistance_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, resistance_d);
  cudaThreadSynchronize(); // we need both matrices to continue;
  cudaCheckError("compute resistance or mobility error");
  //cublasStatus_t status;
  
  //debug
  //printMatrixDev(mobility_d,lda,3*N,"late mobility:");
  //printVectorDev(r_d,3*N,"position: ");
  //printMatrixDev(resistance_d,lda,3*N,"resitstance: ");

                                   
  sd_iterative_solver(cublas, mobility_d, resistance_d, force_d, 3*N,disp_d);
    
  

  /*double alpha=1, beta =0;
  status = cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, DIM*N , DIM*N ,DIM*N, &alpha, mobility_d, lda,resistance_d, lda, &beta,helper_d, lda);
  assert(status == CUBLAS_STATUS_SUCCESS);*/
  /*
  sd_add_identity_matrix<<< numBlocks , numThreadsPerBlock  >>>(helper_d,DIM*N,lda);
  cudaThreadSynchronize();
  
  cudaCheckError("add identity error");
  
  // inverting the matrix 
  int ipiv[DIM*N];
  int info;
  magma_dgetrf_gpu( DIM*N, DIM*N,helper_d, lda, ipiv, &info);
  assert(info==0);
  magma_dgetri_gpu( N*DIM, helper_d, lda, ipiv,resistance_d,ressize, &info);
  assert(info==0);
  // compute the inverse matrix
  // this is an alternative implementation ...
  // be sure to make sure everything else matches, because this one does not overwrite the original matrix
  // GPUGausSeidelDev(helper_d,inverse_d,DIM*N);
  cudaThreadSynchronize();
  cudaCheckError("inversion error");
  // compute total_mobility_d
  status = cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, DIM*N , DIM*N ,DIM*N, &alpha, \
		       mobility_d, lda,helper_d, lda, &beta,total_mobility_d, lda);
		       assert(status == CUBLAS_STATUS_SUCCESS);*/
  // free the two matrices again
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
void sd_iterative_solver(cublasHandle_t cublas, const double * mobility, const double * resistance, const double * force, int size, double * disp)
{
  int lda = ((size+31)/32)*32;
  double * mat_a = NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mat_a, lda*size*sizeof(double) ));  assert(mat_a != NULL);
  sd_set_zero_matrix<<<192,32>>>(mat_a,size);
  //double * mob_force=NULL;
  //cuda_safe_mem(cudaMalloc( (void**)&mob_force, size*sizeof(double) ));  assert(mob_force !=NULL);
  // vars for cuBLAS calls
  double alpha=1;
  double beta=0;
  // mat_a = (1+resistance*mobility)
  cublasStatus_t stat;
  cublasCall(cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, size , size ,size, &alpha, mobility, lda,resistance, lda, &beta,mat_a, lda));
  sd_add_identity_matrix<<<128,10>>>(mat_a,size,lda);// TODO: FIXME:  calculate something to set better values ...
  // mob_force = mobility * force
  cublasCall(cublasDgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility, lda, force, 1, &beta, disp, 1));
  int info;
  magma_dpotrs_gpu('U', size, 1, mat_a, lda, disp, size, &info);
  //magma_int_t magma_dpotrs_gpu( magma_uplo_t uplo,  magma_int_t n, magma_int_t nrhs,
  //				double *dA, magma_int_t ldda,
  //				double *dB, magma_int_t lddb, magma_int_t *info);
  if (info != 0){
    fprintf(stderr, "We do not use the symmetry of the matrix ...\n");
    cublasCall(cublasDgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility, lda, force, 1, &beta, disp, 1)); // reset disp
    int ipiv[size];
    magma_dgetrf_gpu( size, size,mat_a, lda, ipiv, &info);
    assert(info==0);
    magma_dgetrs_gpu('N', size, 1,
		     mat_a, lda, ipiv,
		     disp, 1, &info);
    assert(info==0);
    //magma_int_t magma_dgetrs_gpu( magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    //				  double *dA, magma_int_t ldda, magma_int_t *ipiv,
    //				  double *dB, magma_int_t lddb, magma_int_t *info);
  }
  //assert(info==0);

}
// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
void sd_iterative_solver_own(cublasHandle_t cublas, const double * mobility, const double * resistance, const double * force, int size, double * disp)
{
  int lda = ((size+31)/32)*32;
  double * mat_a = NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mat_a, lda*size*sizeof(double) ));  assert(mat_a != NULL);
  sd_set_zero_matrix<<<192,32>>>(mat_a,size);
  double * mob_force=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mob_force, size*sizeof(double) ));  assert(mob_force !=NULL);
  double * resid=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&resid, size*sizeof(double) ));      assert(resid !=NULL);
  double * p=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&p, size*sizeof(double) ));          assert(p !=NULL);
  double * Ap=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&Ap, size*sizeof(double) ));         assert(Ap !=NULL);
  double rs_old;
  // count how many iterations we need
  int counter=0;
  assert(!hasAnyNanDev(mobility,size*lda));
  assert(!hasAnyNanDev(resistance,size*lda));
  // vars for cuBLAS calls
  double alpha=1;
  double beta=0;
  // mat_a = (1+resistance*mobility)
  cublasStatus_t stat;
  cublasCall(cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, size , size ,size, &alpha, mobility, lda,resistance, lda, &beta,mat_a, lda));
  sd_add_identity_matrix<<<128,10>>>(mat_a,size,lda);// TODO: FIXME:  calculate something to set better values ...
  // mob_force = mobility * force
  cublasCall(cublasDgemv( cublas, CUBLAS_OP_T, size, size, &alpha, mobility, lda, force, 1, &beta, mob_force, 1));
  //printMatrixDev(mat_a,lda,size,"A");
  // use mob_force as initial guess
  cublasCall(cublasDcopy(cublas, size,mob_force,1,disp, 1));
  //resid = mob_force-mat_a * disp; //r = b-A*x
  alpha = -1;
  cublasCall(cublasDgemv(cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, disp, 1, &beta, resid, 1));
  printVectorDev(resid,size,"-A*disp");
  printVectorDev(mob_force,size,"solution");
  alpha = 1;
  cublasCall(cublasDaxpy(cublas, size, &alpha, mob_force, 1, resid, 1));
  printVectorDev(resid,size,"residuum");
  
  // p = resid;                     //p=r
  cublasCall(cublasDcopy(cublas, size,resid,1,p, 1));
  // rsquare_old = r * r;           //rsold=r*r
  cublasCall(cublasDdot( cublas, size, resid, 1, resid, 1, &rs_old));
  std::cerr << counter <<" iterations in integrate_sd::inversion, residuum is "<<rs_old<<std::endl;
  const double req_prec=1e-4;
  if (sqrt(rs_old) < req_prec){
    printf("Converged immediatly\n");
    return;
  }
  while (true){
    // Ap = A * p
    beta = 0;  alpha = 1; cublasCall(cublasDgemv(cublas, CUBLAS_OP_T, size, size, &alpha, mat_a, lda, p, 1, &beta, Ap, 1));
    double pAp;
    cublasCall(cublasDdot( cublas, size, p, 1, Ap, 1, &pAp));
    assert(!isnan(pAp));
    //                              //alpha=rsold/pAp
    alpha = rs_old / pAp;
    // disp += alpha * p            // x=x+alpha * p
    cublasCall(cublasDaxpy(cublas, size, &alpha,  p, 1, disp, 1));
    // resid -= alpha * Ap;         // r=r-alpha * Ap
    double minusalpha=-alpha;
    cublasCall(cublasDaxpy(cublas, size, &minusalpha, Ap, 1, resid, 1));
    double rs_new;
    // rs_new = r * r;              // rsnew = r*r
    cublasCall(cublasDdot( cublas, size, resid, 1, resid, 1, &rs_new));
    if (sqrt(rs_new) < req_prec || counter > 2000){
      break;
    }
    // p=resid+rs_new/rs_old*p      // p = r+rsnew/rsold*p
    alpha = rs_new/rs_old;
    cublasCall(cublasDscal( cublas, size, &alpha, p, 1));
    alpha=1;
    cublasCall(cublasDaxpy( cublas, size, &alpha, resid, 1, p, 1));
    //                              // rsold=rsnew;
    rs_old=rs_new;
    counter++;
    if (counter % 100 == 0){
      std::cerr << counter <<" iterations in integrate_sd::inversion, residuum is "<<rs_new<<std::endl;
    }
  }
  printf("Converged after %d iterations\n",counter);
}


// This computes the farfield contribution of the mobility
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L_d is the boxlength
__global__ void sd_compute_mobility_matrix(double * r, int N, double self_mobility, double a, double * L_g, double * mobility){
  double mypos[3];
  const int lda=((3*N+31)/32)*32;
  __shared__ double L[3];
  __shared__ double cachedPos[3*numThreadsPerBlock];
  __shared__ double writeCache[3*numThreadsPerBlock];
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (threadIdx.x < 3){ // copy L to shared memory
    L[threadIdx.x]=L_g[threadIdx.x];
  }
  __syncthreads();
  // get data for myposition - using coalscaled memory access
  for (int l=0;l<3;l++){
    cachedPos[numThreadsPerBlock*l+threadIdx.x] = r[numThreadsPerBlock*l+i];
  }
  __syncthreads();
  for (int l=0;l<3;l++){
    mypos[l] = cachedPos[threadIdx.x*3+l];
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
    // copy positions to shared memory
#pragma unroll
    for (int l=0;l<3;l++){
      cachedPos[threadIdx.x+l*numThreadsPerBlock] = r[numThreadsPerBlock*l+blockIdx.x*blockDim.x*3+threadIdx.x];
    }
    __syncthreads();
    if (i < N){
      for (int j=offset;j<offset+numThreadsPerBlock;j++){
	// this destroys coascaled memory access ...
	/*if (i==j){
	  j++; //just continue with next particle
	  if (j==offset+numThreadsPerBlock){
	    continue;
	  }
	}*/
	/*if (i==j){
#pragma unroll 3
	  for (int l=0;l<3;l++){
	    writeCache[threadIdx.x*3+l]=0;
	  }
	}*/
	if (j < N ){
	  double dr[DIM];
	  double dr2=0;
#pragma unroll 3
	  for (int k=0;k<DIM;k++){
	    //dr[k]=r[DIM*i+k]-r[DIM*j+k]; // r_ij
	    dr[k]=mypos[k]-cachedPos[DIM*(j-offset)+k]; // r_ij
	    /*if (isnan(dr[k])){
	      dr[k]=1337;
	      }*/
	    dr[k]-=rint(dr[k]/L[k])*L[k]; // fold back
	    dr2+=dr[k]*dr[k];
	  }
	  dr2=max(dr2,0.01);
	  double drn= sqrt(dr2); // length of dr
	  double b = a/drn;
      
	  /*if (0.5 < b){  // drn < 2*a
	    /*double t=3./32./drn/a*self_mobility;
	      double t2=(1-9./32.*drn/a)*self_mobility;
	      for (k=0; k < DIM; k++){
	      for (l=0;l < DIM; l++){
	      mobility[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*t;
	      }
	      mobility[myindex(DIM*i+k,DIM*j+k)]+=t2;
	      }*/ // this should not happen ...
	    // python implementation:
	    //T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
	  //}
	  double t,t2;
	  // this also catches the case i == j
	  if (0.5 < b){  // drn < 2*a
	    t=0;
	    t2=0;
	    if (i==j){
	      t2=self_mobility;
	    }
	  } else {
	    double b2=b*b;
	    t=(0.75-1.5*b2)*b/dr2*self_mobility;
	    t2=(0.75+0.5*b2)*b*self_mobility;
	  }
	    
	  //mobility[threadIdx.x]=3+threadIdx.x;
	  double tmp_el13;
#pragma unroll 3 // line 455
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
	    int max = min(blockDim.x, N-offset);
	    for (int l=0;l<3;l++){
	      //mobility[(DIM*j+k)*3*N+blockIdx.x*blockDim.x+threadIdx.x+blockDim.x*l]=writeCache[threadIdx.x+blockDim.x*l];
	      mobility[(DIM*j+k)*lda+offset+threadIdx.x+max*l]=writeCache[threadIdx.x+max*l];
	    }
	    //mobility[myindex(DIM*i+k,DIM*j+k)]+=t2;
	  }
	  // python implementation:
	  // T=one*(0.75+0.5*b2)*b+(0.75-1.5*b2)*b*drt*dr/dr2;
	} // if (j <N)
      } // for (j = ...
    } // if (i < N)
  }// for offset = ...
}


// this computes the near field
// it calculates the ResistanceMatrix
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// L_d is the boxlength
// resistance is the resistance matrix which will be retruned
__global__ void sd_compute_resistance_matrix(double * r, int N, double self_mobility, double a, double * L_g, double * resistance){
  //__shared__ double myPos[3*numThreadsPerBlock];
  double mypos[3];
  __shared__ double L[3];
  __shared__ double cachedPos[3*numThreadsPerBlock];
  const int lda=(((N*3)+31)/32)*32;
  //__shared__ double myresistance[6*numThreadsPerBlock];
  double myresistance[6];
  //__shared__ double otherresistance[6*numThreadsPerBlock];
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
    cachedPos[threadIdx.x+l*numThreadsPerBlock] = r[threadIdx.x+l*numThreadsPerBlock+blockIdx.x*blockDim.x*3];
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
      cachedPos[threadIdx.x+l*numThreadsPerBlock] = r[threadIdx.x+l*numThreadsPerBlock+offset*3];
    }
    __syncthreads();
    for (int j=offset;j<min(offset+numThreadsPerBlock,N);j++){
      double dr[DIM];
      double dr2=0;
#pragma unroll
      for (int k=0;k<DIM;k++){
	dr[k]=mypos[k]-cachedPos[3*(j-offset)+k]; // r_ij
	dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
	dr2+=dr[k]*dr[k];
      }
      double offdiag_fac=0;
      double diag_fac=0;
      if (i >= N || i ==j || j >= N){
	;
      }
      else if (dr2 < 4*a*4*a && 2*a*2*a < dr2 ){// check whether 2*a < drn < 4*a
	// python code:
	// # Use only singular therms, namely to order O(s_ij^0)                                                                  
	// T=(1./4./s-1/4-9./40.*ls)*dr*drt/dr2
	// #           ^ this additonal constant is so that the mobility is smooth
	// # c.f. N.-Q. Nguyen and A. J. C. Ladd, PHYSICAL REVIEW E 66, 046708 (2002) equation (34)                               
	// T+=1./6.*ls*(-one+dr*drt/dr2)
	// R[3*i:3*i+3,3*j:3*j+3]=-T
	// R[3*i:3*i+3,3*i:3*i+3]+=T
	double drn= sqrt(dr2); // length of dr
	double s = drn/a-2;
	double ls = log(s);
	
	double const t_c=-0.125+9./40.*log(2.)+3./112.*2.*log(2.);
	double const t2_c=2./6.*log(2.);
	offdiag_fac =(-0.25/s+9./40.*ls+3./112.*s*ls-t_c)/dr2/self_mobility;
	diag_fac    =(1./6.*ls-t2_c)/self_mobility;
      }
      if (i < N){
#pragma unroll 3
	for (int k=0; k < DIM; k++){
#pragma unroll 3
	  for (int l=0;l < DIM; l++){
	    resistance[myindex(DIM*i+k,DIM*j+l)]=dr[k]*dr[l]*offdiag_fac;
	    //resistance[myindex(DIM*i+k,DIM*i+l)]-=dr[k]*dr[l]*t;
	  }
	  myresistance[k]-=dr[k]*dr[k]*offdiag_fac;
	  resistance[myindex(DIM*i+k,DIM*j+k)]+=diag_fac;
	  myresistance[k]-=diag_fac;
	}
      }
      myresistance[3]-=offdiag_fac*dr[0]*dr[1];
      myresistance[4]-=offdiag_fac*dr[0]*dr[2];
      myresistance[5]-=offdiag_fac*dr[1]*dr[2];
      
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
}


// this adds the identity matrix to a given matrix of ld=size
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
// block : (ignored) the number of elements to process per thread
//         if this is e.g. 3 and the matrix is 3Nx3N, than N threads have to be started
__global__ void sd_add_identity_matrix(double * matrix, int size, int lda){
  //int lda=((size+31)/32)*32;
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //for (int i = idx*block; i< (idx+1)*block; i++){
  for (int i = idx;i< size; i+=blockDim.x*gridDim.x){
    //#define myindex(i,j) ((i)*(DIM*(N))+(j))
    //if ( i < size)
    matrix[i+i*lda]+=1;
  }
}

// this sets a block to zero
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
__global__ void sd_set_zero_matrix(double * matrix, int size){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //for (int i = idx*block; i< (idx+1)*block; i++){
  int matsize=((size+31)/32)*32;
  matsize*=size;
  for (int i = idx;i< matsize; i+=blockDim.x*gridDim.x){
    //if ( i < size)
    matrix[i]=0;
  }
}

// just a dummy kernel ...
// remember to call the other one too ...
__global__ void sd_dummy1(double * L){
  if (threadIdx.x ==0 && blockIdx.x == 0){
    L[0]+=1;
  }
}

__global__ void sd_dummy2(double * L){
  if (threadIdx.x ==0 && blockIdx.x == 0){
    L[0]-=1;
  }  
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





#define DIST (2+1e-1)
#define DISP_MAX (0.5)

__global__ void sd_real_integrate( double * r_d , double * disp_d, double * L, double a, int N)
{
  
  for (int idx = blockIdx.x*blockDim.x + threadIdx.x;
       idx<N ;
       idx+=blockDim.x*gridDim.x){
    // t is the factor how far of disp_d we will move.
    // in case everything is fine, we will move t, if there is some trouble,
    // we will move less to avoid collision
    double t=1;
    double disp2=0;
    //TODO: FIXME: put to separate kernel to avoid race condition
    for (int d=0;d<DIM;d++){
      disp2+=disp_d[idx*DIM+d]*disp_d[idx*DIM+d];
    }
    if (disp2 > DISP_MAX*DISP_MAX){
      double fac=DISP_MAX/sqrt(disp2);
      for (int d=0;d<DIM;d++){
	disp_d[idx*DIM+d]*=fac;
      }
    }
    double rnew[DIM];
    for (int d=0;d<DIM;d++){
      rnew[d]=r_d[DIM*idx+d]+disp_d[DIM*idx+d];
    }
    const double distmin=(3*a)*(3*a);
    for (int i=0;i<N;i++){
      if (idx==i){
	//continue;
	i++;
	if (i >N){
	  continue;
	}
      }
      double dr2=0;
      for (int d=0;d<DIM;d++){
	double tmp=r_d[i*DIM+d]-rnew[d];
	tmp-=L[d]*rint(tmp/L[d]);
	dr2+=tmp*tmp;
      }
      if (dr2 <distmin){ // possible colision - check better
	dr2=0;
	//double dr2o=0; // or do we need old distance?
	for (int d=0;d<DIM;d++){
	  double tmp=r_d[i*DIM+d]+disp_d[i*DIM+d]-rnew[d];
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
	  double alpha=0,beta=0,gamma=0;
	  for (int d=0;d<DIM;d++){
	    double t1=r_d[i*DIM+d]-r_d[idx*DIM+d];
	    t1-=L[d]*rint(t1/L[d]);
	    double t2=disp_d[i*DIM+d]-disp_d[idx*DIM+d];
	    //t2-=L*rint(t2/L); // we would have a problem if we would need to fold back these ...
	    alpha +=t2*t2;
	    beta  +=2*t1*t2;
	    gamma +=t1*t1;
	  } 
	  // now we want to solve for t: alpha*t**2+beta*t+gamma=DIST*a
	  // we want the solution with the minus in the 'mitternachtsformel'
	  // because the other solution is when the particles moved through each other
	  double tnew = (-beta-sqrt(beta*beta-4*alpha*gamma))/(2*alpha);
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
}








#endif
