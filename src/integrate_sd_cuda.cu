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


#include "assert.h"
//#include "helper_functions.cuh"
#include <iostream>
#include "integrate_sd.hpp" // this includes magma and cublas
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "global.hpp"

void _cudaCheckError(const char *msg, const char * file, const int line);
#define cudaCheckError(msg)  _cudaCheckError((msg),__FILE__,__LINE__)

#define myindex(i,j,N) ((i)*(DIM*(N))+(j))


/* ************************************* *
 * *******   private functions   ******* *
 * ************************************* */
void sd_compute_mobility(cublasHandle_t cublas, double * r_d, int N, double eta, double a, double * L_d, double * total_mobility_d);

// This computes the farfield contribution.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L is the boxlength
__global__ void sd_compute_mobility_matrix(double * r, int N, double self_mobility, double a, double * L, double * mobility);


__global__ void sd_add_identity_matrix(double * matrix, int size, int block);
void _cudaCheckError(const char *msg, const char * file, const int line);
// this computes the near field
// it calculates the ResistanceMatrix
__global__ void sd_compute_resistance_matrix(double * r, int N, double self_mobility, double a, double * L, double * resistance);


__global__ void sd_real_integrate( double * r_d , double * disp_d, double * L, double a, int N);




/* ************************************* *
 * *******     implementation    ******* *
 * ************************************* */
void propagate_pos_sd_cuda(double * box_l_h, int N,double * pos_h, double * force_h){
  double viscosity=sd_viscosity;
  double radius   =sd_radius;
    
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
  cuda_safe_mem(cudaMalloc((void**)&pos_d, DIM*N*sizeof(double)));
  cuda_safe_mem(cudaMemcpy(pos_d,pos_h,N*DIM*sizeof(double),cudaMemcpyHostToDevice));
  double * force_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&force_d, DIM*N*sizeof(double)));
  cuda_safe_mem(cudaMemcpy(force_d,force_h,N*DIM*sizeof(double),cudaMemcpyHostToDevice));
  double * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&mobility_d, DIM*DIM*N*N*sizeof(double)));
  double * disp_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&disp_d, DIM*N*sizeof(double)));
  
  sd_compute_mobility(cublas, pos_d, N, viscosity, radius, box_l_d, mobility_d);
  
  //double alpha=1;
  double beta=0;
  cublasStatus_t stat = cublasDgemv( cublas, CUBLAS_OP_T, DIM*N, DIM*N, &time_step, mobility_d, DIM*N, force_d, 1, &beta, disp_d, 1);
  if (stat != CUBLAS_STATUS_SUCCESS) { std::cerr << "CUBLAS Multiplication failed!\n"; errexit(); }
  
  int numThreadsPerBlock = 32;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  //stat = cublasDaxpy(cublas, DIM*N, &alpha, v_d, 1, xr_d, 1);
  //assert(stat==CUBLAS_STATUS_SUCCESS);
  sd_real_integrate<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
  
  
  cuda_safe_mem(cudaMemcpy(pos_h,pos_d,N*DIM*sizeof(double),cudaMemcpyDeviceToHost));
  

  cuda_safe_mem(cudaFree((void*)box_l_d));
  cuda_safe_mem(cudaFree((void*)pos_d));
  cuda_safe_mem(cudaFree((void*)force_d));
  cuda_safe_mem(cudaFree((void*)mobility_d));
  cuda_safe_mem(cudaFree((void*)disp_d));
}



// we try to use as few mallocs as possible, as they seem to be slow ...
void sd_compute_mobility(cublasHandle_t cublas, double * r_d, int N, double eta, double a, double * L_d, double * total_mobility_d){
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("");
  int numThreadsPerBlock = 32;
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  
  // compute the mobility Matrix
  double * helper_d=NULL;
  int ressize=max(4000,DIM*DIM*N*N); // ressize has to be larger for small matrizes (else magma complains)
  cuda_safe_mem(cudaMalloc( (void**)&helper_d, ressize*sizeof(double) ));
  assert(helper_d);
  double * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mobility_d, DIM*DIM*N*N*sizeof(double) ));
  assert(mobility_d);
  sd_compute_mobility_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, mobility_d);
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("compute mobility error");
  // compute the resistance matrix
  double * resistance_d=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&resistance_d, ressize )); //this needs to be bigger for matrix inversion
  assert(resistance_d !=NULL);
  sd_compute_resistance_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*eta*a), a, L_d, resistance_d);
  cudaThreadSynchronize(); // we need both matrices to continue;
  cudaCheckError("compute resistance or mobility error");
  cublasStatus_t status;

  double alpha=1, beta =0;
  status = cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, DIM*N , DIM*N ,DIM*N, &alpha, mobility_d, DIM*N,resistance_d, DIM*N, &beta,helper_d, DIM*N);
  assert(status == CUBLAS_STATUS_SUCCESS);
  sd_add_identity_matrix<<< numBlocks , numThreadsPerBlock  >>>(helper_d,DIM*N,DIM);
  cudaThreadSynchronize();
  
  cudaCheckError("add identity error");
  
  // inverting the matrix 
  int ipiv[DIM*N];
  int info;
  magma_dgetrf_gpu( DIM*N, DIM*N,helper_d, DIM*N, ipiv, &info);
  assert(info==0);
  magma_dgetri_gpu( N*DIM, helper_d, DIM*N, ipiv,resistance_d,ressize, &info);
  assert(info==0);
  // compute the inverse matrix
  // this is an alternative implementation ...
  // be sure to make sure everything else matches, because this one does not overwrite the original matrix
  // GPUGausSeidelDev(helper_d,inverse_d,DIM*N);
  cudaThreadSynchronize();
  cudaCheckError("inversion error");
  // compute total_mobility_d
  status = cublasDgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, DIM*N , DIM*N ,DIM*N, &alpha, \
		       mobility_d, DIM*N,helper_d, DIM*N, &beta,total_mobility_d, DIM*N);
  assert(status == CUBLAS_STATUS_SUCCESS);
  // free the two matrices again
  cudaFree((void*)resistance_d);
  cudaFree((void*)mobility_d);
  cudaFree((void*)helper_d);
  cudaCheckError("in mobility");
}


// This computes the farfield contribution.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L_d is the boxlength
__global__ void sd_compute_mobility_matrix(double * r, int N, double self_mobility, double a, double * L, double * mobility){
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < N){
    int j;
    // get r to sheared memory ?
    // could improve speed
    for (j=0;j<N;j++){
      int k,l;
      if (i==j){
	for (k=0; k < DIM; k++){
	  for (l=0; l < DIM; l++){
	    mobility[myindex(DIM*i+k,DIM*i+l,N)]=0;
	  }
	  mobility[myindex(DIM*i+k,DIM*i+k,N)]=self_mobility;
	}
      }
      else{
	double dr[DIM];
	double dr2=0;
	for (k=0;k<DIM;k++){
	  dr[k]=r[DIM*i+k]-r[DIM*j+k]; // r_ij
	  dr[k]-=rint(dr[k]/L[k])*L[k]; // fold back
	  dr2+=dr[k]*dr[k];
	}
	double drn= sqrt(dr2); // length of dr
	double b = a/drn;
	if (0.5 < b){  // drn < 2*a
	  double t=3./32./drn/a*self_mobility;
	  double t2=(1-9./32.*drn/a)*self_mobility;
	  for (k=0; k < DIM; k++){
	    for (l=0;l < DIM; l++){
	      mobility[myindex(DIM*i+k,DIM*j+l,N)]=dr[k]*dr[l]*t;
	    }
	    mobility[myindex(DIM*i+k,DIM*j+k,N)]+=t2;
	  }
	  // python implementation:
	  //T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
	}
	else{
	  double b2=b*b;
	  double t=(0.75-1.5*b2)*b/dr2*self_mobility;
	  double t2=(0.75+0.5*b2)*b*self_mobility;
	  for (k=0; k < DIM; k++){
	    for (l=0;l < DIM; l++){
	      mobility[myindex(DIM*i+k,DIM*j+l,N)]=dr[k]*dr[l]*t;
	    }
	    mobility[myindex(DIM*i+k,DIM*j+k,N)]+=t2;
	  }
	  // python implementation:
	  // T=one*(0.75+0.5*b2)*b+(0.75-1.5*b2)*b*drt*dr/dr2;
	}
      }
    }
  }
}

// this computes the near field
// it calculates the ResistanceMatrix
__global__ void sd_compute_resistance_matrix(double * r, int N, double self_mobility, double a, double * L, double * resistance){
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < N){ // this should be checked at some point to avoid writing to place we should not write to
    int j;
    // get r to sheared memory ?
    // could improve speed
    for (int k=0; k < DIM; k++){
      for (int l=0;l < DIM; l++){
	resistance[myindex(DIM*i+k,DIM*i+l,N)]=0; // we will add some terms on the diagonal, so set it to zero before
      }
    }
    for (j=0;j<N;j++){
      if (i!=j){ // no self-contribution
	double dr[DIM];
	double dr2=0;
	for (int k=0;k<DIM;k++){
	  dr[k]=r[DIM*i+k]-r[DIM*j+k]; // r_ij
	  dr[k]-=L[k]*rint(dr[k]/L[k]); // fold back
	  dr2+=dr[k]*dr[k];
	}
	if (dr2 < 4*a*4*a && 2*a*2*a < dr2 ){// check whether 2*a < drn < 4*a
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
	  double t=(-0.25/s+9./40.*ls+3./112.*s*ls-t_c)/dr2/self_mobility;
	  double const t2_c=2./6.*log(2.);
	  double t2=(1./6.*ls-t2_c)/self_mobility;
	  for (int k=0; k < DIM; k++){
	    for (int l=0;l < DIM; l++){
	      resistance[myindex(DIM*i+k,DIM*j+l,N)]=dr[k]*dr[l]*t;
	      resistance[myindex(DIM*i+k,DIM*i+l,N)]-=dr[k]*dr[l]*t;
	    }
	    resistance[myindex(DIM*i+k,DIM*j+k,N)]+=t2;
	    resistance[myindex(DIM*i+k,DIM*i+k,N)]-=t2;
	  }
	  // python implementation:
	  //T=one*(1-9./32.*drn/a)+3./32.*dr*drt/drn/a;
	}
	else{ // set the block to zero
	  // it might be faster to set everything in the beginning to zero ...
	  // or use sparse matrices ...
	  for (int k=0; k < DIM; k++){
	    for (int l=0;l < DIM; l++){
	      resistance[myindex(DIM*i+k,DIM*j+l,N)]=0;
	    }
	  }  
	}
      }
    }
  }
}



__global__ void sd_add_identity_matrix(double * matrix, int size, int block){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int i = idx*block; i< (idx+1)*block; i++){
    //#define myindex(i,j,N) ((i)*(DIM*(N))+(j))
    if ( i < size)
      matrix[i+i*size]+=1;
  }
}





void _cudaCheckError(const char *msg, const char * file, const int line)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err)
    {
      std::cerr <<  "Cuda error:" <<  msg << ": " <<  cudaGetErrorString( err) << " in "<<file << "l. "<<line<<"\n";
      errexit();
    }
}





#define DIST (2+1e-6)
#define DISP_MAX (0.5)

__global__ void sd_real_integrate( double * r_d , double * disp_d, double * L, double a, int N)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  // t is the factor how far of disp_d we will move.
  // in case everything is fine, we will move t, if there is some trouble,
  // we will move less to avoid collision
  double t=1;
  if (idx < N){
    double disp2=0;
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
      if (idx==i)
	continue;
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
  }
}








#endif
