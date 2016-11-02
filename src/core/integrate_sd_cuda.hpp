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

#ifndef __INTEGRATE_SD_CUDA_HPP
#define __INTEGRATE_SD_CUDA_HPP

#include "integrate_sd.hpp"

#ifdef SD
#ifdef CUDA
#include <cublas_v2.h>

void errexit();

#ifndef SD_USE_FLOAT
typedef double real;
#define cublasRgemm(...)                cublasDgemm( __VA_ARGS__)
#define cublasRdot(...)                 cublasDdot( __VA_ARGS__)
//#define cublasRgemv(...)                cublasDgemv( __VA_ARGS__) // always use sd_Rgemv
#define cublasRcopy(...)                cublasDcopy( __VA_ARGS__)
#define cublasRaxpy(...)                cublasDaxpy( __VA_ARGS__)
#define cublasRscal(...)                cublasDscal( __VA_ARGS__)
//#define cublasRnrm2(...)              cublasDnrm2( __VA_ARGS__)
#define cublasRnrm2(a1, a2, a3, a4, a5) cublasDdot(a1,a2, a3, a4, a3, a4, a5));*(a5)=sqrt(*(a5)
#define cublasRtrsv(...)                cublasDtrsv( __VA_ARGS__)
#define curandGenerateNormalReal(...)   curandGenerateNormalDouble(__VA_ARGS__)
#define __real2int_rd(...)              __double2int_rd(__VA_ARGS__)
#define __real2int_ru(...)              __double2int_ru(__VA_ARGS__)
#define rnaupd(...)                     dnaupd_(__VA_ARGS__)
#else // use float
typedef float real;
#define cublasRgemm(...)                cublasSgemm( __VA_ARGS__)
#define cublasRdot(...)                 cublasSdot( __VA_ARGS__)
//#define cublasRgemv(...)                cublasSgemv( __VA_ARGS__) // always use sd_Rgemv
#define cublasRcopy(...)                cublasScopy( __VA_ARGS__)
#define cublasRaxpy(...)                cublasSaxpy( __VA_ARGS__)
#define cublasRscal(...)                cublasSscal( __VA_ARGS__)
//#define cublasRnrm2(...)              cublasSnrm2( __VA_ARGS__)
#define cublasRnrm2(a1, a2, a3, a4, a5) cublasSdot(a1,a2, a3, a4, a3, a4, a5));*(a5)=sqrt(*(a5)
#define cublasRtrsv(...)                cublasStrsv( __VA_ARGS__)
#define curandGenerateNormalReal(...)   curandGenerateNormal(__VA_ARGS__)
#define __real2int_rd(...)              __float2int_rd(__VA_ARGS__)
#define __real2int_ru(...)              __float2int_ru(__VA_ARGS__)
#define rnaupd(...)                     snaupd_(__VA_ARGS__)
#endif //#ifndef SD_USE_FLOAT

#include "cuda_utils.hpp"
#include "integrate_sd.hpp"
#include "integrate_sd_cuda.hpp"
#include "integrate_sd_matrix.hpp"

extern cublasHandle_t cublas;

extern double temperature; ///< this is defined in \file thermostat.cpp

const int numThreadsPerBlock = 32; ///< gives the number of threads per cuda block. should be a multiple of 32 and much smaller then the numbers of particles
#define numThreadsPerBlock_is_power_of_two TRUE

#ifndef SD_FF_ONLY
const int bucket_size_factor = 3;  //TODO: find optimal parameter for given implementation
#endif

void _cudaCheckError(const char *msg, const char * file, const int line);
#define cudaCheckError(msg)  cudaThreadSynchronize();_cudaCheckError((msg),__FILE__,__LINE__)

#define myindex(i,j) ((i)*(lda)+(j))

#define SQR(x) (x)*(x)

#define cublasCall(call) { cublasStatus_t stat=(call);			\
    if (stat != CUBLAS_STATUS_SUCCESS) {				\
      std::cerr << "CUBLAS failed in " << __FILE__			\
		<< " l. " << __LINE__ <<"\n\t"  ;			\
      if (stat == CUBLAS_STATUS_NOT_INITIALIZED){			\
	std::cerr << "the CUDA Runtime initialization failed.\n";	\
      } else if (stat == CUBLAS_STATUS_ALLOC_FAILED) {			\
	std::cerr << "the resources could not be allocated\n";		\
      } else {								\
	std::cerr << "unknown error\n";					\
      }									\
      errexit();							\
    }									\
  }
#define curandCall(call) { curandStatus_t stat =(call);			\
    if (stat != CURAND_STATUS_SUCCESS) {				\
      std::cerr << "CURAND failed in " << __FILE__			\
		<< " l. " << __LINE__ <<"\n\t"  ;			\
      if (stat == CURAND_STATUS_NOT_INITIALIZED){			\
	std::cerr << "the CUDA Runtime initialization failed.\n";	\
      } else if (stat == CURAND_STATUS_VERSION_MISMATCH) {		\
	std::cerr << "Header file and linked library version do not match.\n"; \
      } else if (stat == CURAND_STATUS_ALLOCATION_FAILED) {		\
	std::cerr << "Memory allocation failed.\n";			\
      } else if (stat == CURAND_STATUS_TYPE_ERROR) {			\
	std::cerr << "Generator is wrong type.\n";			\
      } else if (stat == CURAND_STATUS_OUT_OF_RANGE) {			\
	std::cerr << "the argument was out of range.\n";		\
      } else if (stat == CURAND_STATUS_LENGTH_NOT_MULTIPLE) {		\
	std::cerr << "Length requested is not a multple of dimension.\n"; \
      } else if (stat == CURAND_STATUS_DOUBLE_PRECISION_REQUIRED) {	\
	std::cerr << "GPU does not have double precision required by MRG32k3a.\n"; \
      } else if (stat == CURAND_STATUS_LAUNCH_FAILURE) {		\
	std::cerr << "Kernel launch failure.\n";			\
      } else if (stat == CURAND_STATUS_PREEXISTING_FAILURE) {		\
	std::cerr << "Preexisting failure on library entry.\n";		\
      } else if (stat == CURAND_STATUS_INITIALIZATION_FAILED) {		\
	std::cerr << "Initialization of CUDA failed.\n";		\
      } else if (stat == CURAND_STATUS_ARCH_MISMATCH) {			\
	std::cerr << "Architecture mismatch, GPU does not support requested feature.\n"; \
      } else if (stat == CURAND_STATUS_INTERNAL_ERROR) {		\
	std::cerr << "Internal library error.\n";\
      } else {								\
	fprintf(stderr,"unknown error (error number: %d)\n",stat);	\
      }									\
      errexit();							\
    }									\
  }

// headers for ARPACK-NG: http://forge.scilab.org/index.php/p/arpack-ng/
extern "C"
{
  void dnaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, double* TOL, double RESID[], int* NCV, double V[], int* LDV, int IPARAM[],
	       int IPNTR[], double WORKD[], double WORKL[], int* LWORKL, int* INFO);
  void snaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, float* TOL, float RESID[], int* NCV, float V[], int* LDV, int IPARAM[],
	       int IPNTR[], float WORKD[], float WORKL[], int* LWORKL, int* INFO);
}

#define SD_RESISTANCE_CORRECT

/* ************************************* *
 * *******   private functions   ******* *
 * ************************************* */
void sd_compute_displacement(const real * r_d, int N, real viscosity, real radius, const real * L_d, const real * L_h,
			     const real * total_mobility_d, const real * force_d, real * disp_d, int * myInfo);


// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
// recalc whether to mobility has changed since the last call or not
real sd_iterative_solver(const matrix & mobility, const matrix & resistance, const real * force, real * disp, real rel_err, real abs_err, bool recalc);

// description below
// in the case of FF_ONLY this computes directly the displacement
real sd_compute_brownian_force_farfield(const matrix & mobility, const real * gaussian_ff,
					real tol, real * brownian_force_ff );


/// This computes the k-vectors and the wavespace mobility matrices as defined be Beenakker 1986
void sd_compute_mobility_matrix_wave_cpu(int kmax, real kmax2, matrix & mobility, real a, const real * L_h, real xi, real selfmobility, real ew_prec);




void _cuda_check_error(char *file, unsigned int line);
#define cuda_check_error(); _cuda_check_error(__FILE__,__LINE__);
// BICGSTAB-Solver
// implimented as given in `Numerik linearer Gleichungssysteme` by Prof. Dr. Andreas Meister
// additional a preconditioner can be used (compile-time-switch). This overwrites the soluten vector b
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
/*#ifdef SD_PC
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, const real * A,                   int lda,
		       real * b, real tol,real abs_tol,  int maxit, real * x, real * res);
#else
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, const real * A1, const real * A2, int lda,
		       real * b, real tol,real abs_tol,  int maxit, real * x, real * res);
		       #endif*/

/// GMRes-Solver
// implimented as given in Numerik linearer Gleichungssysteme by Prof. Dr. Andreas Meister
// additional a preconditioner can be used (compile-time-switch). This overwrites the soluten vector b
// this solves A*x=b
// cublas a handle for cublas
// size   the size n of the matrix
// A      the given n*n matrix (in) (to use with preconditioning)
// A1,A2  the matrizes that give A=A1*A2+1 (to use without preconditioning)
// lda    the leading demension of A
// b      the given solution vector (in)
// tol    requested tolerance of the solution
// maxit  maximum number of iterations
// x      the requested solution with an initial guess (in/out)
// returns 0 on success, else error code
#ifdef SD_PC
int sd_gmres_solver(int size, const real * A, int lda,
		    real * b, real tol,real abs_tol, int maxit, real * x, real * res);
#else
int sd_gmres_solver(const matrix & A1, const matrix & A2,
		    const real * b, real tol,real abs_tol, int maxit, real * x, real * res);
#endif


// calculates the largest and snalles eigenvalue of the matrix
// mobility    : the mobility matrix (data on the device)      (IN)
// lambda_min  : smalles eigenvalue                            (OUT)
// lambda_max  : largest eigenvalue                            (OUT)
void calculate_maxmin_eigenvalues(const matrix & mobility, real * lamba_min, real * lambda_max, real tol);

/// matrix vector multiplication
/// this solves \[ out = factor * A * in \]
/// where A is a matrix (see struct matrix)
/// in and out are vectors
/// factor is a scalar
void sd_Rgemv(const real * factor,const matrix & A,const real * in, real * out);


// compute the chebyshev coeffeicents which are needed to achieve a given accuracy.
// lambda_min   : the lower boundery
// lambda_max   : the upper boundery of the interval
// tol          : the given tollerance which should be achieved
// coefficents  : the pointer where the data will be stored
real calculate_chebyshev_coefficents(real lambda_min, real lambda_max, real tol,real ** coefficents);

// find the maximal value of an array on the GPU
// data         : the array
// length       : number of elements
// return value : the maximal value
int sd_find_max(int * data, int length);

#endif /* SD */
#endif /* CUDA */

#endif
