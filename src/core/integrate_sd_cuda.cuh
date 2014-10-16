#pragma once

#ifdef CUDA
#include <cublas_v2.h>

#include "cuda_utils.hpp"
#include "integrate_sd.hpp"

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


/// stuct for sparse matrix: (blocked version off ELLPACK-R, see \link http://elib.uni-stuttgart.de/opus/volltexte/2010/5033/pdf/DIP_2938.pdf \endlink
/// for dense matrizes use only data
struct wavepart{
  real * vecs;
  real * matrices;
  real * cosines;
  real * sines;
  int num;
  int max;
  void _free(){
    if (vecs){
      cuda_safe_mem(cudaFree((void*)vecs));
      vecs=NULL;
    }
    if (matrices){
      cuda_safe_mem(cudaFree((void*)matrices));
      matrices=NULL;
    }
    if (cosines){
      cuda_safe_mem(cudaFree((void*)cosines));
      cosines=NULL;
    }
    if (sines){
      cuda_safe_mem(cudaFree((void*)sines));
      sines=NULL;
    }
    max=0;
  }
};
struct matrix{
  real * data;
  int  * col_idx;
  int  * row_l;
  wavepart * wavespace;
  bool is_sparse;
  int size;
  int ldd;
  int ldd_short;
  void _free(){
    if (data){
      cuda_safe_mem(cudaFree((void*)data));
      data=NULL;
    }
    if (col_idx){
      cuda_safe_mem(cudaFree((void*)col_idx));
      col_idx=NULL;
    }
    if (row_l){
      cuda_safe_mem(cudaFree((void*)row_l));
      row_l=NULL;
    }
    if (wavespace){
      wavespace->_free();
      free((void*)wavespace);
    }
  }
};

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
}

#define SD_RESISTANCE_CORRECT

/* ************************************* *
 * *******   private functions   ******* *
 * ************************************* */
void sd_compute_displacement(cublasHandle_t cublas, real * r_d, int N, real viscosity, real radius, real * L_d, real * L_h,
			     real * total_mobility_d, real * force_d, real * disp_d, int * myInfo);


// this solves iteratively using CG
// disp * (1+resistance*mobility) = mobility_d *  force_d 
// and returnes disp
// mobility and resistance are square matrizes with size <size> and lda <((size+31)/32)*32>
// force and disp are vectors of size <size>
// recalc whether to mobility has changed since the last call or not
real sd_iterative_solver(cublasHandle_t cublas, const real * mobility, const matrix resistance, const real * force, int size, real * disp, real rel_err, real abs_err, bool recalc);

// description below
// in the case of FF_ONLY this computes directly the displacement
real sd_compute_brownian_force_farfield(const matrix mobility, real * gaussian_ff,
					real tol, real * brownian_force_ff );


/// This computes the k-vectors and the wavespace mobility matrices as defined be Beenakker 1986
void sd_compute_mobility_matrix_wave_cpu(int kmax, real kmax2, matrix & mobility, real a, real * L_h, real xi, real selfmobility);




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
int sd_gmres_solver(cublasHandle_t cublas ,int size, const real * A,int lda,
		    real * b, real tol,real abs_tol, int maxit, real * x, real * res);
#else
int sd_gmres_solver(const matrix A1, const matrix A2,
		    real * b, real tol,real abs_tol, int maxit, real * x, real * res);
#endif


// calculates the largest and snalles eigenvalue of the matrix
// mobility    : the mobility matrix (data on the device)      (IN)
// lambda_min  : smalles eigenvalue                            (OUT)
// lambda_max  : largest eigenvalue                            (OUT)
void calculate_maxmin_eigenvalues(const matrix mobility, real * lamba_min, real * lambda_max, real tol);

/// matrix vector multiplication
/// this solves \[ out = factor * A * in \]
/// where A is a matrix (see struct matrix)
/// in and out are vectors
/// factor is a scalar
void sd_Rgemv(const real * factor,const matrix A,const real * in, real * out);


// this function should be fast, as the data should fit (mostly) in L1
// lambda_min   : the lower boundery
// lambda_max   : the upper boundery of the interval
// tol          : the given tollerance which should be achieved
// coefficents  : the pointer where the data will be stored
real calculate_chebyshev_coefficents(real lambda_min, real lambda_max, real tol,real ** coefficents);



#endif
