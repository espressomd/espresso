#ifndef EWALDGPUFORCE_HPP
#define EWALDGPUFORCE_HPP
#include "config.hpp"

#ifdef EWALD_GPU

#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "cells.hpp"

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define CUDA_USE_SINGLE_PRECISION
#endif

#include <cuda.h>
#include <cuda_runtime.h>

// the following works around an incompatibility between Eigen and the nvcc preprocessor
#define EIGEN_DONT_VECTORIZE
#include "VectorForce.hpp"

#ifdef ELECTROSTATICS_GPU_DOUBLE_PRECISION
typedef double ewaldgpu_real;
#else
typedef float ewaldgpu_real;
#endif

#ifndef INT64
#define INT64	long long
#endif
#ifndef POW2
#define POW2(val) ((val)*(val))
#endif

class EwaldgpuForce : public VectorForce {
public:
	//KONSTRUCTOR
	EwaldgpuForce(double rcut, int num_kx, int num_ky, int num_kz, double alpha);
	~EwaldgpuForce();
	//BASIC FUNCTIONS
	void init(SystemInterface &systemInterface); // Called by ForceIterator to initalize
	void run(SystemInterface &systemInterface); // Calculate forces via GPU
	bool isReady(); // Wait until all force computations via GPU are done
  void runEnergies(SystemInterface &s); // Calculate forces via GPU
  bool isReadyEnergies(); // Wait until all energy computations via GPU are done
  void destroy(); // Free memory etc
	bool m_ready; // All computations are done

protected:
	//SYSTEM
  double m_box_l[3];  // box length
  double m_V; // Volume
  double m_coulomb_prefactor; // Forces, energy is multiplied with this factor
  int   m_N; // number of particles
  ewaldgpu_real *m_r_i; // particle positions;
  ewaldgpu_real *m_dev_r_i; // particle positions;
  ewaldgpu_real *m_q_i; // charges of the particles
  ewaldgpu_real *m_dev_q_i; // charges of the particles
  ewaldgpu_real *m_k; // k-vectors
  ewaldgpu_real *m_dev_k; // k-vectors
  ewaldgpu_real *m_forces_reci; // Forces k-space
  ewaldgpu_real *m_dev_forces_reci; // Forces k-space
  ewaldgpu_real *m_infl_factor; // Influence factor
  ewaldgpu_real *m_dev_infl_factor; // Influence factor
  ewaldgpu_real *m_rho_hat; // Rho hat
  ewaldgpu_real *m_dev_rho_hat; // Rho hat
  ewaldgpu_real *m_energy_reci; // energy k-space
  ewaldgpu_real *m_dev_energy_reci; // energy k-space
  ewaldgpu_real *m_energy_self; // energy self
  ewaldgpu_real m_q_sqr; // sum of the squares of the particle charges
	//PARAMETERS
  double m_alpha; // separation parameter
  double m_rcut; // cutoff radius
  int m_num_kx, m_num_ky, m_num_kz, m_num_k; // Number of k's in xyz-direction
  //STREAM
  cudaEvent_t *start, *stop;
  cudaStream_t    *stream0;

protected:
  //VARIABLES
  bool initialized; //Memory etc has been initialized
  bool isTuned; // Tuning is over
  SystemInterface *System;
	//COMPUTE K's
	void compute_num_k(); // Compute the number of k's in sphere respectively in ellipsoid
	void compute_k_AND_influence_factor(); // Compute the k-vectors in sphere respectively in ellipsoid and influence factor
	//GPU PROGRAM
	void runGPU_Forces(); // Run GPU forces part
	void runGPU_Energy(); // Run GPU forces part
	int nextPow2(int x); // Determine the next power of x
	bool isPow2(int x); // Determine if x of power 2
	void getNumBlocksAndThreads(int Size, int maxBlocks, int maxThreads, int &blocks, int &threads); // Determine the number of blocks and threads in GPU part
	//TRANSFORM EIGEN VECTORS
	void Transform_Vector3d_to_Array(SystemInterface &systemInterface); // Store the positions, currents, types etc in Vecror3d format from TCL in arrays used in GPU part
	void Transform_ForceArray_to_Vector3d(); // Store the calculated forces in Vector3d format used in TCL
	//OUTPUT
	void Output(); // Output in terminal
	//REAL SPACE
	void EwaldCPU_EnergySelf(); // Calculate the self energy
	void compute_q_sqare(); // Calculate sum q_i**2
};

typedef struct {
  Cell *cell;
  Particle *p;
  int i,c,np;

	double rcut;
	int num_kx;
	int num_ky;
	int num_kz;
	double alpha;
	bool ewaldgpu_is_running;
	//Tuning
	double accuracy;
	double precision;
	bool isTuned; // Flag that tuning is over
	int K_max; // Maximal reciprocal K-vector in tuning
	int time_calc_steps; // Steps in time_force_calc function
} Ewaldgpu_params;

extern Ewaldgpu_params ewaldgpu_params;

//Add forces and energy
inline void add_ewald_gpu_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
  int j;
  double fac;
  double rcut=ewaldgpu_params.rcut;
  if(dist < rcut)
  {
		fac=coulomb.prefactor * p1->p.q * p2->p.q * (  2*ewaldgpu_params.alpha/sqrt(M_PI) * exp(-POW2(ewaldgpu_params.alpha *dist)) + erfc(ewaldgpu_params.alpha*dist)/dist )/POW2(dist);

    for(j=0;j<3;j++)
      force[j] += fac * d[j];

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: EWALD_GPU   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: EWALD_GPU   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}
inline double ewaldgpu_coulomb_pair_energy(double chgfac, double *d,double dist2,double dist)
{
  if (dist < ewaldgpu_params.rcut)
  {
     return chgfac*erfc(ewaldgpu_params.alpha*dist)/dist;
  }
  return 0.0;
}
//Set parameters
int ewaldgpu_set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha);
int ewaldgpu_set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps);
//TUNING r_cut, num_kx, num_ky, num_kz, alpha
int ewaldgpu_adaptive_tune(char **log);
double ewaldgpu_error_estimate_r(double q_sqr, int N, double r_cut, double V, double alpha, double accuracy);
double ewaldgpu_error_estimate_k(double q_sqr, int N, int K, double V, double alpha, double accuracy);
double ewaldgpu_tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N);
double ewaldgpu_tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N);
int ewaldgpu_count_charged_particles();
//Kolaffa compute optimal alpha
double ewaldgpu_compute_E_error_estimate_r(double alpha, double rcut, double q_sqr, double box_l[3]);
double ewaldgpu_compute_E_error_estimate_k(double alpha, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3]);
double ewaldgpu_E_estimate_error(double rcut, int num_kx, int num_ky, int num_kz, double alpha, double q_sqr, double box_l[3]);
double ewaldgpu_compute_optimal_alpha(double rcut, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3], double precision);
double ewaldgpu_compute_q_sqare();

#endif //EWALD_GPU
#endif //EWALDGPUFORCE_HPP
