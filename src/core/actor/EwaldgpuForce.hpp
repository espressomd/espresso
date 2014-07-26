#ifndef EWALDGPUFORCE_HPP
#define EWALDGPUFORCE_HPP
#include "config.hpp"

#ifdef EWALD_GPU

#include "SystemInterface.hpp"
#include "Actor.hpp"
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "cells.hpp"

#ifndef INT64
#define INT64	long long
#endif
#ifndef POW2
#define POW2(val) ((val)*(val))
#endif

typedef float ewaldgpu_real;

void addEwaldgpuForce(double r_cut, int num_kx, int num_ky, int num_kz, double alpha);

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

class EwaldgpuForce : public Actor {
public:
	EwaldgpuForce(SystemInterface &s, double r_cut, int num_kx, int num_ky, int num_kz, double alpha);
	~EwaldgpuForce();
	void setup(SystemInterface &s);
	void cuda_check_error(const dim3 &block, const dim3 &grid,char *function, char *file, unsigned int line);
	void computeForces(SystemInterface &s);
	void computeEnergy(SystemInterface &s);
	//Set parameters
	int set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha);
	int set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps);
	//TUNING r_cut, num_kx, num_ky, num_kz, alpha
	int adaptive_tune(char **log);
	double error_estimate_r(double q_sqr, int N, double r_cut, double V, double alpha, double accuracy);
	double error_estimate_k(double q_sqr, int N, int K, double V, double alpha, double accuracy);
	double tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N);
	double tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N);
	int determine_calc_time_steps();
	//Kolaffa compute optimal alpha
	double compute_E_error_estimate_r(double alpha, double rcut, double q_sqr, double box_l[3]);
	double compute_E_error_estimate_k(double alpha, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3]);
	double E_estimate_error(double rcut, int num_kx, int num_ky, int num_kz, double alpha, double q_sqr, double box_l[3]);
	double compute_optimal_alpha(double rcut, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3], double precision);
	double compute_q_sqare();

protected:
	//SYSTEM
  double m_box_l[3];  // box length
  double m_V; // Volume
  double m_coulomb_prefactor; // Forces, energy is multiplied with this factor
  int   m_N; // number of particles
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
	void GPU_Forces(SystemInterface &s); // Run GPU forces part
	void GPU_Energy(); // Run GPU energy part
	int nextPow2(int x); // Determine the next power of x
	bool isPow2(int x); // Determine if x of power 2
	void getNumBlocksAndThreads(int Size, int maxBlocks, int maxThreads, int &blocks, int &threads); // Determine the number of blocks and threads in GPU part
	//OUTPUT
	void Output(); // Output in terminal
	//REAL SPACE
	void EwaldCPU_EnergySelf(); // Calculate the self energy
};

extern EwaldgpuForce *ewaldgpuForce;

void addEwaldgpuForce(double rcut, int num_kx, int num_ky, int num_kz, double alpha);

#endif
#endif

