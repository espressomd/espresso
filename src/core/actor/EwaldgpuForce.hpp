#ifndef EWALDGPUFORCE_HPP
#define EWALDGPUFORCE_HPP
#include "config.hpp"

#ifdef EWALD_GPU

#include "SystemInterface.hpp"
#include "Actor.hpp"
#include "particle_data.hpp"
#include <math.h>

typedef float ewaldgpu_real;

void addEwaldgpuForce(double r_cut, int num_kx, int num_ky, int num_kz, double alpha);

typedef struct {
	//Ewald parameters
	double rcut;
	int num_kx;
	int num_ky;
	int num_kz;
	double alpha;
	//Tuning
	double accuracy;
	double precision;
	bool isTuned; // Tuning is over
	bool isTunedFlag; // Flag tuning is over
	int K_max; // Maximal reciprocal K-vector in tuning
	int time_calc_steps; // Steps in time_force_calc function
} Ewaldgpu_params;

extern Ewaldgpu_params ewaldgpu_params;

class EwaldgpuForce : public Actor {
public:
  EwaldgpuForce(SystemInterface &s, double r_cut, int num_kx, int num_ky, int num_kz, double alpha);
	~EwaldgpuForce();
	void setup(SystemInterface &s);
	void computeForces(SystemInterface &s);
	void computeEnergy(SystemInterface &s);
	//Set parameters
	int set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha);
	int set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps);
	//Tuning
	int adaptive_tune(char **log, SystemInterface &s);
	double error_estimate_r(double q_sqr, int N, double r_cut, double V, double alpha, double accuracy) const ;
	double error_estimate_k(double q_sqr, int N, int K, double V, double alpha, double accuracy) const;
	double tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N);
	double tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N);
	int determine_calc_time_steps();
	//Kolaffa compute optimal alpha
	double compute_E_error_estimate_r(double alpha, double rcut, double q_sqr, double box_l[3]);
	double compute_E_error_estimate_k(double alpha, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3]);
	double E_estimate_error(double rcut, int num_kx, int num_ky, int num_kz, double alpha, double q_sqr, double box_l[3]);
	double compute_optimal_alpha(double rcut, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3], double precision);
	double compute_q_sqare(Particle *particle);

protected:
	//System
  double m_box_l[3];  //Box length
  double m_V; //Volume
  double m_coulomb_prefactor; //Forces, energy is multiplied with this factor
  int   m_N; //Number of particles
  ewaldgpu_real *m_k; //k-vectors
  ewaldgpu_real *m_dev_k; //k-vectors
  ewaldgpu_real *m_forces_reci; //Forces k-space
  ewaldgpu_real *m_dev_forces_reci; //Forces k-space
  ewaldgpu_real *m_infl_factor; //Influence factor
  ewaldgpu_real *m_dev_infl_factor; //Influence factor
  ewaldgpu_real *m_dev_rho_hat; //Rho hat
  ewaldgpu_real *m_energy_reci; //Energy k-space
  ewaldgpu_real *m_dev_energy_reci; //Energy k-space
  ewaldgpu_real *m_q_sqr; //Sum of the squares of the particle charges
  ewaldgpu_real *m_dev_q_sqr; //Sum of the squares of the particle charges
  ewaldgpu_real m_energy_self; //Energy self
  ewaldgpu_real m_energy_tot; //Energy total
	//Parameters
  double m_alpha; //Separation parameter
  double m_rcut; //Cutoff radius
  int m_num_kx, m_num_ky, m_num_kz, m_num_k; //Number of k's in xyz-direction

protected:
  //Variables
  bool m_isTuned; //Tuning is over
  SystemInterface *System;
  //Ccompute k's
  void compute_num_k(); //Compute the number of k's in sphere respectively in ellipsoid
  void compute_k_AND_influence_factor(); //Compute the k-vectors in sphere respectively in ellipsoid and influence factor
  //GPU program
  void GPU_Forces(SystemInterface &s); //Run GPU forces part
  void GPU_Energy(SystemInterface &s); //Run GPU energy part
  void GPU_q_sqr(SystemInterface &s); // Compute q squared
  int nextPow2(int x); //Determine the next power of x
  bool isPow2(int x); //Determine if x of power 2
  void getNumBlocksAndThreads(int Size, int maxBlocks, int maxThreads, int &blocks, int &threads); // Determine the number of blocks and threads in GPU part
  //Output
  void Output(); //Output in terminal
  //Real space
  void EwaldCPU_EnergySelf(); //Calculate the self energy
  size_t shared_mem_per_block;
};

extern EwaldgpuForce *ewaldgpuForce;

void addEwaldgpuForce(double rcut, int num_kx, int num_ky, int num_kz, double alpha);

#endif
#endif

