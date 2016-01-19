#include "actor/EwaldGPU.hpp"

#ifdef EWALD_GPU

#include <iostream>

#include "tuning.hpp"
#include "integrate.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"
#include "cells.hpp"

#define DEBUG(A) std::cout << #A << ": " << A << std::endl;

typedef ewaldgpu_real real;
extern Ewaldgpu_params ewaldgpu_params;

//Compute reciprocal space vectors k
void EwaldgpuForce::compute_num_k()
{
  int Index=0;//Arrayindex

  for (long long ix=0; ix<=m_num_kx; ix++)
    {
      for (long long iy=-m_num_ky; iy<=m_num_ky; iy++)
	{
	  for (long long iz=-m_num_kz; iz<=m_num_kz; iz++)
	    {
	      if(ix*(2*m_num_ky+1)*(2*m_num_kz+1)+iy*(2*m_num_kz+1)+iz >= 0//Half m_k-space
		 and SQR(ix)*SQR(m_num_ky)*SQR(m_num_kz)
		 + SQR(iy)*SQR(m_num_kx)*SQR(m_num_kz)
		 + SQR(iz)*SQR(m_num_kx)*SQR(m_num_ky)
		 <=SQR(m_num_kx)*SQR(m_num_ky)*SQR(m_num_kz))//m_k-space ellipsoid
		{
		  Index+=1;
		}
	    }
	}
    }
  m_num_k=Index;
}
void EwaldgpuForce::compute_k_AND_influence_factor()
{
  real k_sqr;//Absolute square of m_k-vector
  int Index=0;//Arrayindex

	for (long long ix=0; ix<=m_num_kx; ix++)//Half m_k-space
	{
		for (long long iy=-m_num_ky; iy<=m_num_ky; iy++)
		{
			for (long long iz=-m_num_kz; iz<=m_num_kz; iz++)
			{
				if(ix*(2*m_num_ky+1)*(2*m_num_kz+1)+iy*(2*m_num_kz+1)+iz >= 0//Half m_k-space
				   and SQR(ix)*SQR(m_num_ky)*SQR(m_num_kz)
				   	 + SQR(iy)*SQR(m_num_kx)*SQR(m_num_kz)
				   	 + SQR(iz)*SQR(m_num_kx)*SQR(m_num_ky)
				   	 <=SQR(m_num_kx)*SQR(m_num_ky)*SQR(m_num_kz))//m_k-space ellipsoid
				{
					//m_k vectors
					m_k[Index] = 2*M_PI/m_box_l[0]*ix;
					m_k[Index+m_num_k] = 2*M_PI/m_box_l[1]*iy;
					m_k[Index+2*m_num_k] = 2*M_PI/m_box_l[2]*iz;
					//Influence factor
					k_sqr= SQR(m_k[Index]) + SQR(m_k[Index+m_num_k]) + SQR(m_k[Index+2*m_num_k]);
					m_infl_factor[Index] = 8*M_PI/m_V*expf(-k_sqr/(4*SQR(m_alpha)))/k_sqr;
					//Index
					Index+=1;
				}
			}
		}
	}	m_infl_factor[0] = 0;//Influence factor at m_k=(0,0,0)
}

//Real space
void EwaldgpuForce::EwaldCPU_EnergySelf()
{
	m_energy_self=coulomb.prefactor * (-m_alpha/sqrt(M_PI) * m_q_sqr[0]);
}

//Parameters
int EwaldgpuForce::set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha)
{
	ewaldgpu_params.rcut = rcut;
	ewaldgpu_params.num_kx = num_kx;
	ewaldgpu_params.num_ky = num_ky;
	ewaldgpu_params.num_kz = num_kz;
	ewaldgpu_params.alpha = alpha;

	return 0;
}
int EwaldgpuForce::set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps)
{
	ewaldgpu_params.accuracy = accuracy;
	ewaldgpu_params.precision = precision;
	ewaldgpu_params.K_max = K_max;
	ewaldgpu_params.time_calc_steps = time_calc_steps;

	return 0;
}

//  Tuning
int EwaldgpuForce::adaptive_tune(char **log, SystemInterface &s)
{
  ewaldgpu_params.isTuned = false;
  int Kmax = ewaldgpu_params.K_max;
  double alpha_array[Kmax]; //  All computed alpha in dependence of K
  double rcut_array[Kmax]; //  All computed r_cut in dependence of all computed alpha

  //Squared charge
  Particle *particle;
  particle = (Particle*)Utils::malloc(n_part*sizeof(Particle));
  mpi_get_particles(particle, NULL);
  double q_sqr = compute_q_sqare(particle);

  char b[3*ES_INTEGER_SPACE + 3*ES_DOUBLE_SPACE + 128];

  if (skin == -1) {
    *log = strcat_alloc(*log, "ewaldgpu cannot be tuned, since the skin is not yet set");
    return ES_ERROR;
  }

  //Compute alpha for all reciprocal k-sphere radius K
  for(int K = 0; K < Kmax ;K++)
  {
    alpha_array[K] = tune_alpha(ewaldgpu_params.accuracy/sqrt(2), ewaldgpu_params.precision, K+1, box_l[0]*box_l[1]*box_l[2], q_sqr, n_part);
  }
  //Compute r_cut for all computed alpha
  for(int K = 0; K < Kmax ;K++)
  {
    rcut_array[K] = tune_rcut(ewaldgpu_params.accuracy/sqrt(2), ewaldgpu_params.precision, alpha_array[K], box_l[0]*box_l[1]*box_l[2], q_sqr, n_part);
  }
  //Test if accuracy was reached
  if(rcut_array[Kmax-1]<0)
  {
    return ES_ERROR;
  }

  /***********************************************************************************
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	PERFORMANCE TIME
  ***********************************************************************************/

  //Test performance time for the diverent (K, rcut, alpha)
  double int_time_best = std::numeric_limits<double>::max();
  int K_best = Kmax;
  
  for(int K = 0; (K < Kmax); K++)
  {
    if(alpha_array[K]>0 and rcut_array[K]>0 and rcut_array[K]<(std::min(box_l[0],std::min(box_l[1],box_l[2])))/2.0-skin)
    {
      set_params(rcut_array[K], K+1, K+1, K+1, alpha_array[K]);
      mpi_bcast_coulomb_params();
      const double int_time = time_force_calc(ewaldgpu_params.time_calc_steps);
      if(int_time<int_time_best)
      {
        int_time_best = int_time;
        K_best = K;
      } else {
        break;
      }
    }                
  }

  set_params(rcut_array[K_best], K_best+1, K_best+1, K_best+1, alpha_array[K_best]);
  ewaldgpu_params.isTuned = true;
  mpi_bcast_coulomb_params();

  //Print Status
  sprintf(b, "ewaldgpu tune parameters: Accuracy goal = %f\n", ewaldgpu_params.accuracy);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: Alpha = %f\n", ewaldgpu_params.alpha);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: r_cut = %f\n", ewaldgpu_params.rcut);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_kx = %i\n", ewaldgpu_params.num_kx);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_ky = %i\n", ewaldgpu_params.num_ky);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_kz = %i\n", ewaldgpu_params.num_kz);
  *log = strcat_alloc(*log, b);

  return ES_OK;
}

inline double EwaldgpuForce::error_estimate_r(double q_sqr, int N, double r_cut, double V, double alpha, double accuracy) const
{
	return 2 * q_sqr * exp(-alpha*alpha*r_cut*r_cut) / sqrt(N*V*r_cut) - accuracy ;  // Kolafa-Perram, eq. 18
}

inline double EwaldgpuForce::error_estimate_k(double q_sqr, int N, int K, double V, double alpha, double accuracy) const
{
	return sqrt(q_sqr/N) * alpha / (pow(V,1/3.0) * M_PI) * sqrt(8*q_sqr / K) * exp(-M_PI*M_PI * K*K / (alpha*alpha * pow(V,2/3.0))) - accuracy;  // Kolafa-Perram, eq. 32
}

double EwaldgpuForce::tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N)
{
  double alpha_low=0.01;
  double alpha_high=100;
  double alpha_guess;
  double fkt_low;
  double fkt_high;
  double fkt_guess;

  // Find alpha with given K in k-space error estimate via bisection
  fkt_low = error_estimate_k(q_sqr, N, K, V, alpha_low, accuracy);
  fkt_high = error_estimate_k(q_sqr, N, K, V, alpha_high, accuracy);

  if (fkt_low*fkt_high > 0.0)
    {
      return -1; // Value unusable
    }

  do
    {
      alpha_guess = 0.5 *(alpha_low + alpha_high);
      fkt_guess = error_estimate_k(q_sqr, N, K, V, alpha_guess, accuracy);

      if (fkt_low*fkt_guess < 0.0) 
	alpha_high = alpha_guess;
      else 
	alpha_low = alpha_guess;

    } while (fabs(alpha_low-alpha_high) > precision);

  return 0.5 *(alpha_low + alpha_high);
}

double EwaldgpuForce::tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N)
{
  double rcut_low=std::numeric_limits<double>::epsilon();
  /* Limit maximal cutoff to 10000 particles in rcut sphere. */
  
  double rcut_high=std::min(std::min(box_l[0],std::min(box_l[1],box_l[2])), pow((3000.0*V)/(4.0*3.14159*N), 1.0/3.0));
  double rcut_guess;
  double fkt_low;
  double fkt_high;
  double fkt_guess;
  
  // Find rcut with given K in k-space error estimate via bisection
  fkt_low = error_estimate_r(q_sqr,  N, rcut_low, V, alpha, accuracy);
  fkt_high = error_estimate_r(q_sqr, N, rcut_high, V, alpha, accuracy);
  
  if (fkt_low*fkt_high > 0.0)
  {
    return -1; // Value unusable
  }

  do
  {
    rcut_guess = 0.5 *(rcut_low + rcut_high);
    fkt_guess = error_estimate_r(q_sqr, N, rcut_guess, V, alpha, accuracy);
    if (fkt_low*fkt_guess < 0.0) rcut_high = rcut_guess;
    else rcut_low = rcut_guess;

  } while (fabs(rcut_low-rcut_high) > precision);

  return 0.5 *(rcut_low + rcut_high);
}

int EwaldgpuForce::determine_calc_time_steps()
{
  Cell *cell;
  Particle *part;
  int i,c,np;
  int sum_qpart =0;

  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++)
    {
      if( part[i].p.q != 0.0 )
      {
      	sum_qpart += 1.0;
      }
    }
  }
  sum_qpart    = (int)(sum_qpart+0.1);
  return (1999 + sum_qpart)/sum_qpart;
}

//Kolaffa compute optimal alpha
double EwaldgpuForce::compute_E_error_estimate_r(double alpha, double rcut, double q_sqr, double box_l[3])
{
	//Compute the r space part of the force error estimate
  double std_deviation;
  std_deviation = q_sqr*pow(rcut/(2.0*box_l[0]*box_l[1]*box_l[2]),0.5) * exp(-SQR(alpha)*SQR(rcut)) / (SQR(alpha*rcut));  // Kolafa-Perram, eq. 18

  return std_deviation;
}
double EwaldgpuForce::compute_E_error_estimate_k(double alpha, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3])
{
  //Compute the r space part of the force error estimate
  double std_deviation;
	std_deviation = q_sqr * alpha * pow(M_PI,-2.0) * pow(num_kx,-1.5) * exp(-(SQR(M_PI*num_kx/(alpha*(box_l[0]+box_l[1]+box_l[2])/3.0)))); //Kolafa Perram, eq. 32
  return std_deviation;
}
double EwaldgpuForce::E_estimate_error(double rcut, int num_kx, int num_ky, int num_kz, double alpha, double q_sqr, double box_l[3])
{
	//Compute the energy_reci error estimate
  return sqrt(pow(compute_E_error_estimate_r(alpha, rcut, q_sqr, box_l),2) + pow(compute_E_error_estimate_k(alpha, num_kx, num_ky, num_kz, q_sqr, box_l),2));
}
double EwaldgpuForce::compute_optimal_alpha(double rcut, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3], double precision)
{
  //Use bisectional method to get optimal alpha value
  double alpha_low, f_low;
  double alpha_high, f_high;
  double alpha_guess, f_guess;
	int counter=0;
  alpha_low = 0.01;
  alpha_high = 10.0;

  f_low = compute_E_error_estimate_r(alpha_low, rcut, q_sqr, box_l) - compute_E_error_estimate_k(alpha_low, num_kx, num_ky, num_kz, q_sqr, box_l);
  f_high = compute_E_error_estimate_r(alpha_high, rcut, q_sqr, box_l) - compute_E_error_estimate_k(alpha_high, num_kx, num_ky, num_kz, q_sqr, box_l);

  if (f_low*f_high > 0.0)
  {
    fprintf(stderr, "Error: Could not init method to find optimal alpha!\n");
    exit(1);
  }

  do
  {
    alpha_guess = 0.5 *(alpha_low + alpha_high);
    f_guess = compute_E_error_estimate_r(alpha_guess, rcut, q_sqr, box_l) - compute_E_error_estimate_k(alpha_guess, num_kx, num_ky, num_kz, q_sqr, box_l);
    if (f_low*f_guess < 0.0) alpha_high = alpha_guess;
    else alpha_low = alpha_guess;
    counter++;
		if(counter>10000)
		{
			fprintf(stderr, "Find optimal alpha: Maximum number of loops exceeded: %i loops",counter);
			exit(1);
		}
  } while (fabs(alpha_low-alpha_high) > precision);

  return 0.5 *(alpha_low + alpha_high);
}
double EwaldgpuForce::compute_q_sqare(Particle *particle)
{
	double q_sqr=0;

	for(int i = 0; i < n_part; i++)
	{
		q_sqr += SQR(particle[i].p.q);
	}
  return q_sqr;
}

#endif
