#include "config.hpp"

#ifdef MMM1D_GPU

#include "SystemInterface.hpp"
#include "Actor.hpp"
#include <iostream>

typedef float mmm1dgpu_real;

void addMmm1dgpuForce(double maxPWerror, double far_switch_radius, int bessel_cutoff);

class Mmm1dgpuForce : public Actor
{
public:
	// constructor
	Mmm1dgpuForce(SystemInterface &s, mmm1dgpu_real coulomb_prefactor, mmm1dgpu_real maxPWerror, mmm1dgpu_real far_switch_radius = -1, int bessel_cutoff = -1);
	~Mmm1dgpuForce();
	// interface methods
	void computeForces(SystemInterface &s);
	void computeEnergy(SystemInterface &s);
	// configuration methods
	void setup(SystemInterface &s);
	void tune(SystemInterface &s, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff);
	void set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff);

private:
	// CUDA parameters
	unsigned int numThreads;
	unsigned int numBlocks(SystemInterface &s);

	// the box length currently set on the GPU
	// Needed to make sure it hasn't been modified using setmd after inter coulomb was used.
	mmm1dgpu_real host_boxz;
	// the number of particles we had during the last run. Needed to check if we have to realloc dev_forcePairs
	int host_npart;
	bool need_tune;

	// pairs==0: return forces using atomicAdd
	// pairs==1: return force pairs
	// pairs==2: return forces using a global memory reduction
	int pairs;
	// variables for forces and energies calculated pre-reduction
	mmm1dgpu_real *dev_forcePairs, *dev_energyBlocks;

	// MMM1D parameters
	mmm1dgpu_real coulomb_prefactor, maxPWerror, far_switch_radius;
	int bessel_cutoff;

	// run a single force calculation and return the time it takes using high-precision CUDA timers
	float force_benchmark(SystemInterface &s);
};

#endif
