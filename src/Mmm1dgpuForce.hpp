#include "config.hpp"

#ifdef MMM1D_GPU

#include "SystemInterface.hpp"
#include <iostream>

typedef float mmm1dgpu_real;

class Mmm1dgpuForce
{
public:
	// constructor
	Mmm1dgpuForce(SystemInterface &s, mmm1dgpu_real coulomb_prefactor, mmm1dgpu_real maxPWerror, mmm1dgpu_real far_switch_radius = -1, int bessel_cutoff = -1);
	~Mmm1dgpuForce();
	// interface methods
	void calc(SystemInterface &s);
	void calc_energy(SystemInterface &s);
	// configuration methods
	void tune(SystemInterface &s, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff);
	void set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff);

private:
	// the box length currently set on the GPU
	// Needed to make sure it hasn't been modified using setmd after inter coulomb was used.
	mmm1dgpu_real host_boxz;

	// pairs==0: return forces using atomicAdd
	// pairs==1: return force pairs
	// pairs==2: return forces using a global memory reduction
	int pairs;

	mmm1dgpu_real coulomb_prefactor, maxPWerror, far_switch_radius;
	int bessel_cutoff;

	// run a single force calculation and return the time it takes using high-precision CUDA timers
	float force_benchmark(SystemInterface &s);
};

extern Mmm1dgpuForce *mmm1dgpuForce;

#endif
