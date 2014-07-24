#include "actor/Mmm1dgpuForce.hpp"
#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "interaction_data.hpp"
#include "forces.hpp"
#include "grid.hpp"

#ifdef MMM1D_GPU

void addMmm1dgpuForce(double maxPWerror, double switch_rad, int bessel_cutoff)
{
	/* coulomb.prefactor apparently is not available yet at this point */
	static Mmm1dgpuForce *mmm1dgpuForce = NULL;
	if (!mmm1dgpuForce) // inter coulomb mmm1dgpu was never called before
	{
		mmm1dgpuForce = new Mmm1dgpuForce(espressoSystemInterface, 0, maxPWerror, switch_rad, bessel_cutoff);
		potentials.push_back(mmm1dgpuForce);
	}
	/* set_params needs to be called both upon create and upon update because it is responsible for writing
		to the struct from which the TCL command "inter coulomb" retrieves the current parameter set */
	mmm1dgpuForce->set_params(0, 0, maxPWerror, switch_rad, bessel_cutoff, true);

	// turn on MMM1DGPU
	coulomb.method = COULOMB_MMM1D_GPU;
	mpi_bcast_coulomb_params();
}

void Mmm1dgpuForce::disable()
{
	if (coulomb.method == COULOMB_MMM1D_GPU)
	{
		coulomb.method = COULOMB_NONE;
		mpi_bcast_coulomb_params();
	}
}

void Mmm1dgpuForce::check_periodicity()
{
	if (PERIODIC(0) || PERIODIC(1) || !PERIODIC(2))
	{
		std::cerr << "MMM1D requires periodicity (0,0,1)" << std::endl;
		exit(EXIT_FAILURE);
	}
}

#endif
