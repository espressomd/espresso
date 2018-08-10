#include "actor/Mmm1dgpuForce.hpp"
#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "forces.hpp"
#include "grid.hpp"

#ifdef MMM1D_GPU

void Mmm1dgpuForce::check_periodicity()
{
	if (PERIODIC(0) || PERIODIC(1) || !PERIODIC(2))
	{
		std::cerr << "MMM1D requires periodicity (0,0,1)" << std::endl;
		exit(EXIT_FAILURE);
	}
}

#endif
