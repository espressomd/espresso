#include "mmm-common.hpp"
#include "specfunc_cuda.hpp"

// order hardcoded. mmm1d_recalcTables() typically does order less than 30.
// As the coefficients are stored in __constant__ memory, the array needs to be sized in advance.
// We don't know exactly how many coefficients per order, so we size plentiful.
const int modpsi_order = 30;
const int modpsi_constant_size = modpsi_order*modpsi_order*2;

// linearized array on host
int *linModPsi_offsets = NULL, *linModPsi_lengths = NULL;
mmm1dgpu_real *linModPsi = NULL;

// linearized array on device
__constant__ int device_n_modPsi = 0;
__constant__ int device_linModPsi_offsets[2*modpsi_order], device_linModPsi_lengths[2*modpsi_order];
__constant__ mmm1dgpu_real device_linModPsi[modpsi_constant_size];

int modpsi_init()
{
	if (n_modPsi < modpsi_order)
	{
		create_mod_psi_up_to(modpsi_order);
	}
	
	// linearize the coefficients array
	linModPsi_offsets = (int*) realloc(linModPsi_offsets, sizeof(int) * 2*n_modPsi);
	linModPsi_lengths = (int*) realloc(linModPsi_lengths, sizeof(int) * 2*n_modPsi);
	for (int i = 0; i < 2*n_modPsi; i++)
	{
		if (i == 0)
			linModPsi_offsets[i] = 0;
		else
			linModPsi_offsets[i] = linModPsi_offsets[i-1] + linModPsi_lengths[i-1];
		linModPsi_lengths[i] = modPsi[i].n;
	}
	linModPsi = (mmm1dgpu_real*) realloc(linModPsi, sizeof(mmm1dgpu_real) * (linModPsi_offsets[2*n_modPsi-1] + linModPsi_lengths[2*n_modPsi-1]));
	for (int i = 0; i < 2*n_modPsi; i++)
	{
		for (int j = 0; j < modPsi[i].n; j++)
		{
			linModPsi[linModPsi_offsets[i] + j] = (mmm1dgpu_real) modPsi[i].e[j]; // cast to single-precision if necessary
		}
	}

	for (int d = 0; d < deviceCount; d++)
	{
		cudaSetDevice(d);
		
		// copy to GPU
		int linModPsiSize = linModPsi_offsets[2*n_modPsi-1] + linModPsi_lengths[2*n_modPsi-1];
		if (linModPsiSize > modpsi_constant_size)
		{
			printf("ERROR: __constant__ device_linModPsi[] is not large enough\n");
			exit(EXIT_FAILURE);
		}
		HANDLE_ERROR( cudaMemcpyToSymbol(device_linModPsi_offsets, linModPsi_offsets, 2*n_modPsi*sizeof(int)) );
		HANDLE_ERROR( cudaMemcpyToSymbol(device_linModPsi_lengths, linModPsi_lengths, 2*n_modPsi*sizeof(int)) );
		HANDLE_ERROR( cudaMemcpyToSymbol(device_linModPsi, linModPsi, linModPsiSize*sizeof(mmm1dgpu_real)) );
		HANDLE_ERROR( cudaMemcpyToSymbol(device_n_modPsi, &n_modPsi, sizeof(int)) );
	}

	return 0;
}

int modpsi_destroy()
{
	// no need to delete the arrays off the device, they're in constant memory
	// free arrays on host
	free(linModPsi_offsets);
	free(linModPsi_lengths);
	free(linModPsi);
	linModPsi_offsets = NULL;
	linModPsi_lengths = NULL;
	linModPsi = NULL;
	return 0;
}

__device__ mmm1dgpu_real dev_mod_psi_even(int n, mmm1dgpu_real x)
{
	return evaluateAsTaylorSeriesAt(&device_linModPsi[device_linModPsi_offsets[2*n]],
		device_linModPsi_lengths[2*n], x*x);
}

__device__ mmm1dgpu_real dev_mod_psi_odd(int n, mmm1dgpu_real x)
{
	return x*evaluateAsTaylorSeriesAt(&device_linModPsi[device_linModPsi_offsets[2*n+1]],
		device_linModPsi_lengths[2*n+1], x*x);
}
