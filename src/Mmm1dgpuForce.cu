#include "Mmm1dgpuForce.hpp"
#include "cuda_utils.hpp"

#ifdef MMM1D_GPU

Mmm1dgpuForce *mmm1dgpuForce = 0;

// the code is mostly multi-GPU capable, but Espresso is not yet
const int deviceCount = 1;
float multigpu_factors[] = {1.0};
#define cudaSetDevice(d)
#define HANDLE_ERROR(a) cuda_safe_mem(a) // TODO: inline

#include "mmm-common_cuda.hpp"
#include "mmm1d.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"

const mmm1dgpu_real C_GAMMAf = C_GAMMA;
const mmm1dgpu_real C_2PIf = C_2PI;

__constant__ mmm1dgpu_real far_switch_radius_2 = 0.05*0.05;
__constant__ mmm1dgpu_real boxz;
__constant__ mmm1dgpu_real uz;
__constant__ mmm1dgpu_real coulomb_prefactor = 1.0;
__constant__ int bessel_cutoff = 5;
__constant__ mmm1dgpu_real maxPWerror = 1e-5;

Mmm1dgpuForce::Mmm1dgpuForce(SystemInterface &s, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror,
	mmm1dgpu_real _far_switch_radius, int _bessel_cutoff)
:coulomb_prefactor(_coulomb_prefactor), maxPWerror(_maxPWerror), far_switch_radius(_far_switch_radius), bessel_cutoff(_bessel_cutoff)
{
	// interface sanity checks
	if(!s.requestFGpu())
		std::cerr << "Mmm1dgpuForce needs access to forces on GPU!" << std::endl;

	if(!s.requestRGpu())
		std::cerr << "Mmm1dgpuForce needs access to positions on GPU!" << std::endl;

	// system sanity checks
	if (PERIODIC(0) || PERIODIC(1) || !PERIODIC(2))
	{
		printf("MMM1D requires periodicity (0,0,1)\n");
		exit(EXIT_FAILURE);
	}

	if (s.box()[2] <= 0)
	{
		printf("Error: Please set box length before initializing MMM1D!\n");
		exit(EXIT_FAILURE);
	}

	if (s.npart() <= 0 && far_switch_radius < 0)
	{
		std::cerr << "Warning: Please add particles to system before intializing MMM1D." << std::endl;
		std::cerr << "Tuning will not be performed! Setting far switch radius to half box length." << std::endl;
		far_switch_radius = s.box()[2]/2;
	}

	// turn on MMM1DGPU
	modpsi_init();
	set_params(s.box()[2], coulomb_prefactor, maxPWerror, far_switch_radius, bessel_cutoff);
	coulomb.method = COULOMB_MMM1D_GPU;
	mpi_bcast_coulomb_params();
	tune(s, maxPWerror, far_switch_radius, bessel_cutoff);

	// for all but the largest systems, it is faster to store force pairs and then sum them up
	// so unless we're limited by memory, do the latter
	pairs = 2;
	for (int d = 0; d < deviceCount; d++)
	{
		cudaSetDevice(d);

		size_t freeMem, totalMem;
		cudaMemGetInfo(&freeMem, &totalMem);
		if (freeMem/2 < s.npart()*s.npart()*sizeof(mmm1dgpu_real)) // don't use more than half the device's memory
		{
			std::cerr << "Switching to atomicAdd due to memory constraints." << std::endl;
			pairs = 0;
			break;
		}
	}
}

Mmm1dgpuForce::~Mmm1dgpuForce()
{
	modpsi_destroy();

	// TODO: unset coulomb.method
}

__forceinline__ __device__ mmm1dgpu_real sqpow(mmm1dgpu_real x)
{
	return pow(x,2);
}
__forceinline__ __device__ mmm1dgpu_real cbpow(mmm1dgpu_real x)
{
	return pow(x,3);
}

__device__ void sumReduction(mmm1dgpu_real *input, mmm1dgpu_real *sum)
{
	int tid = threadIdx.x;
	for (int i = blockDim.x/2; i > 0; i /= 2)
	{
		__syncthreads();
		if (tid < i)
			input[tid] += input[i+tid];
	}
	__syncthreads();
	if (tid == 0)
		sum[0] = input[0];
}

__global__ void sumKernel(mmm1dgpu_real *data, int N)
{
	extern __shared__ mmm1dgpu_real partialsums[];
	if (blockIdx.x != 0) return;
	int tid = threadIdx.x;
	mmm1dgpu_real result = 0;
	
	for (int i = 0; i < N; i += blockDim.x)
	{
		if (i+tid >= N)
			partialsums[tid] = 0;
		else
			partialsums[tid] = data[i+tid];
		
		sumReduction(partialsums, &result);
		if (tid == 0)
		{
			if (i == 0) data[0] = 0;
			data[0] += result;
		}
	}
}

__global__ void besselTuneKernel(int *result, mmm1dgpu_real far_switch_radius, int maxCut)
{
	mmm1dgpu_real arg = C_2PIf*uz*far_switch_radius;
	mmm1dgpu_real pref = 4*uz*max(1.0f, C_2PIf*uz);
	mmm1dgpu_real err;
	int P = 1;
	do
	{
		err = pref*dev_K1(arg*P)*exp(arg)/arg*(P-1 + 1/arg);
		P++;
	} while (err > maxPWerror && P <= maxCut);
	P--;

	result[0] = P;
}

void Mmm1dgpuForce::tune(SystemInterface &s, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff)
{
	mmm1dgpu_real far_switch_radius = _far_switch_radius;
	int bessel_cutoff = _bessel_cutoff;
	mmm1dgpu_real maxrad = host_boxz;

	if (_far_switch_radius < 0 && _bessel_cutoff < 0)
	// autodetermine switching and bessel cutoff radius
	{
		mmm1dgpu_real bestrad = 0, besttime = INFINITY;

		for (far_switch_radius = 0.05*maxrad; far_switch_radius < maxrad; far_switch_radius += 0.05*maxrad)
		{
			set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
			tune(s, _maxPWerror, far_switch_radius, -2); // tune bessel cutoff
			int runtime = force_benchmark(s);
			if (runtime < besttime)
			{
				besttime = runtime;
				bestrad = far_switch_radius;
			}
		}
		far_switch_radius = bestrad;

		set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
		tune(s, _maxPWerror, far_switch_radius, -2); // tune bessel cutoff
	}

	else if (_bessel_cutoff < 0)
	// autodetermine bessel cutoff
	{
		int *dev_cutoff;
		int maxCut = 30;
		HANDLE_ERROR( cudaMalloc((void**)&dev_cutoff, sizeof(int)) );
		besselTuneKernel<<<1,1>>>(dev_cutoff, far_switch_radius, maxCut);
		HANDLE_ERROR( cudaMemcpy(&bessel_cutoff, dev_cutoff, sizeof(int), cudaMemcpyDeviceToHost) );
		cudaFree(dev_cutoff);
		if (_bessel_cutoff != -2 && bessel_cutoff >= maxCut) // we already have our switching radius and only need to determine the cutoff, i.e. this is the final tuning round
		{
			printf("No reasonable Bessel cutoff could be determined.\n");
			exit(EXIT_FAILURE);
		}

		set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
	}
}

void Mmm1dgpuForce::set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff)
{
	if (_boxz > 0 && _far_switch_radius > _boxz)
	{
		printf("Far switch radius (%f) must not be larger than the box length (%f).\n", _far_switch_radius, _boxz);
		exit(EXIT_FAILURE);
	}
	mmm1dgpu_real _far_switch_radius_2 = _far_switch_radius*_far_switch_radius;
	mmm1dgpu_real _uz = 1.0/_boxz;
	for (int d = 0; d < deviceCount; d++)
	{
		// double colons are needed to access the constant memory variables because they
		// are file globals and we have identically named class variables
		cudaSetDevice(d);
		if (_far_switch_radius >= 0)
		{
			HANDLE_ERROR( cudaMemcpyToSymbol(::far_switch_radius_2, &_far_switch_radius_2, sizeof(mmm1dgpu_real)) );
			mmm1d_params.far_switch_radius_2 = _far_switch_radius*_far_switch_radius;
		}
		if (_boxz > 0)
		{
			host_boxz = _boxz;
			HANDLE_ERROR( cudaMemcpyToSymbol(::boxz, &_boxz, sizeof(mmm1dgpu_real)) );
			HANDLE_ERROR( cudaMemcpyToSymbol(::uz, &_uz, sizeof(mmm1dgpu_real)) );
		}
		if (_coulomb_prefactor != 0)
		{
			HANDLE_ERROR( cudaMemcpyToSymbol(::coulomb_prefactor, &_coulomb_prefactor, sizeof(mmm1dgpu_real)) );
		}
		if (_bessel_cutoff > 0)
		{
			HANDLE_ERROR( cudaMemcpyToSymbol(::bessel_cutoff, &_bessel_cutoff, sizeof(int)) );
			mmm1d_params.bessel_cutoff = _bessel_cutoff;
		}
		if (_maxPWerror > 0)
		{
			HANDLE_ERROR( cudaMemcpyToSymbol(::maxPWerror, &_maxPWerror, sizeof(mmm1dgpu_real)) );
			mmm1d_params.maxPWerror = _maxPWerror;
		}
	}

	if (_far_switch_radius >= 0 && _bessel_cutoff > 0)
		printf("@@@ Using far switch radius %f and bessel cutoff %d\n", _far_switch_radius, _bessel_cutoff);
}

void Mmm1dgpuForce::calc(SystemInterface &s)
{
	
}

void Mmm1dgpuForce::calc_energy(SystemInterface &s)
{
	
}

float Mmm1dgpuForce::force_benchmark(SystemInterface &s)
{
	cudaEvent_t eventStart, eventStop;
	cudaStream_t stream;
	float elapsedTime;

	cudaStreamCreate(&stream);
	HANDLE_ERROR( cudaEventCreate(&eventStart) );
	HANDLE_ERROR( cudaEventCreate(&eventStop) );
	HANDLE_ERROR( cudaEventRecord(eventStart, stream) );
	Mmm1dgpuForce::calc(s);
	HANDLE_ERROR( cudaEventRecord(eventStop, stream) );
	HANDLE_ERROR( cudaEventSynchronize(eventStop) );
	HANDLE_ERROR( cudaEventElapsedTime(&elapsedTime, eventStart, eventStop) );
	printf(">>> Calculated in %3.3f ms\n", elapsedTime);
	HANDLE_ERROR( cudaEventDestroy(eventStart) );
	HANDLE_ERROR( cudaEventDestroy(eventStop) );

	return elapsedTime;
}

#endif /* MMM1D_GPU */
