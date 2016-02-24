/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#include "actor/Mmm1dgpuForce.hpp"
#include "cuda_utils.hpp"

#ifdef MMM1D_GPU

// the code is mostly multi-GPU capable, but Espresso is not yet
const int deviceCount = 1;
float multigpu_factors[] = {1.0};
#define cudaSetDevice(d)

#include "mmm-common_cuda.hpp"
#include "mmm1d.hpp"
#include "interaction_data.hpp"
#include "EspressoSystemInterface.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

__device__ inline void atomicadd(float* address, float value)
{
#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110
	int oldval, newval, readback;
	oldval = __float_as_int(*address);
	newval = __float_as_int(__int_as_float(oldval) + value);
	while ((readback=atomicCAS((int *)address, oldval, newval)) != oldval)
	{
		oldval = readback;
		newval = __float_as_int(__int_as_float(oldval) + value);
	}
#else
#error atomicAdd needs compute capability 1.1 or higher
#endif
}

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
:coulomb_prefactor(_coulomb_prefactor), maxPWerror(_maxPWerror), far_switch_radius(_far_switch_radius),
	bessel_cutoff(_bessel_cutoff), host_boxz(0), host_npart(0), pairs(-1), dev_forcePairs(NULL), dev_energyBlocks(NULL),
	numThreads(64), need_tune(true)
{
	// interface sanity checks
	if(!s.requestFGpu())
		std::cerr << "Mmm1dgpuForce needs access to forces on GPU!" << std::endl;

	if(!s.requestRGpu())
		std::cerr << "Mmm1dgpuForce needs access to positions on GPU!" << std::endl;

	if(!s.requestQGpu())
		std::cerr << "Mmm1dgpuForce needs access to charges on GPU!" << std::endl;

	// system sanity checks
	check_periodicity();

	modpsi_init();
}

void Mmm1dgpuForce::setup(SystemInterface &s)
{
	if (s.box()[2] <= 0)
	{
		std::cerr << "Error: Please set box length before initializing MMM1D!" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (need_tune == true && s.npart_gpu() > 0)
	{
		set_params(s.box()[2], coulomb.prefactor, maxPWerror, far_switch_radius, bessel_cutoff);
		tune(s, maxPWerror, far_switch_radius, bessel_cutoff);
	}
	if (s.box()[2] != host_boxz)
	{
		set_params(s.box()[2], 0,-1,-1,-1);
	}
	if (s.npart_gpu() == host_npart) // unchanged
	{
		return;
	}

	// For all but the largest systems, it is faster to store force pairs and then sum them up.
	// Atomics are just so slow: so unless we're limited by memory, do the latter.
	pairs = 2;
	for (int d = 0; d < deviceCount; d++)
	{
		cudaSetDevice(d);

		size_t freeMem, totalMem;
		cudaMemGetInfo(&freeMem, &totalMem);
		if (freeMem/2 < 3*s.npart_gpu()*s.npart_gpu()*sizeof(mmm1dgpu_real)) // don't use more than half the device's memory
		{
			std::cerr << "Switching to atomicAdd due to memory constraints." << std::endl;
			pairs = 0;
			break;
		}
	}
	if (dev_forcePairs)
		cudaFree(dev_forcePairs);
	if (pairs) // we need memory to store force pairs
	{
		cuda_safe_mem( cudaMalloc((void**)&dev_forcePairs, 3*s.npart_gpu()*s.npart_gpu()*sizeof(mmm1dgpu_real)) );
	}
	if (dev_energyBlocks)
		cudaFree(dev_energyBlocks);
	cuda_safe_mem( cudaMalloc((void**)&dev_energyBlocks, numBlocks(s)*sizeof(mmm1dgpu_real)) );
	host_npart = s.npart_gpu();
}

unsigned int Mmm1dgpuForce::numBlocks(SystemInterface &s)
{
	int b = s.npart_gpu()*s.npart_gpu()/numThreads+1;
	if (b > 65535)
		b = 65535;
	return b;
}

Mmm1dgpuForce::~Mmm1dgpuForce() {
	modpsi_destroy();
	cudaFree(dev_forcePairs);
}

__forceinline__ __device__ mmm1dgpu_real sqpow(mmm1dgpu_real x)
{
  return x*x;
}
__forceinline__ __device__ mmm1dgpu_real cbpow(mmm1dgpu_real x)
{
  return x*x*x;
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
	// autodetermine switching radius and bessel cutoff
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
		cuda_safe_mem( cudaMalloc((void**)&dev_cutoff, sizeof(int)) );
		besselTuneKernel<<<1,1>>>(dev_cutoff, far_switch_radius, maxCut);
		cuda_safe_mem( cudaMemcpy(&bessel_cutoff, dev_cutoff, sizeof(int), cudaMemcpyDeviceToHost) );
		cudaFree(dev_cutoff);
		if (_bessel_cutoff != -2 && bessel_cutoff >= maxCut) // we already have our switching radius and only need to determine the cutoff, i.e. this is the final tuning round
		{
			std::cerr << "No reasonable Bessel cutoff could be determined." << std::endl;
			exit(EXIT_FAILURE);
		}

		set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
	}
}

void Mmm1dgpuForce::set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff, bool manual)
{
	if (_boxz > 0 && _far_switch_radius > _boxz)
	{
		std::cout << "switching radius must not be larger than box length" << std::endl;
		exit(EXIT_FAILURE);
	}
	mmm1dgpu_real _far_switch_radius_2 = _far_switch_radius*_far_switch_radius;
	mmm1dgpu_real _uz = 1.0/_boxz;
	for (int d = 0; d < deviceCount; d++)
	{
		// double colons are needed to access the constant memory variables because they
		// are file globals and we have identically named class variables
		cudaSetDevice(d);
		if (manual) // tuning needs to be performed again
		{
			far_switch_radius = _far_switch_radius;
			bessel_cutoff = _bessel_cutoff;
		}
		if (_far_switch_radius >= 0)
		{
			mmm1d_params.far_switch_radius_2 = _far_switch_radius*_far_switch_radius;
			cuda_safe_mem( cudaMemcpyToSymbol(::far_switch_radius_2, &_far_switch_radius_2, sizeof(mmm1dgpu_real)) );
			far_switch_radius = _far_switch_radius;
		}
		if (_boxz > 0)
		{
			host_boxz = _boxz;
			cuda_safe_mem( cudaMemcpyToSymbol(::boxz, &_boxz, sizeof(mmm1dgpu_real)) );
			cuda_safe_mem( cudaMemcpyToSymbol(::uz, &_uz, sizeof(mmm1dgpu_real)) );
		}
		if (_coulomb_prefactor != 0)
		{
			cuda_safe_mem( cudaMemcpyToSymbol(::coulomb_prefactor, &_coulomb_prefactor, sizeof(mmm1dgpu_real)) );
			coulomb_prefactor = _coulomb_prefactor;
		}
		if (_bessel_cutoff > 0)
		{
			mmm1d_params.bessel_cutoff = _bessel_cutoff;
			cuda_safe_mem( cudaMemcpyToSymbol(::bessel_cutoff, &_bessel_cutoff, sizeof(int)) );
			bessel_cutoff = _bessel_cutoff;
		}
		if (_maxPWerror > 0)
		{
			mmm1d_params.maxPWerror = _maxPWerror;
			cuda_safe_mem( cudaMemcpyToSymbol(::maxPWerror, &_maxPWerror, sizeof(mmm1dgpu_real)) );
			maxPWerror = _maxPWerror;
		}
	}
	need_tune = true;
	
	// The changed parameters in mmm1d_params do not need to be broadcast: they are only accessed by the TCL print function (on node 0) when you call inter coulomb. The CUDA code only runs on node 0, so other nodes do not need the parameters. We couldn't broadcast from here anyway because set_params() might be called from inside computeForces() which is not a time at which the MPI loop on the slave nodes is waiting for broadcasts.
}

__global__ void forcesKernel(const __restrict__ mmm1dgpu_real *r, const __restrict__ mmm1dgpu_real *q, __restrict__ mmm1dgpu_real *force, int N, int pairs, int tStart = 0, int tStop = -1)
{
	if (tStop < 0)
		tStop = N*N;

	for (int tid = threadIdx.x + blockIdx.x * blockDim.x + tStart; tid < tStop; tid += blockDim.x * gridDim.x)
	{
		int p1 = tid%N, p2 = tid/N;
		mmm1dgpu_real x = r[3*p2] - r[3*p1], y = r[3*p2+1] - r[3*p1+1], z = r[3*p2+2] - r[3*p1+2];
		mmm1dgpu_real rxy2 = sqpow(x) + sqpow(y);
		mmm1dgpu_real rxy = sqrt(rxy2);
		mmm1dgpu_real sum_r = 0, sum_z = 0;
		
//		if (boxz <= 0.0) return; // otherwise we'd get into an infinite loop if we're not initialized correctly

		while (fabs(z) > boxz/2) // make sure we take the shortest distance
			z -= (z > 0? 1 : -1)*boxz;

		if (p1 == p2) // particle exerts no force on itself
		{
			rxy = 1; // so the division at the end doesn't fail with NaN (sum_r is 0 anyway)
		}
		else if (rxy2 <= far_switch_radius_2) // near formula
		{
			mmm1dgpu_real uzz = uz*z;
			mmm1dgpu_real uzr = uz*rxy;
			sum_z = dev_mod_psi_odd(0, uzz);
			mmm1dgpu_real uzrpow = uzr;
			for (int n = 1; n < device_n_modPsi; n++)
			{
				mmm1dgpu_real sum_r_old = sum_r;
				mmm1dgpu_real mpe = dev_mod_psi_even(n, uzz);
     			mmm1dgpu_real mpo = dev_mod_psi_odd(n, uzz);

     			sum_r += 2*n*mpe * uzrpow;
     			uzrpow *= uzr;
     			sum_z += mpo * uzrpow;
     			uzrpow *= uzr;

     			if (fabs(sum_r_old - sum_r) < maxPWerror)
					break;
			}

			sum_r *= sqpow(uz);
			sum_z *= sqpow(uz);

			sum_r += rxy*cbpow(rsqrt(rxy2+sqpow(z)));
			sum_r += rxy*cbpow(rsqrt(rxy2+sqpow(z+boxz)));
			sum_r += rxy*cbpow(rsqrt(rxy2+sqpow(z-boxz)));

			sum_z += z*cbpow(rsqrt(rxy2+sqpow(z)));
			sum_z += (z+boxz)*cbpow(rsqrt(rxy2+sqpow(z+boxz)));
			sum_z += (z-boxz)*cbpow(rsqrt(rxy2+sqpow(z-boxz)));

			if (rxy == 0) // particles at the same radial position only exert a force in z direction
			{
				rxy = 1;  // so the division at the end doesn't fail with NaN (sum_r is 0 anyway)
			}
		}
		else // far formula
		{
			for (int p = 1; p < bessel_cutoff; p++)
			{
				mmm1dgpu_real arg = C_2PIf*uz*p;
				sum_r += p*dev_K1(arg*rxy)*cos(arg*z);
				sum_z += p*dev_K0(arg*rxy)*sin(arg*z);
			}
			sum_r *= sqpow(uz)*4*C_2PIf;
			sum_z *= sqpow(uz)*4*C_2PIf;
			sum_r += 2*uz/rxy;
		}

		mmm1dgpu_real pref = coulomb_prefactor*q[p1]*q[p2];
		if (pairs)
		{
			force[3*(p1+p2*N-tStart)] = pref*sum_r/rxy*x;
			force[3*(p1+p2*N-tStart)+1] = pref*sum_r/rxy*y;
			force[3*(p1+p2*N-tStart)+2] = pref*sum_z;
		}
		else
		{
#ifdef ELECTROSTATICS_GPU_DOUBLE_PRECISION
			atomicadd8(&force[3*p2], pref*sum_r/rxy*x);
			atomicadd8(&force[3*p2+1], pref*sum_r/rxy*y);
			atomicadd8(&force[3*p2+2], pref*sum_z);
#else
			atomicadd(&force[3*p2], pref*sum_r/rxy*x);
			atomicadd(&force[3*p2+1], pref*sum_r/rxy*y);
			atomicadd(&force[3*p2+2], pref*sum_z);
#endif
		}
	}
}

__global__ void energiesKernel(const __restrict__ mmm1dgpu_real *r, const __restrict__ mmm1dgpu_real *q, __restrict__ mmm1dgpu_real *energy, int N, int pairs, int tStart = 0, int tStop = -1)
{
	if (tStop < 0)
		tStop = N*N;

	extern __shared__ mmm1dgpu_real partialsums[];
	if (!pairs)
	{
		partialsums[threadIdx.x] = 0;
		__syncthreads();
	}
	for (int tid = threadIdx.x + blockIdx.x * blockDim.x + tStart; tid < tStop; tid += blockDim.x * gridDim.x)
	{
		int p1 = tid%N, p2 = tid/N;
		mmm1dgpu_real z = r[3*p2+2] - r[3*p1+2];
		mmm1dgpu_real rxy2 = sqpow(r[3*p2] - r[3*p1]) + sqpow(r[3*p2+1] - r[3*p1+1]);
		mmm1dgpu_real rxy = sqrt(rxy2);
		mmm1dgpu_real sum_e = 0;

//		if (boxz <= 0.0) return; // otherwise we'd get into an infinite loop if we're not initialized correctly

		while (fabs(z) > boxz/2) // make sure we take the shortest distance
			z -= (z > 0? 1 : -1)*boxz;

		if (p1 == p2) // particle exerts no force on itself
		{
		}
		else if (rxy2 <= far_switch_radius_2) // near formula
		{
			mmm1dgpu_real uzz = uz*z;
			mmm1dgpu_real uzr2 = sqpow(uz*rxy);
			mmm1dgpu_real uzrpow = uzr2;
			sum_e = dev_mod_psi_even(0, uzz);
			for (int n = 1; n < device_n_modPsi; n++)
			{
				mmm1dgpu_real sum_e_old = sum_e;
				mmm1dgpu_real mpe = dev_mod_psi_even(n, uzz);
     			sum_e += mpe * uzrpow;
     			uzrpow *= uzr2;
				
				if (fabs(sum_e_old - sum_e) < maxPWerror)
					break;
			}

			sum_e *= -1*uz;
			sum_e -= 2*uz*C_GAMMAf;
			sum_e += rsqrt(rxy2+sqpow(z));
			sum_e += rsqrt(rxy2+sqpow(z+boxz));
			sum_e += rsqrt(rxy2+sqpow(z-boxz));
		}
		else // far formula
		{
			sum_e = -(log(rxy*uz/2) + C_GAMMAf)/2;
			for (int p = 1; p < bessel_cutoff; p++)
			{
				mmm1dgpu_real arg = C_2PIf*uz*p;
				sum_e += dev_K0(arg*rxy)*cos(arg*z);
			}
			sum_e *= uz*4;
		}

		if (pairs)
		{
			energy[p1+p2*N-tStart] = coulomb_prefactor*q[p1]*q[p2]*sum_e;
		}
		else
		{
			partialsums[threadIdx.x] += coulomb_prefactor*q[p1]*q[p2]*sum_e;
		}
	}
	if (!pairs)
	{
		sumReduction(partialsums, &energy[blockIdx.x]);
	}
}

__global__ void vectorReductionKernel(mmm1dgpu_real *src, mmm1dgpu_real *dst, int N, int tStart = 0, int tStop = -1)
{
	if (tStop < 0)
		tStop = N*N;

	for (int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < N; tid += blockDim.x * gridDim.x)
	{
		int offset = ((tid + (tStart % N)) % N);
		
		for (int i = 0; tid+i*N < (tStop - tStart); i++)
		{
			#pragma unroll 3
			for (int d = 0; d<3; d++)
			{
				dst[3*offset+d] -= src[3*(tid+i*N)+d];
			}
		}
	}
}

void Mmm1dgpuForce::computeForces(SystemInterface &s)
{
	if (coulomb.method != COULOMB_MMM1D_GPU) // MMM1DGPU was disabled. nobody cares about our calculations anymore
	{
		std::cerr << "MMM1D: coulomb.method has been changed, skipping calculation" << std::endl;
		return;
	}
	setup(s);

	if (pairs < 0)
	{
		std::cerr << "MMM1D was not initialized correctly" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (pairs) // if we calculate force pairs, we need to reduce them to forces
	{
		int blocksRed = s.npart_gpu()/numThreads+1;
		KERNELCALL(forcesKernel,numBlocks(s),numThreads,(s.rGpuBegin(), s.qGpuBegin(), dev_forcePairs, s.npart_gpu(), pairs))
		KERNELCALL(vectorReductionKernel,blocksRed,numThreads,(dev_forcePairs, s.fGpuBegin(), s.npart_gpu()))
	}
	else
	{
		KERNELCALL(forcesKernel,numBlocks(s),numThreads,(s.rGpuBegin(), s.qGpuBegin(), s.fGpuBegin(), s.npart_gpu(), pairs))
	}
}

__global__ void scaleAndAddKernel(mmm1dgpu_real *dst, mmm1dgpu_real *src, int N, mmm1dgpu_real factor)
{
	for (int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < N; tid += blockDim.x * gridDim.x)
	{
		dst[tid] += src[tid]*factor;
	}
}

void Mmm1dgpuForce::computeEnergy(SystemInterface &s)
{
	if (coulomb.method != COULOMB_MMM1D_GPU) // MMM1DGPU was disabled. nobody cares about our calculations anymore
	{
		std::cerr << "MMM1D: coulomb.method has been changed, skipping calculation" << std::endl;
		return;
	}
	setup(s);

	if (pairs < 0)
	{
		std::cerr << "MMM1D was not initialized correctly" << std::endl;
		exit(EXIT_FAILURE);
	}
	int shared = numThreads*sizeof(mmm1dgpu_real);

	KERNELCALL_shared(energiesKernel,numBlocks(s),numThreads,shared,(s.rGpuBegin(), s.qGpuBegin(), dev_energyBlocks, s.npart_gpu(), 0));
	KERNELCALL_shared(sumKernel,1,numThreads,shared,(dev_energyBlocks, numBlocks(s)));
	KERNELCALL(scaleAndAddKernel,1,1,(&(((CUDA_energy*)s.eGpu())->coulomb), &dev_energyBlocks[0], 1, 0.5)); // we have counted every interaction twice, so halve the total energy
}

float Mmm1dgpuForce::force_benchmark(SystemInterface &s)
{
	cudaEvent_t eventStart, eventStop;
	float elapsedTime;
	mmm1dgpu_real *dev_f_benchmark;

	cuda_safe_mem( cudaMalloc((void**)&dev_f_benchmark, 3*s.npart_gpu()*sizeof(mmm1dgpu_real)) );
	cuda_safe_mem( cudaEventCreate(&eventStart) );
	cuda_safe_mem( cudaEventCreate(&eventStop) );
	cuda_safe_mem( cudaEventRecord(eventStart, stream[0]) );
	KERNELCALL(forcesKernel,numBlocks(s),numThreads,(s.rGpuBegin(), s.qGpuBegin(), dev_f_benchmark, s.npart_gpu(), 0))
	cuda_safe_mem( cudaEventRecord(eventStop, stream[0]) );
	cuda_safe_mem( cudaEventSynchronize(eventStop) );
	cuda_safe_mem( cudaEventElapsedTime(&elapsedTime, eventStart, eventStop) );
	cuda_safe_mem( cudaEventDestroy(eventStart) );
	cuda_safe_mem( cudaEventDestroy(eventStop) );
	cuda_safe_mem( cudaFree(dev_f_benchmark));

	return elapsedTime;
}

#endif /* MMM1D_GPU */
