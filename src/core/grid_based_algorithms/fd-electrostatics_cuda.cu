#include <hip/hip_runtime.h>

// TODO: throw exceptions upon errors initialization

#include <hipfft.h>
#include "grid_based_algorithms/fd-electrostatics.cuh"
#include "cuda_utils.hpp"
#include <string>
//#include <cuda_interface.hpp>
#include <cstdio>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

__global__ void createGreensfcn();
__global__ void multiplyGreensfcn(hipfftComplex *charge_potential);

__device__ __constant__ FdElectrostatics::Parameters fde_parameters_gpu[1];

__device__ unsigned int fde_getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
         threadIdx.x;
}

__device__ hipfftReal fde_getNode(int x, int y, int z) {
  hipfftReal *field =
      reinterpret_cast<hipfftReal *>(fde_parameters_gpu->charge_potential);
  return field[fde_parameters_gpu->dim_y * fde_parameters_gpu->dim_x_padded * z +
               fde_parameters_gpu->dim_x_padded * y + x];
}

__device__ void fde_setNode(int x, int y, int z, hipfftReal value) {
  hipfftReal *field =
      reinterpret_cast<hipfftReal *>(fde_parameters_gpu->charge_potential);
  field[fde_parameters_gpu->dim_y * fde_parameters_gpu->dim_x_padded * z +
        fde_parameters_gpu->dim_x_padded * y + x] = value;
}

__device__ hipfftReal fde_getNode(int i) {
  int x = i % fde_parameters_gpu->dim_x_padded;
  i /= fde_parameters_gpu->dim_x_padded;
  int y = i % fde_parameters_gpu->dim_y;
  int z = i / fde_parameters_gpu->dim_y;
  return fde_getNode(x, y, z);
}

__device__ void fde_setNode(int i, hipfftReal value) {
  int x = i % fde_parameters_gpu->dim_x_padded;
  i /= fde_parameters_gpu->dim_x_padded;
  int y = i % fde_parameters_gpu->dim_y;
  int z = i / fde_parameters_gpu->dim_y;
  fde_setNode(x, y, z, value);
}

FdElectrostatics::~FdElectrostatics() {
  hipfftDestroy(plan_ifft);
  hipfftDestroy(plan_fft);

  cuda_safe_mem(hipFree(parameters.greensfcn));
  cuda_safe_mem(hipFree(parameters.charge_potential));
}

FdElectrostatics::FdElectrostatics(InputParameters inputParameters,
                                   hipStream_t stream)
    : parameters(inputParameters), cuda_stream(stream) {
  cuda_safe_mem(hipMalloc((void **)&parameters.charge_potential,
                           sizeof(hipfftComplex) * parameters.dim_z *
                               parameters.dim_y * (parameters.dim_x / 2 + 1)));

  cuda_safe_mem(hipMalloc((void **)&parameters.greensfcn,
                           sizeof(hipfftReal) * parameters.dim_z *
                               parameters.dim_y * (parameters.dim_x / 2 + 1)));

  if (hipGetLastError() != hipSuccess) {
    throw "Failed to allocate\n";
  }

  cuda_safe_mem(
      hipMemcpyToSymbol(HIP_SYMBOL(fde_parameters_gpu), &parameters, sizeof(Parameters)));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (parameters.dim_z * parameters.dim_y * (parameters.dim_x / 2 + 1) +
       threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
  KERNELCALL_stream(createGreensfcn, dim_grid, threads_per_block, stream, ());

  /* create 3D FFT plans */

  if (hipfftPlan3d(&plan_fft, parameters.dim_z, parameters.dim_y,
                  parameters.dim_x, HIPFFT_R2C) != HIPFFT_SUCCESS) {
    throw std::string("Unable to create fft plan");
  }

  if (hipfftSetStream(plan_fft, cuda_stream) != HIPFFT_SUCCESS) {
    throw std::string("Unable to assign FFT to cuda stream");
  }

  if (hipfftPlan3d(&plan_ifft, parameters.dim_z, parameters.dim_y,
                  parameters.dim_x, HIPFFT_C2R) != HIPFFT_SUCCESS) {
    throw std::string("Unable to create ifft plan");
  }

  if (hipfftSetStream(plan_ifft, cuda_stream) != HIPFFT_SUCCESS) {
    throw std::string("Unable to assign FFT to cuda stream");
  }

  initialized = true;
}

__global__ void createGreensfcn() {
  unsigned int index = fde_getThreadIndex();
  unsigned int tmp;
  unsigned int coord[3];

  coord[0] = index % (fde_parameters_gpu->dim_x / 2 + 1);
  tmp = index / (fde_parameters_gpu->dim_x / 2 + 1);
  coord[1] = tmp % fde_parameters_gpu->dim_y;
  coord[2] = tmp / fde_parameters_gpu->dim_y;

  if (index < fde_parameters_gpu->dim_z * fde_parameters_gpu->dim_y *
                  (fde_parameters_gpu->dim_x / 2 + 1)) {

    if (index == 0) {
      // setting 0th Fourier mode to 0 enforces charge neutrality
      fde_parameters_gpu->greensfcn[index] = 0.0f;
    } else {
      fde_parameters_gpu->greensfcn[index] =
          -4.0f * PI_FLOAT * fde_parameters_gpu->prefactor *
          fde_parameters_gpu->agrid * fde_parameters_gpu->agrid * 0.5f /
          (cos(2.0f * PI_FLOAT * coord[0] /
               (hipfftReal)fde_parameters_gpu->dim_x) +
           cos(2.0f * PI_FLOAT * coord[1] /
               (hipfftReal)fde_parameters_gpu->dim_y) +
           cos(2.0f * PI_FLOAT * coord[2] /
               (hipfftReal)fde_parameters_gpu->dim_z) -
           3.0f) /
          (fde_parameters_gpu->dim_x * fde_parameters_gpu->dim_y *
           fde_parameters_gpu->dim_z);
    }

    // fde_parameters_gpu->greensfcn[index] = 0.0f; //TODO delete
  }
}

__global__ void multiplyGreensfcn(hipfftComplex *charge_potential) {

  unsigned int index = fde_getThreadIndex();

  if (index < fde_parameters_gpu->dim_z * fde_parameters_gpu->dim_y *
                  (fde_parameters_gpu->dim_x / 2 + 1)) {
    charge_potential[index].x *= fde_parameters_gpu->greensfcn[index];
    charge_potential[index].y *= fde_parameters_gpu->greensfcn[index];
  }
}

void FdElectrostatics::calculatePotential() {
  calculatePotential(parameters.charge_potential);
}

void FdElectrostatics::calculatePotential(hipfftComplex *charge_potential) {

  if (hipfftExecR2C(plan_fft, (hipfftReal *)charge_potential, charge_potential) !=
      HIPFFT_SUCCESS) {

    fprintf(stderr, "ERROR: Unable to execute FFT plan\n");
  }

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (parameters.dim_z * parameters.dim_y * (parameters.dim_x / 2 + 1) +
       threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(multiplyGreensfcn, dim_grid, threads_per_block,
             (charge_potential));

  if (hipfftExecC2R(plan_ifft, charge_potential,
                   (hipfftReal *)charge_potential) != HIPFFT_SUCCESS) {

    fprintf(stderr, "ERROR: Unable to execute iFFT plan\n");
  }
}

FdElectrostatics::Grid FdElectrostatics::getGrid() {
  Grid g = {(float *)parameters.charge_potential, parameters.dim_x,
            parameters.dim_y, parameters.dim_z, parameters.agrid};
  return g;
}
