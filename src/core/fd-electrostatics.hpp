#ifndef _FD_ELECTROSTATICS_HPP
#define _FD_ELECTROSTATICS_HPP


#ifdef __CUDACC__

#include <cuda.h>
#include <cufft.h>

#else

typedef void cufftComplex;
typedef void cufftReal;

#endif

#define PI_FLOAT 3.14159265358979323846f


class FdElectrostatics {
  public:
    struct InputParameters {
      float prefactor;
      int dim_x, dim_y, dim_z;
      float agrid;
    };

    struct Parameters : public InputParameters {
      Parameters() {}
      Parameters(InputParameters& inputParameters) : InputParameters(inputParameters)
      {
        charge_potential = 0;
        greensfcn = 0;
        dim_x_padded = (inputParameters.dim_x/2+1)*2;
      }

      cufftComplex *charge_potential;
      cufftReal *greensfcn;
      int dim_x_padded;
    };

    struct Grid {
      float* grid;
      int dim_x;
      int dim_y;
      int dim_z;
      float agrid;
    };

    ~FdElectrostatics();
    FdElectrostatics(InputParameters inputParameters, cudaStream_t stream);
    void calculatePotential();
    void calculatePotential(cufftComplex *charge_potential);
    Grid getGrid();

  private:
    Parameters parameters;
    cudaStream_t cuda_stream;
    cufftHandle plan_fft;
    cufftHandle plan_ifft;
    bool initialized;
};

#ifdef __CUDACC__

//extern __device__ __constant__ FdElectrostatics::Parameters fde_parameters_gpu;

__device__ cufftReal fde_getNode(int x, int y, int z);
__device__ cufftReal fde_getNode(int i);
__device__ void fde_setNode(int x, int y, int z, cufftReal value);
__device__ void fde_setNode(int i, cufftReal value);

#endif //__CUDACC__


#endif
