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
      float bjerrum_length, kT;
      int dim_x, dim_y, dim_z;
      float agrid;
    };

    struct Parameters : public InputParameters {
      Parameters() {}
      Parameters(InputParameters& inputParameters) : InputParameters(inputParameters)
      {
        charge_potential = 0;
        greensfcn = 0;
      }

      cufftComplex *charge_potential;
      cufftReal *greensfcn;
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


#endif
