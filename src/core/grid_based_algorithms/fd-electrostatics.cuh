/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef _FD_ELECTROSTATICS_HPP
#define _FD_ELECTROSTATICS_HPP

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
    Parameters(InputParameters &inputParameters)
        : InputParameters(inputParameters) {
      charge_potential = 0;
      greensfcn = 0;
      dim_x_padded = (inputParameters.dim_x / 2 + 1) * 2;
    }

    hipfftComplex *charge_potential;
    hipfftReal *greensfcn;
    int dim_x_padded;
  };

  struct Grid {
    float *grid;
    int dim_x;
    int dim_y;
    int dim_z;
    float agrid;
  };

  ~FdElectrostatics();
  FdElectrostatics(InputParameters inputParameters, hipStream_t stream);
  void calculatePotential();
  void calculatePotential(hipfftComplex *charge_potential);
  Grid getGrid();

private:
  Parameters parameters;
  hipStream_t cuda_stream;
  hipfftHandle plan_fft;
  hipfftHandle plan_ifft;
  bool initialized;
};

// extern __device__ __constant__ FdElectrostatics::Parameters
// fde_parameters_gpu;

__device__ hipfftReal fde_getNode(int x, int y, int z);
__device__ hipfftReal fde_getNode(int i);
__device__ void fde_setNode(int x, int y, int z, hipfftReal value);
__device__ void fde_setNode(int i, hipfftReal value);

#endif
