#ifndef CUFFT_WRAPPER_HPP
#define CUFFT_WRAPPER_HPP

#if defined(__HIPCC__) // AMD or Nvidia via HIP

#include <hipfft.h>

#define cufftComplex hipfftComplex
#define cufftDestroy hipfftDestroy
#define cufftDoubleComplex hipfftDoubleComplex
#define cufftExecC2R hipfftExecC2R
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecR2C hipfftExecR2C
#define cufftExecZ2D hipfftExecZ2D
#define cufftHandle hipfftHandle
#define cufftPlan3d hipfftPlan3d
#define cufftReal hipfftReal
#define cufftSetStream hipfftSetStream
#define CUFFT_R2C HIPFFT_R2C
#define CUFFT_C2R HIPFFT_C2R
#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_SUCCESS HIPFFT_SUCCESS

#else // Nvidia via CUDA

#include <cufft.h>

#endif

#endif // CUFFT_WRAPPER_HPP
