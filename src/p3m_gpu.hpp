#ifndef _P3M_GPU_H
#define _P3M_GPU_H

//NOTE :if one wants to use doubles it requires cuda compute capability 1.3
#define _P3M_GPU_FLOAT
//#define _P3M_GPU_REAL_DOUBLE

#ifdef _P3M_GPU_FLOAT
#define REAL_TYPE float
#define CUFFT_TYPE_COMPLEX cufftComplex
#define CUFFT_FFT cufftExecC2C
#define CUFFT_PLAN_FLAG CUFFT_C2C
#endif

#ifdef _P3M_GPU_REAL_DOUBLE
#define REAL_TYPE double
#define CUFFT_TYPE_COMPLEX cufftDoubleComplex
#define CUFFT_FFT cufftExecZ2Z
#define CUFFT_PLAN_FLAG CUFFT_Z2Z
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void p3m_gpu_init(int cao, int mesh, REAL_TYPE alpha, REAL_TYPE box);
  void p3m_gpu_add_farfield_force();
#ifdef __cplusplus
}
#endif
#endif
