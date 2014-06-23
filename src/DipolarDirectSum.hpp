#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

#include "SystemInterface.hpp"
#include <iostream>

void DipolarDirectSum_kernel_wrapper(float x, float y, float z, float k,
		     int n, float *pos, float *dip, float *f, float* torque);

class DipolarDirectSum {
public:
  DipolarDirectSum(SystemInterface &s) {
    if(!s.requestFGpu())
      std::cerr << "DipolarDirectSum needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "DipolarDirectSum needs access to positions on GPU!" << std::endl;

  }; 
  void calc(SystemInterface &s) {
    DipolarDirectSum_kernel_wrapper(x,y,z,k,s.npart_gpu(),
					 s.rGpuBegin(), s.fGpuBegin());
  };
protected:
  float x,y,z;
  float k;
};

extern DipolarDirectSum *dipolarDirectSum;

#endif
