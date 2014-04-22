#include "config.hpp"

#ifdef HARMONICFORCE

#include "SystemInterface.hpp"
#include <iostream>

void HarmonicForce_kernel_wrapper(float x, float y, float z, float k,
		     int n, float *pos, float *f);

class HarmonicForce {
public:
  HarmonicForce(float x1, float x2, float x3, float _k, SystemInterface &s) {
    x = x1;
    y = x2;
    z = x3;
    k = _k;

    if(!s.requestFGpu())
      std::cerr << "HarmonicForce needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "HarmonicForce needs access to positions on GPU!" << std::endl;

  }; 
  void calc(SystemInterface &s) {
    HarmonicForce_kernel_wrapper(x,y,z,k,s.npart_gpu(),
					 s.rGpuBegin(), s.fGpuBegin());
  };
protected:
  float x,y,z;
  float k;
};

extern HarmonicForce *harmonicForce;

#endif
