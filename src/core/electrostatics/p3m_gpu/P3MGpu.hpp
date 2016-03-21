#ifndef __ELECTROSTATICS_P3M_GPU_P3M_GPU_HPP
#define __ELECTROSTATICS_P3M_GPU_P3M_GPU_HPP

#include <thrust/device_vector.h>

#include "actor/Actor.hpp"

namespace Electrostatics {
namespace P3MGpu {

template<typename Scalar = double, typename cufft_complex>
class P3MGpu : public Actor {
  void computeForces(SystemInterface &s);
  void computeEnergy(SystemInterface &s);  
 private:
  thrust::device_vector<Scalar> G_hat;
  thrust::device_vector<cufft_complex> charge_mesh; 
  thrust::device_vector<cufft_complex> force_mesh_x;
  thrust::device_vector<cufft_complex> force_mesh_y;
  thrust::device_vector<cufft_complex> force_mesh_z; 
};

}

}

#endif
