#ifndef CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP
#define CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP

#include "virtual_sites/lb_inertialess_tracers_cuda_interface.hpp"

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataInput &d,
               const unsigned int /* version */) {
  ar &d.pos;
  ar &d.f;
  ar &d.is_virtual;
}

template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataOutput &d,
               const unsigned int /* version */) {
  ar &d.v;
}
}
}

#endif
