#ifndef CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP
#define CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP

#include "ibm_cuda_interface.hpp"

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataInput &d,
               const unsigned int /* version */) {
  ar &d.pos;
  ar &d.f;
  ar &d.isVirtual;
}

template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataOutput &d,
               const unsigned int /* version */) {
  ar &d.v;
}
}
}

#endif
