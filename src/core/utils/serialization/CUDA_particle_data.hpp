#ifndef UTILS_SERIALIZATION_CUDA_PARTICLE_DATA_HPP
#define UTILS_SERIALIZATION_CUDA_PARTICLE_DATA_HPP

#include "core/cuda_interface.hpp"

#include <boost/serialization/array.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/split_free.hpp>

BOOST_IS_BITWISE_SERIALIZABLE(CUDA_particle_data)

 namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, CUDA_particle_data &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(CUDA_particle_data));
}

template <typename Archive>
void save(Archive &ar, CUDA_particle_data const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar << make_array(reinterpret_cast<char const *>(&p),
                   sizeof(CUDA_particle_data));
}

template <class Archive>
void serialize(Archive &ar, CUDA_particle_data &p,
               const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
