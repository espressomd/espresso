#ifndef CORE_UTILS_SERIALIZATION_PARTICLE_LIST_HPP
#define CORE_UTILS_SERIALIZATION_PARTICLE_LIST_HPP

#include "core/utils/serialization/Particle.hpp"

namespace boost {
namespace serialization {
template <class Archive>
void load(Archive &ar, ParticleList &pl,
          const unsigned int /* version */) {
  int size;
  ar >> size;

  realloc_particlelist(&pl, pl.n = size);
  for (int i = 0; i < size; i++) {
    ar >> pl.part[i];
  }
}

template <class Archive>
void save(Archive &ar, ParticleList const &pl,
          const unsigned int /* version */) {
  ar << pl.n;

  for (int i = 0; i < pl.n; i++) {
    ar << pl.part[i];
  }
}

template <typename Archive>
void serialize(Archive &ar, ParticleList &pl, unsigned int file_version) {
  split_free(ar, pl, file_version);
}
}
}

#endif
