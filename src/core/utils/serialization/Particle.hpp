#ifndef CORE_UTILS_SERIALIZATION_PARTICLE_HPP
#define CORE_UTILS_SERIALIZATION_PARTICLE_HPP

#include "core/particle_data.hpp"

namespace boost {
namespace serialization {
/* Pod serialize for Particle */
template <typename Archive>
void serialize(Archive &ar, Particle &p, unsigned int) {
  /* Cruel but effective */
  ar &make_array(reinterpret_cast<char *>(&p), sizeof(Particle));
}
}
}

#endif
