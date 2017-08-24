#ifndef CORE_UTILS_SERIALIZATION_PARTICLE_HPP
#define CORE_UTILS_SERIALIZATION_PARTICLE_HPP

#include "core/particle_data.hpp"
#include "core/utils/serialization/List.hpp"
#include <boost/serialization/vector.hpp>

namespace boost {
namespace serialization {
/* Pod serialize for Particle */
template <typename Archive>
void load(Archive &ar, Particle &p, const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(Particle));
  new (&(p.bl)) IntList(p.bl.size());
  ar >> p.bl;

#ifdef EXCLUSIONS
  new (&(p.el)) IntList(p.el.size());
  ar >> p.el;
#endif
}

template <typename Archive>
void save(Archive &ar, Particle const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar << make_array(reinterpret_cast<char const *>(&p), sizeof(Particle));
  ar << p.bl;

#ifdef EXCLUSIONS
  ar << p.el;
#endif
}

template <class Archive>
void serialize(Archive &ar, Particle &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
