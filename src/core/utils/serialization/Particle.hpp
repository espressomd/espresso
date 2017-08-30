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
  //  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(Particle));
  ar >> make_array(&p, 1);
  new (&(p.bl)) IntList(p.bl.size());
  ar >> p.bl;
#ifdef EXCLUSIONS
  new (&(p.el)) IntList(p.el.size());
  ar >> p.el;
#endif

#ifdef LB_BOUNDARIES_GPU
  new (&(p.p.anchors)) std::vector<float>();
  ar >> p.p.anchors;
  new (&(p.p.anchors_out)) std::vector<float>();
  ar >> p.p.anchors_out;
#endif
}

template <typename Archive>
void save(Archive &ar, Particle const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  //ar << make_array(reinterpret_cast<char const *>(&p), sizeof(Particle));
  ar << make_array(&p, 1);
  ar << p.bl;
#ifdef EXCLUSIONS
  ar << p.el;
#endif

#ifdef LB_BOUNDARIES_GPU
  ar << p.p.anchors;
  ar << p.p.anchors_out;
#endif
}

template <class Archive>
void serialize(Archive &ar, Particle &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
