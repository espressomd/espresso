#ifndef SERIALIZATION_BONDED_IA_PARAMETERS_HPP
#define SERIALIZATION_BONDED_IA_PARAMETERS_HPP

#include "core/interaction_data.hpp"

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, Bonded_ia_parameters &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(Bonded_ia_parameters));
}

template <typename Archive>
void save(Archive &ar, Bonded_ia_parameters const &p,
          const unsigned int /* file_version */) {
  ar << make_array(reinterpret_cast<char const *>(&p), sizeof(Bonded_ia_parameters));
}

template <class Archive>
void serialize(Archive &ar, Bonded_ia_parameters &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
