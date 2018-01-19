#ifndef SERIALIZATION_IA_PARAMETERS_HPP
#define SERIALIZATION_IA_PARAMETERS_HPP

#include "core/interaction_data.hpp"

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, IA_parameters &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(IA_parameters));

#ifdef TABULATED
  TabulatedPotential tab;
  ar >> tab;

  new (&(p.TAB)) TabulatedPotential(std::move(tab));
#endif
}

template <typename Archive>
void save(Archive &ar, IA_parameters const &p,
          const unsigned int /* file_version */) {
  ar << make_array(reinterpret_cast<char const *>(&p), sizeof(IA_parameters));

#ifdef TABULATED
  ar << p.TAB;
#endif
}

template <class Archive>
void serialize(Archive &ar, IA_parameters &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
