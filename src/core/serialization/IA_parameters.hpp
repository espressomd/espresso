/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SERIALIZATION_IA_PARAMETERS_HPP
#define SERIALIZATION_IA_PARAMETERS_HPP

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

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
} // namespace serialization
} // namespace boost

#endif
