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
#ifndef CORE_UTILS_SERIALIZATION_PARTICLE_LIST_HPP
#define CORE_UTILS_SERIALIZATION_PARTICLE_LIST_HPP

#include "Particle.hpp"

namespace boost {
namespace serialization {
template <class Archive>
void load(Archive &ar, ParticleList &pl, const unsigned int /* version */) {
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
} // namespace serialization
} // namespace boost

#endif
