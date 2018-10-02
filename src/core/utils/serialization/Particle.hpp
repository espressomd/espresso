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
  p.~Particle();
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(Particle));
  p.p = new ParticleProperties;
  ar >> make_array(reinterpret_cast<char *>(p.p), sizeof(ParticleProperties));
  p.m = new ParticleMomentum;
  ar >> make_array(reinterpret_cast<char *>(p.m), sizeof(ParticleMomentum));
  p.f = new ParticleForce;
  ar >> make_array(reinterpret_cast<char *>(p.f), sizeof(ParticleForce));
  p.l = new ParticleLocal;
  ar >> make_array(reinterpret_cast<char *>(p.l), sizeof(ParticleLocal));
#ifdef LB
  p.lc = new ParticleLatticeCoupling;
  ar >> make_array(reinterpret_cast<char *>(p.lc), sizeof(ParticleLatticeCoupling));
#endif
#ifdef ENGINE
  p.swim = new ParticleParametersSwimming;
  ar >> make_array(reinterpret_cast<char *>(p.swim), sizeof(ParticleParametersSwimming));
#endif
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
  ar << make_array(reinterpret_cast<char const *>(p.p), sizeof(ParticleProperties));
  ar << make_array(reinterpret_cast<char const *>(p.m), sizeof(ParticleMomentum));
  ar << make_array(reinterpret_cast<char const *>(p.f), sizeof(ParticleForce));
  ar << make_array(reinterpret_cast<char const *>(p.l), sizeof(ParticleLocal));
#ifdef LB
  ar << make_array(reinterpret_cast<char const *>(p.lc), sizeof(ParticleLatticeCoupling));
#endif
#ifdef ENGINE
  ar << make_array(reinterpret_cast<char const *>(p.swim), sizeof(ParticleParametersSwimming));
#endif
  ar << p.bl;
#ifdef EXCLUSIONS
  ar << p.el;
#endif
}

template <class Archive>
void serialize(Archive &ar, Particle &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
} // namespace serialization
} // namespace boost

#endif
