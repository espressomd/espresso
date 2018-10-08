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
  ar >> make_array(reinterpret_cast<char *>(&p.r), sizeof(ParticlePosition));
  p.p.reset(new ParticleProperties{});
  ar >> make_array(reinterpret_cast<char *>(p.p.get()), sizeof(ParticleProperties));
  p.m.reset(new ParticleMomentum{});
  ar >> make_array(reinterpret_cast<char *>(p.m.get()), sizeof(ParticleMomentum));
  ar >> make_array(reinterpret_cast<char *>(&p.f), sizeof(ParticleForce));
  p.l.reset(new ParticleLocal{});
  ar >> make_array(reinterpret_cast<char *>(p.l.get()), sizeof(ParticleLocal));
#ifdef LB
  p.lc.reset(new ParticleLatticeCoupling{});
  ar >> make_array(reinterpret_cast<char *>(p.lc.get()), sizeof(ParticleLatticeCoupling));
#endif
#ifdef ENGINE
  p.swim.reset(new ParticleParametersSwimming{});
  ar >> make_array(reinterpret_cast<char *>(p.swim.get()), sizeof(ParticleParametersSwimming));
#endif
  uint32_t n;
  ar >> n;
  new (&(p.bl)) IntList(n);
  ar >> p.bl;
#ifdef EXCLUSIONS
  ar >> n;
  new (&(p.el)) IntList(n);
  ar >> p.el;
#endif
}

template <typename Archive>
void save(Archive &ar, Particle const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar << make_array(reinterpret_cast<char const *>(&p.r), sizeof(ParticlePosition));
  ar << make_array(reinterpret_cast<char const *>(p.p.get()), sizeof(ParticleProperties));
  ar << make_array(reinterpret_cast<char const *>(p.m.get()), sizeof(ParticleMomentum));
  ar << make_array(reinterpret_cast<char const *>(&p.f), sizeof(ParticleForce));
  ar << make_array(reinterpret_cast<char const *>(p.l.get()), sizeof(ParticleLocal));
#ifdef LB
  ar << make_array(reinterpret_cast<char const *>(p.lc.get()), sizeof(ParticleLatticeCoupling));
#endif
#ifdef ENGINE
  ar << make_array(reinterpret_cast<char const *>(p.swim.get()), sizeof(ParticleParametersSwimming));
#endif
  ar << p.bl.n;
  ar << p.bl;
#ifdef EXCLUSIONS
  ar << p.el.n;
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
