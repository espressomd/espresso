#ifndef OBSERVABLES_PARTICLE_TRAITS
#define OBSERVABLES_PARTICLE_TRAITS

#include "Particle.hpp"
#include "config.hpp"

namespace GenObs {
template <> struct traits<Particle> {
  static auto position(Particle const &p) { return p.r.p; }
  static auto velocity(Particle const &p) { return p.m.v; }
  static auto mass(Particle const &p) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      return decltype(p.p.mass){};
#endif
    return p.p.mass;
  }
  static auto charge(Particle const &p) { return p.p.q; }
  static auto force(Particle const &p) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      return decltype(p.f.f){};
#endif
    return p.f.f;
  }
  static auto dipole_moment(Particle const &p) {
#if defined(ROTATION) && defined(DIPOLES)
    return p.calc_dip();
#else
    return Utils::Vector3d{};
#endif
  }
};

} // namespace GenObs

#endif
