#ifndef OBSERVABLES_PARTICLE_TRAITS
#define OBSERVABLES_PARTICLE_TRAITS

#include "Particle.hpp"
#include "config.hpp"

namespace GenObs {
template <> struct traits<const ::Particle *> {
  using Particle = ::Particle;

  auto position(Particle const *const p) const { return p->r.p; }
  auto velocity(Particle const *const p) const { return p->m.v; }
  auto mass(Particle const *const p) const {
#ifdef VIRTUAL_SITES
    if (p->p.is_virtual)
      return decltype(p->p.mass){};
#endif
    return p->p.mass;
  }
  auto charge(Particle const *const p) const { return p->p.q; }
  auto force(Particle const *const p) const {
#ifdef VIRTUAL_SITES
    if (p->p.is_virtual)
      return decltype(p->f.f){};
#endif
    return p->f.f;
  }
  auto dipole_moment(Particle const *const p) const { return p->calc_dip(); }
};

} // namespace GenObs

#endif
