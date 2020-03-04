#ifndef TRAITS_HPP
#define TRAITS_HPP

#include <particle/particle.hpp>

namespace Traits {
namespace Particle {

struct Charge {
  auto operator()(::Particle::Particle const &p) { return p.q; }
};

struct Position {
  auto operator()(::Particle::Particle const &p) { return p.x; }
};

struct Velocity {
  auto operator()(::Particle::Particle const &p) { return p.v; }
};

struct Mass {
  auto operator()(::Particle::Particle const &p) { return p.m; }
};
} // namespace Particle
} // namespace Traits
#endif
