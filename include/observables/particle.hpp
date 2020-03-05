#ifndef OBSERVABLES_PARTICLE_HPP
#define OBSERVABLES_PARTICLE_HPP

namespace Particle {

struct Particle {
  double q; // charge
  double x; // position
  double v; // velocity
  double m; // mass
  bool is_virtual = false;
};
} // namespace Particle

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
  auto operator()(::Particle::Particle const &p) {
    return p.is_virtual ? 0.0 : p.m;
  }
};
} // namespace Particle
} // namespace Traits

#endif // OBSERVABLES_PARTICLE_HPP
