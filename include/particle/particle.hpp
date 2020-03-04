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
#endif // OBSERVABLES_PARTICLE_HPP
