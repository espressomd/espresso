#ifndef OBSERVABLES_PROPERTIES_HPP
#define OBSERVABLES_PROPERTIES_HPP

namespace GenObs {
template <class DoF> struct traits;

struct Position {
  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.position(p);
  }
};

struct Velocity {
  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.velocity(p);
  }
};

struct Mass {
  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.mass(p);
  }
};

struct Charge {
  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.charge(p);
  }
};

struct Force {
  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.force(p);
  }
};
} // namespace GenObs

#endif // OBSERVABLES_PROPERTIES_HPP
