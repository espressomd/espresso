#ifndef OBSERVABLES_PROPERTIES_HPP
#define OBSERVABLES_PROPERTIES_HPP

namespace GenObs {
template <class DoF> struct traits;

namespace detail {
template<class T>
struct decay{
  using type = typename std::decay<T>;
};

template<class U>
struct decay<std::reference_wrapper<U>> {
  using type = typename std::decay_t<U>;
};

template<class T>
using decay_t = typename decay<T>::type;
}

template<class Particle>
using default_traits = traits<detail::decay_t<Particle>>;



struct Position {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.position(p);
  }
};

struct Velocity {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.velocity(p);
  }
};

struct Mass {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.mass(p);
  }
};

struct Charge {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.charge(p);
  }
};

struct Force {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.force(p);
  }
};

struct DipoleMoment {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.dipole_moment(p);
  }
};
} // namespace GenObs

#endif // OBSERVABLES_PROPERTIES_HPP
