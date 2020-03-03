#ifndef TRAITS_HPP
#define TRAITS_HPP

namespace Traits {

struct Charge {
  template<class Particle>
  auto operator()(Particle const& p) {
    return p.q;
  }
};

struct Position {
  template<class Particle>
  auto operator()(Particle const& p) {
    return p.x;
  }
};

struct Velocity {
  template<class Particle>
  auto operator()(Particle const& p) {
    return p.v;
  }
};

struct Mass {
  template<class Particle>
  auto operator()(Particle const& p) {
    return p.m;
  }
};
}

#endif
