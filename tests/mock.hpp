#ifndef TESTS_MOCK_HPP
#define TESTS_MOCK_HPP

namespace Testing {
struct Particle {
  double position = 1.;
  double m_vel = 2.;
  double m_force = 2.1;
  double mass = 3.;
  double charge = 0.0;
  auto const &force() const { return m_force; }
  auto const &velocity() const { return m_vel; }
};
} // namespace Testing

namespace Observables {
template <> struct traits<Testing::Particle> {
  using Particle = Testing::Particle;

  double position(Particle const &p) const { return p.position; }
  double velocity(Particle const &p) const { return p.velocity(); }
  double mass(Particle const &p) const { return p.mass; }
  double charge(Particle const &p) const { return p.charge; }
  double force(Particle const &p) const { return p.force(); }
};
} // namespace Observables

#endif
