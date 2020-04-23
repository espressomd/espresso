#ifndef TESTS_MOCK_HPP
#define TESTS_MOCK_HPP

namespace Testing {

struct Particle {
  double position = 1.;
  double m_vel = 2.;
  double m_force = 2.1;
  double m_dipole_moment = 1.2;
  double mass = 3.;
  double charge = 0.0;
  auto const &force() const { return m_force; }
  auto const &velocity() const { return m_vel; }
  auto const &dipole_moment() const { return m_dipole_moment; }
};
} // namespace Testing

namespace GenObs {

template <> struct traits<Testing::Particle> {
  using Particle = Testing::Particle;

  double position(Particle const &p) const { return p.position; }
  double velocity(Particle const &p) const { return p.velocity(); }
  double mass(Particle const &p) const { return p.mass; }
  double charge(Particle const &p) const { return p.charge; }
  double force(Particle const &p) const { return p.force(); }
  double dipole_moment(Particle const &p) const { return p.dipole_moment(); }
};
} // namespace GenObs

#endif
