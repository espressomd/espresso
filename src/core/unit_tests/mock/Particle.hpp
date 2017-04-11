#ifndef CORE_UNIT_TESTS_MOCK_PARTICLE_HPP
#define CORE_UNIT_TESTS_MOCK_PARTICLE_HPP

namespace Testing {
class Particle {
  unsigned m_id;

public:
  explicit Particle(int id) : m_id(id) {}

  int identity() const { return m_id; }
};
}

#endif
