#ifndef CORE_CONSTRAINTS_CONSTRAINTS_HPP
#define CORE_CONSTRAINTS_CONSTRAINTS_HPP

#include "statistics.hpp"

#include <memory>
#include <vector>

void on_constraint_change();

namespace Constraints {
template <class ParticleRange, class Constraint> class Constraints {
  using container_type = std::vector<std::shared_ptr<Constraint>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  void reset_foces() const {
    for (auto const &c : *this) {
      c->reset_force();
    }
  }

  container_type m_constraints;

public:
  void add(std::shared_ptr<Constraint> const &c) {
    if (not c->fits_in_box(Vector3d{box_l})) {
      throw std::runtime_error("Constraint not compatible with box size.");
    }

    m_constraints.emplace_back(c);
    on_constraint_change();
  }
  void remove(std::shared_ptr<Constraint> const &c) {
    m_constraints.erase(
        std::remove(m_constraints.begin(), m_constraints.end(), c),
        m_constraints.end());
    on_constraint_change();
  }

  iterator begin() { return m_constraints.begin(); }
  iterator end() { return m_constraints.end(); }
  const_iterator begin() const { return m_constraints.begin(); }
  const_iterator end() const { return m_constraints.end(); }

  void add_forces(ParticleRange &particles) const {
    reset_foces();

    for (auto &p : particles) {
      auto const pos = folded_position(p);
      ParticleForce force{};
      for (auto const &c : *this) {
        force += c->force(p, pos);
      }

      p.f += force;
    }
  }

  void add_energy(ParticleRange &particles, Observable_stat &energy) const {
    for (auto &p : particles) {
      auto const pos = folded_position(p);

      for (auto const &c : *this) {
        c->add_energy(p, pos, energy);
      }
    }
  }

  void on_boxl_change() const {
    if (not this->empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active constraints.");
    }
  }
};
}

#endif
