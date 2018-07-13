#ifndef CORE_CONSTRAINTS_CONSTRAINTS_HPP
#define CORE_CONSTRAINTS_CONSTRAINTS_HPP

#include "ObjectRegistry.hpp"
#include "statistics.hpp"

#include <memory>
#include <vector>

namespace Constraints {
template <class ParticleRange, class Constraint>
class Constraints
    : public ObjectRegistry<std::vector<std::shared_ptr<Constraint>>> {

  void reset_foces() const {
    for (auto const &c : *this) {
      c->reset_force();
    }
  }

public:
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
