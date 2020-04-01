/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef ESPRESSO_CORE_PARTICLE_LIST_HPP
#define ESPRESSO_CORE_PARTICLE_LIST_HPP

#include "Particle.hpp"

#include <utils/Span.hpp>
#include <utils/memory.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

/**
 * @brief List of particles.
 */
class ParticleList {
  using storage_type = std::vector<Particle>;

public:
  using value_type = Particle;
  using iterator = Particle *;
  using const_iterator = const Particle *;
  using pointer = Particle *;
  using reference = Particle &;

private:
  std::vector<Particle> m_storage;

  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &m_storage;
  }

public:
  Particle *data() { return m_storage.data(); }
  const Particle *data() const { return m_storage.data(); }

  Particle *begin() { return m_storage.data(); }
  Particle *end() { return m_storage.data() + size(); }
  const Particle *begin() const { return m_storage.data(); }
  const Particle *end() const { return m_storage.data() + size(); }

  /**
   * @brief Resize storage.
   *
   * Newly added Particles are default-initialized.
   *
   *     @param new_size Size to resize to.
   */
  void resize(size_t new_size) { m_storage.resize(new_size); }

  /**
   * @brief Resize the List to zero.
   */
  void clear() { m_storage.clear(); }

  /**
   * @brief Number of entries.
   */
  size_t size() const { return m_storage.size(); }
  bool empty() const { return m_storage.empty(); }

  /**
   * @brief Emplace a particle in the list.
   *
   * @param args Constructor arguments.
   * @return Reference to the added particle.
   */
  template <class... Args> Particle &emplace(Args &&... args) {
    m_storage.emplace_back(std::forward<Args>(args)...);

    return m_storage.back();
  }

  /**
   * @brief Remove element from the list.
   *
   * @param it Iterator pointing to the element to remove.
   * @return An iterator past the element that was removed.
   */
  iterator erase(iterator it) {
    *it = std::move(m_storage.back());

    m_storage.pop_back();

    return it;
  }
};

#endif
