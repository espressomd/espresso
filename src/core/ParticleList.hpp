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

/** List of particles. The particle array is resized using a sophisticated
 *  (we hope) algorithm to avoid unnecessary resizes.
 */
struct ParticleList {
  ParticleList() : part{nullptr}, n{0}, max{0} {}
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;

private:
  /** granularity of the particle buffers in particles */
  static constexpr int INCREMENT = 8;

  /** Number of particles that fit in until a resize is needed */
  int max;

  int realloc(int size) {
    assert(size >= 0);
    int old_max = max;
    Particle *old_start = part;

    if (size < max) {
      if (size == 0)
        /* to be able to free an array again */
        max = 0;
      else
        /* shrink not as fast, just lose half, rounded up */
        max = INCREMENT * (((max + size + 1) / 2 + INCREMENT - 1) / INCREMENT);
    } else
      /* round up */
      max = INCREMENT * ((size + INCREMENT - 1) / INCREMENT);
    if (max != old_max)
      part = Utils::realloc(part, sizeof(Particle) * max);
    /* If there are new particles, default initialize them */
    if (max > old_max)
      std::uninitialized_fill(part + old_max, part + max, Particle());

    return part != old_start;
  }

public:
  /** Current allocation size. */
  auto capacity() const { return max; }
  Particle *data() { return part; }

  Particle *begin() { return data(); }
  Particle *end() { return data() + size(); }

  Utils::Span<Particle> particles() { return {part, static_cast<size_t>(n)}; }
  Utils::Span<const Particle> particles() const {
    return {part, static_cast<size_t>(n)};
  }

  Particle &front() { return assert(not empty()), *begin(); }
  Particle &back() { return assert(not empty()), *(std::prev(end())); }

  /**
   * @brief Resize storage for local particles and ghosts.
   *
   * This version does \em not care for the bond information to be freed if
   * necessary.
   *     @param size the size to provide at least.
   *     @return true iff particle addresses have changed
   */
  int resize(int size) { return realloc(this->n = size); }

  /**
   * @brief Resize the List to zero.
   */
  void clear() { resize(0); }

  /**
   * @brief Number of entries.
   */
  int size() const { return n; }
  bool empty() const { return size() <= 0; }

  /**
   * @brief Add a particle at the end of the list.
   *
   * @param p Particle to add.
   */
  void push_back(Particle &&p) {
    resize(size() + 1);
    part[n - 1] = std::move(p);
  }

  /**
   * @brief Move out the last particle in the list.
   *
   * Reduces the size of the list by one.
   *
   * @return Last particle in the list.
   */
  Particle &&extract_back() { return std::move(part[--n]); }

  /**
   * @brief Move out the i-th particle in the list.
   *
   *  This can change the order of particles
   *  from position i onward (including i).
   *
   *  @param i Position  to extract from.
   */
  Particle extract(int i) {
    using std::swap;

    swap(part[i], part[n - 1]);
    return extract_back();
  }
};

#endif
