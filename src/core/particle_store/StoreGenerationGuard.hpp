/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "particle_store/ParticleStore.hpp"

#include <cstdint>

/**
 * @brief Reusable debug guard for row-based caches held across steps.
 *
 * Any structure that records @ref ParticleStore ROW indices (rather than
 * @c Particle* pointers) and reuses them across integration steps -- the Verlet
 * list (@c CellStructure::m_verlet_list) and the collision-detection /
 * bond-handler held references -- must stamp the store
 * @ref ParticleStore::generation() when it records the rows and re-check it
 * before resolving them back to views. A rebuild renumbers (or drops) rows and
 * bumps the generation; a cache tagged with an older generation is STALE and a
 * resolved row would alias the wrong particle. The identity gate is blind to
 * this class of bug (positions may coincidentally match after a permuting
 * sort), so the check lives here as a debug canary.
 *
 * Two overloads are provided so a caller can stamp either the raw generation
 * counter alone (when a single store is implied by context) or additionally the
 * store identity (address) -- the latter guards against a fresh store object
 * reusing the same address with a coincidentally-equal generation counter, the
 * same failure mode the columnar ghost cache used to guard against before the
 * ghost layer moved to @ref GhostComm::HaloPlan.
 *
 * The check compiles to nothing in release builds (@c NDEBUG); it is a debug /
 * @c ESPRESSO_ADDITIONAL_CHECKS-time invariant, never a production branch.
 */
namespace ParticleStoreGuard {

/** @brief Assert that @p recorded_generation still matches @p store's current
 *  @ref ParticleStore::generation(). No-op in release builds. */
inline void assert_generation(ParticleStore const &store,
                              std::uint64_t const recorded_generation) {
  assert(recorded_generation == store.generation() and
         "stale store-row cache: recorded generation does not match the "
         "current ParticleStore generation (a rebuild renumbered the rows "
         "without invalidating this cache)");
  static_cast<void>(store);
  static_cast<void>(recorded_generation);
}

/** @brief Assert that both the recorded generation AND the recorded store
 *  identity (address) still match @p store. No-op in release builds. */
inline void assert_generation(ParticleStore const &store,
                              std::uint64_t const recorded_generation,
                              void const *const recorded_store) {
  assert(recorded_store == &store and
         "stale store-row cache: recorded store identity does not match the "
         "current ParticleStore (a different store object is in use)");
  assert_generation(store, recorded_generation);
  static_cast<void>(recorded_store);
}

} // namespace ParticleStoreGuard
