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

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <cstdlib>

/** Attach hand-made particles to a standalone store so force/torque
 *  accessors work in unit tests.
 *
 *  The ParticleStore allocates Kokkos Views, which require an initialized
 *  runtime. Tests that use this fixture but do not otherwise bring up a full
 *  ESPResSo system need not add any boilerplate: the fixture initializes
 *  Kokkos lazily (once) and registers finalize() at program exit. */
struct ParticleStoreTestFixture {
  ParticleStore store{};
  int next_row = 0;

  explicit ParticleStoreTestFixture(std::size_t capacity = 8u) {
    ensure_kokkos_initialized();
    store.begin_rebuild(capacity, 0u);
    store.finish_rebuild();
  }
  /** Attach an existing (freshly default-constructed) view to the next store
   *  row: seeds the row to the new-particle defaults and binds @p p to it
   *  (assign_row's not-preserve branch). Any field the caller wants must be
   *  set AFTER this call, through the now-attached view. */
  void attach(Particle &p) { store.assign_row(p, next_row++); }
  /** Reserve and default-seed the next store row and return a view bound to it
   *  (the store owns the data; the returned view reads/writes it). */
  Particle make() {
    auto const row = next_row++;
    store.seed_default_row(row);
    return store.make_view(row);
  }

private:
  static void ensure_kokkos_initialized() {
    if (not Kokkos::is_initialized() and not Kokkos::is_finalized()) {
      Kokkos::initialize();
      std::atexit([]() {
        if (Kokkos::is_initialized() and not Kokkos::is_finalized()) {
          Kokkos::finalize();
        }
      });
    }
  }
};
